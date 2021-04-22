"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.
"""


import itertools

from snakemake.utils import min_version

min_version('6.1.1')

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------

configfile: 'config.yaml'

#----------------------------------------------------------------------------
# helper functions
#----------------------------------------------------------------------------

def genome_fasta(wc):
    """Get genome FASTA (trimmed or untrimmed)."""
    if config['genome_trim3_polyA']:
        return rules.trim3_polyA.output.fasta
    else:
        return rules.get_genome_fasta.output.fasta

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

rule all:
    input:
        'results/pileup/merged.csv',
        expand('results/consensus_to_genbank_alignments/' +
               "{aligner}/{genome}/{accession}.fa",
               aligner=config['aligners'],
               genome=config['genomes'],
               accession=config['accessions']),

rule get_genome_fasta:
   """Download reference genome fasta."""
   output: fasta="results/genomes/untrimmed/{genome}.fa"
   params: ftp=lambda wildcards: config['genomes'][wildcards.genome]['fasta']
   conda: 'environment.yml'
   shell:
        """
        wget -O - {params.ftp} | gunzip -c > {output}
        python scripts/strip_fasta_head_to_id.py --fasta {output.fasta}
        """

rule trim3_polyA:
    """Trim 3' polyA nucleotides from FASTA."""
    input: fasta=rules.get_genome_fasta.output.fasta
    output: fasta="results/genomes/trim3_polyA/{genome}.fa"
    conda: 'environment.yml'
    script:
        "scripts/trim3_polyA.py"

rule download_accession:
    """Download SRA accession to gzipped FASTQ."""
    output:
        fastq=temp("results/sra_downloads/{accession}.fastq"),
        fastq_gz="results/sra_downloads/{accession}.fastq.gz"
    threads: config['max_cpus']
    params: tempdir=os.path.join(config['scratch_dir'], 'fasterq-dump-temp')
    conda: 'environment.yml'
    shell:
        """
        fasterq-dump \
            {wildcards.accession} \
            --skip-technical \
            --split-spot \
            --outfile {output.fastq} \
            --threads {threads} \
            --force \
            --temp {params.tempdir}
        gzip --keep {output.fastq}
        """

rule bbmap_genome:
    """Build ``bbmap`` reference genome."""
    input: fasta=genome_fasta
    output: path=directory("results/genomes/bbmap_{genome}")
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "bbmap.sh ref={input.fasta} path={output.path} threads={threads}"

rule bwa_mem2_genome:
    """Build ``bwa-mem2`` reference genome."""
    input: fasta=genome_fasta
    output: prefix=directory("results/genomes/bwa-mem2_{genome}/")
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        mkdir -p {output.prefix}
        bwa-mem2 index -p {output.prefix}/index {input.fasta}
        """

rule align_bbmap:
    """Align using ``bbmap``."""
    input:
        fastq=rules.download_accession.output.fastq_gz,
        path=rules.bbmap_genome.output.path,
        ref=genome_fasta
    output:
        sam=temp("results/alignments/bbmap/{genome}/{accession}.sam"),
        bamscript=temp("results/alignments/bbmap/{genome}/{accession}_bamscript.sam"),
        bam="results/alignments/bbmap/{genome}/{accession}_sorted.bam",
    conda: 'environment.yml'
    threads: config['max_cpus']
    shell:
        """
        bbmap.sh \
            in={input.fastq} \
            ref={input.ref} \
            path={input.path} \
            minid=0.8 \
            maxlen=500 \
            threads={threads} \
            outm={output.sam} \
            idtag=t \
            mdtag=t \
            nmtag=t \
            ignorebadquality=t \
            nullifybrokenquality=t \
            ignorejunk=t \
            overwrite=t \
            bamscript={output.bamscript}
        source {output.bamscript}
        """

rule align_bwa_mem2:
    """Align using ``bwa-mem2``."""
    input:
        fastq=rules.download_accession.output.fastq_gz,
        prefix=rules.bwa_mem2_genome.output.prefix,
    output:
        sam=temp("results/alignments/bwa-mem2/{genome}/{accession}.sam"),
        unsorted_bam=temp("results/alignments/bwa-mem2/{genome}/{accession}.bam"),
        bam="results/alignments/bwa-mem2/{genome}/{accession}_sorted.bam",
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        # https://www.biostars.org/p/395057/
        """
        bwa-mem2 mem \
            -t {threads} \
            {input.prefix}/index \
            {input.fastq} > {output.sam}
        samtools view -b -F 4 -o {output.unsorted_bam} {output.sam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        """

rule index_bam:
    """Create BAI file for BAMs."""
    input: bam="{bampath}_sorted.bam"
    output: bai="{bampath}_sorted.bam.bai"
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "samtools index -b -m {threads} {input.bam} {output.bai}"

rule bam_pileup:
    """Make BAM pileup CSVs with mutations."""
    output:
        pileup_csv="results/pileup/{aligner}/{genome}/{accession}.csv",
    input:
        bam=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        }[wc.aligner] + '.bai',
        ref_fasta=lambda wc: {'bbmap': rules.bbmap_genome.input.fasta,
                              'bwa-mem2': rules.bwa_mem2_genome.input.fasta,
                              }[wc.aligner],
    params:
        sample_description=lambda wc: (config['accessions'][wc.accession]
                                       ['sample_description']),
        ref=lambda wc: config['genomes'][wc.genome]['name']
    conda: 'environment.yml'
    shell:
        """
        python scripts/bam_pileup.py \
            --bam {input.bam} \
            --bai {input.bai} \
            --ref {params.ref} \
            --ref_fasta {input.ref_fasta} \
            --pileup_csv {output.pileup_csv} \
            --add_cols aligner {wildcards.aligner} \
            --add_cols genome {wildcards.genome} \
            --add_cols accession {wildcards.accession} \
            --add_cols sample_description "{params.sample_description}"
        """

rule consensus_from_pileup:
    """Make consensus sequence from BAM pileup."""
    input:
        pileup=rules.bam_pileup.output.pileup_csv
    output:
        consensus="results/consensus/{aligner}/{genome}/{accession}.fa"
    params:
        fasta_header = "{aligner}_{genome}_{accession}",
        min_coverage=config['consensus_min_coverage'],
        min_frac=config['consensus_min_frac']
    conda: 'environment.yml'
    shell:
        """
        python scripts/consensus_from_pileup.py \
            --pileup {input.pileup} \
            --consensus {output.consensus} \
            --fasta_header {params.fasta_header} \
            --min_coverage {params.min_coverage} \
            --min_frac {params.min_frac}
        """

rule get_genbank_fasta:
    """Get fasta from Genbank."""
    output: fasta="results/genbank/{genbank}.fa"
    conda: 'environment.yml'
    shell:
        """
        efetch \
            -format fasta \
            -db nuccore \
            -id {wildcards.genbank} \
            > {output.fasta}
        """

rule align_consensus_to_genbank:
    """Align pileup consensus to its Genbank."""
    input:
        consensus=rules.consensus_from_pileup.output.consensus,
        genbank=lambda wc: ('results/genbank/' +
                            config['accessions'][wc.accession]['genbank'] +
                            '.fa')
    output:
        concat_fasta=temp('results/consensus_to_genbank_alignments/' +
                          "{aligner}/{genome}/_{accession}_to_align.fa"),
        alignment=('results/consensus_to_genbank_alignments/' +
                   "{aligner}/{genome}/{accession}.fa")
    conda: 'environment.yml'
    shell:
        # insert newline between FASTA files when concatenating:
        # https://stackoverflow.com/a/23549826
        """
        cat {input.consensus} <(echo) {input.genbank} > {output.concat_fasta}
        mafft {output.concat_fasta} > {output.alignment}
        """

rule analyze_consensus:
    """Analyze consensus sequences from pileup versus Genbank."""
    input:
        consensus_seqs=expand(rules.consensus_from_pileup.output.consensus,
                              aligner=config['aligners'],
                              genome=config['genomes'],
                              accession=config['accessions'],
                              ),
        genbanks=[f"results/genbank/{config['accessions'][accession]['genbank']}.fa"
                  for _, _, accession
                  in itertools.product(config['aligners'],
                                       config['genomes'],
                                       config['accessions'])
                  ]
    params:
        descriptors=[{'aligner': aligner, 'genome': genome, 'accession': accession}
                     for aligner, genome, accession
                     in itertools.product(config['aligners'],
                                          config['genomes'],
                                          config['accessions'])
                     ]
    log:
        notebook='results/logs/notebooks/analyze_consensus.ipynb'
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_consensus.py.ipynb'

rule merge_pileup_csv:
    output:
        merged_csv='results/pileup/merged.csv'
    input:
        expand("results/pileup/{aligner}/{genome}/{accession}.csv",
               aligner=config['aligners'],
               accession=config['accessions'],
               genome=config['genomes'],
               ),
    conda: 'environment.yml'
    script:
        'scripts/merge_csvs.py'

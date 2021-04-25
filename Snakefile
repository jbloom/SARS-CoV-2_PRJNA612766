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

samples = config['samples']  # read samples from config

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
        #'results/consensus_to_genbank_alignments/stats.csv',
        #'results/consensus_to_genbank_alignments/chart.html',
        expand("results/pileup/{sample}/interactive_pileup_chart.html",
               sample=samples),
        expand("results/pileup/{sample}/diffs_from_ref.csv",
               sample=samples)

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
        fastqs=lambda wc: [f"results/sra_downloads/{accession}.fastq.gz"
                           for accession in samples[wc.sample]['accessions']],
        path=rules.bbmap_genome.output.path,
        ref=genome_fasta
    output:
        concat_fastq=temp("results/alignments/bbmap/{genome}/_{sample}" +
                          '_concat.fastq.gz'),
        sam=temp("results/alignments/bbmap/{genome}/{sample}.sam"),
        bamscript=temp("results/alignments/bbmap/{genome}/{sample}" +
                       '_bamscript.sam'),
        bam="results/alignments/bbmap/{genome}/{sample}_sorted.bam",
    conda: 'environment.yml'
    threads: config['max_cpus']
    shell:
        """
        echo "concatenating {input.fastqs}"
        cat {input.fastqs} > {output.concat_fastq}
        echo "done concatenating"
        ls -lh {output.concat_fastq}
        echo "mapping {input.fastqs}"
        bbmap.sh \
            in={output.concat_fastq} \
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
        echo "done mapping"
        ls -lh {output.concat_fastq}
        """

rule align_bwa_mem2:
    """Align using ``bwa-mem2``."""
    input:
        fastqs=lambda wc: [f"results/sra_downloads/{accession}.fastq.gz"
                           for accession in samples[wc.sample]['accessions']],
        prefix=rules.bwa_mem2_genome.output.prefix,
    output:
        concat_fastq=temp("results/alignments/bwa-mem2/{genome}/_{sample}" +
                          '_concat.fastq.gz'),
        sam=temp("results/alignments/bwa-mem2/{genome}/{sample}.sam"),
        unsorted_bam=temp("results/alignments/bwa-mem2/{genome}/{sample}.bam"),
        bam="results/alignments/bwa-mem2/{genome}/{sample}_sorted.bam",
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        # https://www.biostars.org/p/395057/
        """
        cat {input.fastqs} > {output.concat_fastq}
        bwa-mem2 mem \
            -t {threads} \
            {input.prefix}/index \
            {output.concat_fastq} > {output.sam}
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
        pileup_csv="results/pileup/{sample}/pileup_{genome}_{aligner}.csv",
    input:
        bam=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        }[wc.aligner] + '.bai',
        ref_fasta=genome_fasta
    params:
        ref=lambda wc: config['genomes'][wc.genome]['name'],
        minq=config['minq'],
    conda: 'environment.yml'
    shell:
        """
        python scripts/bam_pileup.py \
            --bam {input.bam} \
            --bai {input.bai} \
            --ref {params.ref} \
            --ref_fasta {input.ref_fasta} \
            --minq {params.minq} \
            --pileup_csv {output.pileup_csv} \
        """

rule consensus_from_pileup:
    """Make consensus sequence from BAM pileup."""
    input:
        pileup=rules.bam_pileup.output.pileup_csv
    output:
        consensus="results/consensus/{sample}/consensus_{genome}_{aligner}.fa"
    params:
        fasta_header = "{sample}_{genome}_{aligner}",
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
                            samples[wc.sample]['genbank'] +
                            '.fa')
    output:
        concat_fasta=temp('results/consensus_to_genbank_alignments/' +
                          "{sample}/_{genome}_{aligner}_to_align.fa"),
        alignment=('results/consensus_to_genbank_alignments/' +
                   "{sample}/alignment_{genome}_{aligner}.fa")
    conda: 'environment.yml'
    shell:
        # insert newline between FASTA files when concatenating:
        # https://stackoverflow.com/a/23549826
        """
        cat {input.consensus} <(echo) {input.genbank} > {output.concat_fasta}
        mafft {output.concat_fasta} > {output.alignment}
        """

rule analyze_consensus_vs_genbank:
    """Analyze consensus sequences from pileup versus Genbank."""
    output:
        csv='results/consensus_to_genbank_alignments/stats.csv',
        chart='results/consensus_to_genbank_alignments/chart.html',
    input:
        alignments=expand(rules.align_consensus_to_genbank.output.alignment,
                          aligner=config['aligners'],
                          genome=config['genomes'],
                          sample=config['samples'],
                          ),
    params:
        descriptors=[{'aligner': aligner,
                      'genome': genome,
                      'sample': sample}
                     for aligner, genome, sample
                     in itertools.product(config['aligners'],
                                          config['genomes'],
                                          samples)
                     ]
    log:
        notebook='results/logs/notebooks/analyze_consensus_vs_genbank.ipynb'
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_consensus_vs_genbank.py.ipynb'

rule analyze_pileups:
    """Analyze and plot BAM pileups."""
    input:
        pileups=expand(rules.bam_pileup.output.pileup_csv,
                       genome=config['genomes'],
                       aligner=config['aligners'],
                       allow_missing=True)
    output:
        chart="results/pileup/{sample}/interactive_pileup_chart.html",
        diffs_from_ref="results/pileup/{sample}/diffs_from_ref.csv",
    params:
        consensus_min_frac=config['consensus_min_frac'],
        consensus_min_coverage=config['consensus_min_coverage'],
        descriptors=[{'genome': genome,
                      'aligner': aligner}
                     for genome, aligner
                     in itertools.product(config['genomes'],
                                          config['aligners'])
                     ],
        chart_title="{sample}"
    log:
        notebook="results/logs/notebooks/analyze_pileups_{sample}.ipynb"
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_pileups.py.ipynb'

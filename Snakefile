"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.
"""

configfile: 'config.yaml'

rule all:
    input:
        'results/pileup/merged.csv'

rule get_genome_fasta:
   """Download reference genome fasta."""
   output: fasta="results/genomes/untrimmed/{genome}.fa"
   params: ftp=lambda wildcards: config['genomes'][wildcards.genome]['fasta']
   conda: 'environment.yml'
   shell:
        # rename FASTA to {genome} wildcard name
        """
        wget -O - {params.ftp} | gunzip -c > {output}
        if [ $(grep '^>' {output} | wc -l) != 1 ]; then
            echo "more than one entry in {output}" 1>&2
            exit 1
        fi
        sed -i 's/^>.*$/>{wildcards.genome}/g' {output}
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
    input:
        fasta=(rules.trim3_polyA.output.fasta
               if config['genome_trim3_polyA'] else
               rules.get_genome_fasta.output.fasta)
    output: path=directory("results/genomes/bbmap_{genome}")
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "bbmap.sh ref={input.fasta} path={output.path} threads={threads}"

rule bwa_mem2_genome:
    """Build ``bwa-mem2`` reference genome."""
    input:
        fasta=(rules.trim3_polyA.output.fasta
               if config['genome_trim3_polyA'] else
               rules.get_genome_fasta.output.fasta)
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
        ref=rules.bbmap_genome.input.fasta,
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
                                       ['sample_description'])
    conda: 'environment.yml'
    shell:
        """
        python scripts/bam_pileup.py \
            --bam {input.bam} \
            --bai {input.bai} \
            --ref {wildcards.genome} \
            --ref_fasta {input.ref_fasta} \
            --pileup_csv {output.pileup_csv} \
            --add_cols aligner {wildcards.aligner} \
            --add_cols genome {wildcards.genome} \
            --add_cols accession {wildcards.accession} \
            --add_cols sample_description "{params.sample_description}"
        """

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

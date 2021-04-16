"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.
"""

configfile: 'config.yaml'

rule all:
    input:
        expand("results/genomes/{aligner}_{genome}",
               aligner=config['aligners'],
               genome=config['genomes'],
               ),
        expand("results/sra_downloads/{accession}.fastq.gz",
               accession=config['accessions'])

rule get_genome_fasta:
   """Download reference genome fasta."""
   output: fasta="results/genomes/untrimmed/{genome}.fa"
   params: ftp=lambda wildcards: config['genomes'][wildcards.genome]['fasta']
   conda: 'environment.yml'
   shell: "wget -O - {params.ftp} | gunzip -c > {output}"

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
        ref=(rules.trim3_polyA.output.fasta
             if config['genome_trim3_polyA'] else
             rules.get_genome_fasta.output.fasta)
    output: path=directory("results/genomes/bbmap_{genome}")
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "bbmap.sh ref={input.ref} path={output.path} threads={threads}"

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

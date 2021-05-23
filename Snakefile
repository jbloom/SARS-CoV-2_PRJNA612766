"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.
"""


import itertools
import os

from snakemake.utils import min_version

min_version('6.3.0')

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------

configfile: 'config.yaml'

# get all samples
samples = {key: val for key, val in
           list(config['samples_fasterq_dump'].items()) +
           list(config['samples_wget'].items())
           }

#----------------------------------------------------------------------------
# helper functions
#----------------------------------------------------------------------------

def use_wget(wc):
    """For a given accession, do we use `wget`?"""
    for sample_d in config['samples_wget'].values():
        if wc.accession in sample_d['accessions']:
            return 'use_wget'
    return 'no_wget'

def comparator_fastas(wc):
    """Get FASTA files for all comparator genomes."""
    comparators = []
    for key, d in config['comparator_genomes'].items():
        if 'genbank' in d:
            comparators.append(f"results/genbank/{d['genbank']}.fa")
        elif 'gisaid' in d:
            comparators.append(f"results/gisaid/{d['gisaid']}.fa")
        else:
            raise ValueError(f"comparator {key} lacks genbank and gisaid")
    return comparators

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

rule all:
    input:
        expand("results/pileup/{sample}/interactive_pileup.html",
               sample=samples),
        'results/pileup/frac_coverage.csv',
        'results/pileup/frac_coverage.html',
        'results/pileup/diffs_from_ref.csv',
        'results/pileup/diffs_from_ref.html',
        'results/consensus/consensus_seqs.csv',
        'results/comparator_annotated_gisaid_muts/muts.csv.gz',

rule get_ref_genome_fasta:
    """Download reference genome fasta."""
    output: fasta="results/ref_genome/ref_genome.fa"
    params:
        url=config['ref_genome']['fasta'],
        name=config['ref_genome']['name'],
        add_mutations=config['ref_genome']['add_mutations']
    conda: 'environment.yml'
    script:
        'scripts/get_ref_genome_fasta.py'

rule get_genome_gff:
    """Download reference genome GFF."""
    output: gff="results/ref_genomes/ref_genome.gff"
    params: url=config['ref_genome']['gff']
    conda: 'environment.yml'
    shell:
        "wget -O - {params.url} | gunzip -c > {output}"

rule download_sra:
    """Download SRA accession to gzipped FASTQ, concat when multiple FASTQs.

    Code is complicated because if sample is in `samples_wget` then we
    use `wget` to get it, and otherwise `fasterq-dump`.
    """
    output:
        fastq_dir=temp(directory("results/sra_downloads/{accession}/")),
        fastq_gz=protected("results/sra_downloads/{accession}.fastq.gz"),
        temp_dir=temp(directory(os.path.join(config['scratch_dir'],
                                             "fasterq-dump/{accession}"))),
        sra_file=protected("results/sra_downloads/{accession}.sra"),
    params:
        use_wget=use_wget,
        wget_paths=['https://storage.googleapis.com/nih-sequence-read-archive/run',
                    'https://sra-pub-sars-cov2.s3.amazonaws.com/run',
                    'https://sra-pub-run-odp.s3.amazonaws.com/sra']
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        if [[ "{params.use_wget}" == "use_wget" ]]
        then
            echo "using wget for {wildcards.accession}"
            wget {params.wget_paths[0]}/{wildcards.accession}/{wildcards.accession} \
                -O {output.sra_file} || \
            wget {params.wget_paths[1]}/{wildcards.accession}/{wildcards.accession} -O \
                {output.sra_file} || \
            wget {params.wget_paths[2]}/{wildcards.accession}/{wildcards.accession} -O \
                {output.sra_file}
            acc="{output.sra_file}"
        else
            echo "not using wget for {wildcards.accession}"
            touch {output.sra_file}
            acc="{wildcards.accession}"
        fi
        fasterq-dump \
            $acc \
            --skip-technical \
            --split-spot \
            --outdir {output.fastq_dir} \
            --threads {threads} \
            --force \
            --temp {output.temp_dir}
        pigz -c -p {threads} {output.fastq_dir}/*.fastq > {output.fastq_gz}
        """

rule preprocess_fastq:
    """Pre-process the FASTQ files by trimming adaptors etc with ``fastp``."""
    input:
        fastq_gz=rules.download_sra.output.fastq_gz
    output:
        fastq_gz=temp("results/preprocessed_fastqs/{accession}.fastq.gz"),
        html="results/preprocessed_fastqs/{accession}.html",
        json="results/preprocessed_fastqs/{accession}.json",
    params:
        minq=config['minq'],
        min_read_length=config['min_read_length'],
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        # filter if >40% of read has quality < minq, if read shorter
        # than minimum read-length, and polyG / polyX tails.
        """
        fastp \
            -i {input.fastq_gz} \
            -q {params.minq} \
            -u 40 \
            -l 20 {params.min_read_length} \
            --trim_poly_g \
            --trim_poly_x \
            -o {output.fastq_gz} \
            --html {output.html} \
            --json {output.json}
        """

rule bwa_mem2_genome:
    """Build ``bwa-mem2`` reference genome."""
    input: fasta=rules.get_ref_genome_fasta.output.fasta
    output: prefix=directory("results/genomes/bwa-mem2/")
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        mkdir -p {output.prefix}
        bwa-mem2 index -p {output.prefix}/index {input.fasta}
        """

rule minimap2_genome:
    """Build ``minimap2`` reference genome."""
    input: fasta=rules.get_ref_genome_fasta.output.fasta
    output: mmi="results/genomes/minimap2.mmi"
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "minimap2 -t {threads} -d {output.mmi} {input.fasta}"

rule align_bwa_mem2:
    """Align using ``bwa-mem2``."""
    input:
        fastqs=lambda wc: expand(rules.preprocess_fastq.output.fastq_gz,
                                 accession=samples[wc.sample]['accessions']),
        prefix=rules.bwa_mem2_genome.output.prefix,
    output:
        concat_fastq=temp("results/alignments/bwa-mem2/_{sample}" +
                          '_concat.fastq.gz'),
        sam=temp("results/alignments/bwa-mem2/{sample}.sam"),
        unsorted_bam=temp("results/alignments/bwa-mem2/{sample}.bam"),
        bam="results/alignments/bwa-mem2/{sample}_sorted.bam",
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

rule align_minimap2:
    """Align using ``minimap2``."""
    input:
        fastqs=lambda wc: expand(rules.preprocess_fastq.output.fastq_gz,
                                 accession=samples[wc.sample]['accessions']),
        mmi=rules.minimap2_genome.output.mmi,
    output:
        concat_fastq=temp("results/alignments/minimap2/_{sample}" +
                          '_concat.fastq.gz'),
        sam=temp("results/alignments/minimap2/{sample}.sam"),
        unsorted_bam=temp("results/alignments/minimap2/{sample}.bam"),
        bam="results/alignments/minimap2/{sample}_sorted.bam",
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        cat {input.fastqs} > {output.concat_fastq}
        minimap2 -a {input.mmi} {input.fastqs} > {output.sam}
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
        pileup_csv="results/pileup/{sample}/pileup_{aligner}.csv",
    input:
        bam=lambda wc: {'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner] + '.bai',
        ref_fasta=rules.get_ref_genome_fasta.output.fasta
    params:
        ref=config['ref_genome']['name'],
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
        consensus="results/consensus/{sample}/consensus_{aligner}.fa"
    params:
        fasta_header = "{sample}_{aligner}",
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

rule get_gisaid_fasta:
    """Get FASTA from GISAID downloads."""
    output: fasta="results/gisaid/{gisaid}.fa"
    params: gisaid_dirs=config['gisaid_dirs']
    conda: 'environment.yml'
    script:
        'scripts/get_gisaid_fasta.py'

rule genome_comparator_alignment:
    """Align genome to comparators."""
    input:
        genome=rules.get_ref_genome_fasta.output.fasta,
        comparators=comparator_fastas
    output:
        concat_fasta=temp('results/genome_to_comparator/to_align.fa'),
        alignment="results/genome_to_comparator/alignment.fa"
    conda: 'environment.yml'
    shell:
        # insert newline between FASTA files when concatenating:
        # https://stackoverflow.com/a/25030513/4191652
        """
        awk 1 {input} > {output.concat_fasta}
        mafft {output.concat_fasta} > {output.alignment}
        """

rule genome_comparator_map:
    """Map sites in viral genome to comparator identities."""
    input:
        alignment=rules.genome_comparator_alignment.output.alignment
    output:
        site_map="results/genome_to_comparator/site_identity_map.csv"
    params:
        comparators=list(config['comparator_genomes'])
    conda: 'environment.yml'
    script:
        'scripts/genome_comparator_map.py'

rule annotate_gisaid_muts_by_comparators:
    """Annotate GISAID mutations by comparator genomes."""
    input:
        comparator_map=rules.genome_comparator_map.output.site_map,
        genome_fasta=rules.get_ref_genome_fasta.output.fasta,
        gisaid_metadata='data/gisaid_mutations/metadata.tsv.gz',
        gisaid_muts='data/gisaid_mutations/mut_summary.tsv.gz',
    output:
        annotated_muts='results/comparator_annotated_gisaid_muts/muts.csv.gz'
    log:
        notebook='results/logs/notebooks/annotate_gisaid_muts_by_comparators.py.ipynb'
    params:
        add_mutations=config['ref_genome']['add_mutations']        
    conda: 'environment.yml'
    notebook:
        'notebooks/annotate_gisaid_muts_by_comparators.py.ipynb'

rule analyze_pileups:
    """Analyze and plot BAM pileups per sample."""
    input:
        pileups=expand(rules.bam_pileup.output.pileup_csv,
                       aligner=config['aligners'],
                       allow_missing=True)
    output:
        chart="results/pileup/{sample}/interactive_pileup.html",
        diffs_from_ref="results/pileup/{sample}/diffs_from_ref.csv",
        frac_coverage="results/pileup/{sample}/frac_coverage.csv",
    params:
        consensus_min_frac=config['consensus_min_frac'],
        consensus_min_coverage=config['consensus_min_coverage'],
        report_frac_coverage=config['report_frac_coverage'],
        descriptors=[{'aligner': aligner} for aligner in config['aligners']],
        chart_title="{sample}"
    log:
        notebook="results/logs/notebooks/analyze_pileups_{sample}.ipynb"
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_pileups.py.ipynb'

rule aggregate_pileup_analysis:
    """Analyze viral deep sequencing aggregated across samples."""
    input:
        diffs_from_ref=expand(rules.analyze_pileups.output.diffs_from_ref,
                              sample=samples),
        frac_coverage=expand(rules.analyze_pileups.output.frac_coverage,
                             sample=samples),
        comparator_map=rules.genome_comparator_map.output.site_map,
    output:
        frac_coverage_stats='results/pileup/frac_coverage.csv',
        frac_coverage_chart='results/pileup/frac_coverage.html',
        diffs_from_ref_stats='results/pileup/diffs_from_ref.csv',
        diffs_from_ref_chart='results/pileup/diffs_from_ref.html',
    params:
        samples=list(samples),
    conda: 'environment.yml'
    log:
        notebook='results/logs/notebooks/aggregate_pileup_analysis.ipynb'
    notebook:
        'notebooks/aggregate_pileup_analysis.py.ipynb'

rule aggregate_consensus_seqs:
    """Aggregate the consensus sequences from the pileup."""
    input:
        consensus_seqs=expand(rules.consensus_from_pileup.output.consensus,
                              aligner=config['aligners'],
                              sample=samples),
    output:
        csv='results/consensus/consensus_seqs.csv'
    params:
        descriptors=[{'aligner': aligner, 'sample': sample} for
                     aligner, sample in itertools.product(config['aligners'],
                                                          samples)]
    conda: 'environment.yml'
    log: notebook='results/logs/notebooks/aggregate_consensus_seqs.ipynb'
    notebook: 'notebooks/aggregate_consensus_seqs.py.ipynb'

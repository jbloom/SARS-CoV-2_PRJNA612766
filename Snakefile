"""``snakemake`` file that runs analysis.

Written by Jesse Bloom.
"""


import itertools
import os

from snakemake.utils import min_version

min_version('6.1.1')

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------

configfile: 'config.yaml'

# get all samples
samples = {key: val for key, val in
           list(config['samples_fasterq_dump'].items()) +
           list(config['samples_wget'].items())
           }

report: 'report/workflow.rst'

#----------------------------------------------------------------------------
# helper functions
#----------------------------------------------------------------------------

def genome_fasta(wc):
    """Get genome FASTA (trimmed or untrimmed)."""
    if wc.genome in config['genomes'] and config['genome_trim3_polyA']:
        return rules.trim3_polyA.output.fasta
    else:
        return rules.get_genome_fasta.output.fasta

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
        'results/consensus_vs_genbank/stats.csv',
        'results/consensus_vs_genbank/chart.html',
        'results/consensus_vs_genbank/mismatches.csv',
        expand("results/pileup/{sample}/interactive_pileup.html",
               sample=samples),
        'results/pileup/frac_coverage.csv',
        'results/pileup/frac_coverage.html',
        'results/pileup/diffs_from_ref.csv',
        'results/pileup/diffs_from_ref.html',
        'results/sex_chromosome/stats.csv',
        'results/sex_chromosome/chart.html',
        expand("results/consensus/{sample}/consensus_{genome}_{aligner}.fa",
               sample=samples,
               genome=config['genomes'],
               aligner=config['aligners']),
        'results/ivar_variants/aggregated_ivar_variants.csv',
        expand("results/comparator_annotated_gisaid_muts/{genome}.csv.gz",
               genome=config['genomes']),

rule get_genome_fasta:
    """Download reference genome fasta."""
    output: fasta="results/genomes/untrimmed/{genome}.fa"
    params:
        ftp=lambda wildcards: (config['genomes'][wildcards.genome]['fasta']
                               if wildcards.genome in config['genomes'] else
                               config['host_genomes'][wildcards.genome]['fasta'])
    conda: 'environment.yml'
    shell:
        """
        wget -O - {params.ftp} | gunzip -c > {output}
        python scripts/strip_fasta_head_to_id.py --fasta {output.fasta}
        """

rule get_genome_gff:
    """Download reference genome GFF."""
    output: gff="results/genomes/untrimmed/{genome}.gff"
    params: ftp=lambda wc: config['genomes'][wc.genome]['gff']
    conda: 'environment.yml'
    shell:
        "wget -O - {params.ftp} | gunzip -c > {output}"

rule trim3_polyA:
    """Trim 3' polyA nucleotides from FASTA."""
    input: fasta=rules.get_genome_fasta.output.fasta
    output: fasta="results/genomes/trim3_polyA/{genome}.fa"
    conda: 'environment.yml'
    script:
        "scripts/trim3_polyA.py"

rule download_sra:
    """Download SRA accession to gzipped FASTQ, concat when multiple FASTQs.

    Code is complicated because if sample is in `samples_wget` then we
    use `wget` to get it, and otherwise `fasterq-dump`.
    """
    output:
        fastq_dir=temp(directory("results/sra_downloads/{accession}/")),
        fastq_gz="results/sra_downloads/{accession}.fastq.gz",
        temp_dir=temp(directory(os.path.join(config['scratch_dir'],
                                             "fasterq-dump/{accession}"))),
        sra_file=temp(os.path.join(config['scratch_dir'], "{accession}.sra"))
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

rule minimap2_genome:
    """Build ``minimap2`` reference genome."""
    input: fasta=genome_fasta
    output: mmi="results/genomes/minimap2_{genome}.mmi"
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        "minimap2 -t {threads} -d {output.mmi} {input.fasta}"

rule align_bbmap:
    """Align using ``bbmap``."""
    input:
        fastqs=lambda wc: expand(rules.preprocess_fastq.output.fastq_gz,
                                 accession=samples[wc.sample]['accessions']),
        prefix=rules.bwa_mem2_genome.output.prefix,
        path=rules.bbmap_genome.output.path,
        ref=genome_fasta
    output:
        concat_fastq=temp("results/alignments/bbmap/{genome}/_{sample}" +
                          '_concat.fastq.gz'),
        concat_fasta=temp("results/alignments/bbmap/{genome}/_{sample}" +
                          '_concat.fasta.gz'),
        sam=temp("results/alignments/bbmap/{genome}/{sample}.sam"),
        bamscript=temp("results/alignments/bbmap/{genome}/{sample}" +
                       '_bamscript.sam'),
        bam="results/alignments/bbmap/{genome}/{sample}_sorted.bam",
    params:
        # use perfect mode only for host reads
        perfectmode=lambda wc: 'f' if wc.genome in config['genomes'] else 't'
    conda: 'environment.yml'
    threads: config['max_cpus']
    shell:
        # first try to align FASTQ, then FASTA if that fails. Convert like this:
        # https://bioinformaticsworkbook.org/dataWrangling/fastaq-manipulations/converting-fastq-format-to-fasta.html#gsc.tab=0
        """
        echo "concatenating {input.fastqs}"
        cat {input.fastqs} > {output.concat_fastq}
        (
            (echo "mapping {output.concat_fastq}" &&
             bbmap.sh \
                in={output.concat_fastq} \
                ref={input.ref} \
                path={input.path} \
                minid=0.8 \
                perfectmode={params.perfectmode} \
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
                bamscript={output.bamscript} &&
             touch {output.concat_fasta}
             ) ||
            (echo "making {output.concat_fasta}" &&
             zcat {output.concat_fastq} | sed -n '1~4s/^@/>/p;2~4p' - | \
                gzip > {output.concat_fasta} &&
             echo "mapping {output.concat_fasta}" &&
             bbmap.sh \
                in={output.concat_fasta} \
                ref={input.ref} \
                path={input.path} \
                minid=0.8 \
                perfectmode={params.perfectmode} \
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
             )
        )
        echo "running {output.bamscript}"
        source {output.bamscript}
        """

rule align_bwa_mem2:
    """Align using ``bwa-mem2``."""
    input:
        fastqs=lambda wc: expand(rules.preprocess_fastq.output.fastq_gz,
                                 accession=samples[wc.sample]['accessions']),
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

rule align_minimap2:
    """Align using ``minimap2``."""
    input:
        fastqs=lambda wc: expand(rules.preprocess_fastq.output.fastq_gz,
                                 accession=samples[wc.sample]['accessions']),
        mmi=rules.minimap2_genome.output.mmi,
    output:
        concat_fastq=temp("results/alignments/minimap2/{genome}/_{sample}" +
                          '_concat.fastq.gz'),
        sam=temp("results/alignments/minimap2/{genome}/{sample}.sam"),
        unsorted_bam=temp("results/alignments/minimap2/{genome}/{sample}.bam"),
        bam="results/alignments/minimap2/{genome}/{sample}_sorted.bam",
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

rule ivar_variants:
    """Call variants using ``ivar``."""
    input:
        bam=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner] + '.bai',
        ref_fasta=genome_fasta,
        ref_gff=rules.get_genome_gff.output.gff
    output:
        tsv="results/ivar_variants/{sample}/{aligner}_{genome}.tsv"
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.tsv)[0],
        min_threshold=config['consensus_min_frac'],
        min_depth=config['consensus_min_coverage'],
        minq=config['minq'],
    conda: 'environment.yml'
    shell:
        """
        samtools mpileup -aa -A -d 0 -B -Q 0 \
            --reference {input.ref_fasta} {input.bam} |
        ivar variants -p {params.prefix} \
            -q {params.minq} \
            -t {params.min_threshold} \
            -m {params.min_depth} \
            -r {input.ref_fasta} \
            -g {input.ref_gff}
        """

rule aggregate_ivar_variants:
    """Aggregate ``ivar`` analysis of variants."""
    input:
        ref_gffs=expand(rules.get_genome_gff.output.gff, genome=config['genomes']),
        tsvs=expand("results/ivar_variants/{sample}/{aligner}_{genome}.tsv",
                    aligner=config['aligners'],
                    genome=config['genomes'],
                    sample=samples),
    output:
        agg_csv='results/ivar_variants/aggregated_ivar_variants.csv',
    params:
        descriptors=[{'aligner': aligner,
                      'genome': genome,
                      'sample': sample}
                     for aligner, genome, sample
                     in itertools.product(config['aligners'],
                                          config['genomes'],
                                          samples)
                     ],
        genomes=config['genomes'],
    conda: 'environment.yml'
    script:
        'scripts/aggregate_ivar_variants.py'

rule bam_pileup:
    """Make BAM pileup CSVs with mutations."""
    output:
        pileup_csv="results/pileup/{sample}/pileup_{genome}_{aligner}.csv",
    input:
        bam=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
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
        genome=genome_fasta,
        comparators=comparator_fastas
    output:
        concat_fasta=temp('results/genome_to_comparator/' +
                          "{genome}/to_align.fa"),
        alignment="results/genome_to_comparator/{genome}/alignment.fa"
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
        site_map="results/genome_to_comparator/{genome}/site_identity_map.csv"
    params:
        comparators=list(config['comparator_genomes'])
    conda: 'environment.yml'
    script:
        'scripts/genome_comparator_map.py'

rule annotate_gisaid_muts_by_comparators:
    """Annotate GISAID mutations by comparator genomes."""
    input:
        comparator_map=rules.genome_comparator_map.output.site_map,
        genome_fasta=genome_fasta,
        gisaid_metadata='data/gisaid_mutations/metadata.tsv.gz',
        gisaid_muts='data/gisaid_mutations/mut_summary.tsv.gz',
    output:
        annotated_muts="results/comparator_annotated_gisaid_muts/{genome}.csv.gz"
    log:
        notebook="results/logs/notebooks/{genome}_annotate_gisaid_muts_by_comparators.py.ipynb"
    conda: 'environment.yml'
    notebook:
        'notebooks/annotate_gisaid_muts_by_comparators.py.ipynb'

rule align_consensus_to_genbank:
    """Align pileup consensus to its Genbank and other comparators."""
    input:
        consensus=rules.consensus_from_pileup.output.consensus,
        genbank=lambda wc: ('results/genbank/' +
                            samples[wc.sample]['genbank'] +
                            '.fa'),
        comparators=comparator_fastas,
    output:
        concat_fasta=temp('results/consensus_vs_genbank/' +
                          "{sample}/_{genome}_{aligner}_to_align.fa"),
        alignment=('results/consensus_vs_genbank/' +
                   "{sample}/alignment_{genome}_{aligner}.fa")
    conda: 'environment.yml'
    shell:
        # insert newline between FASTA files when concatenating:
        # https://stackoverflow.com/a/25030513/4191652
        """
        awk 1 {input} > {output.concat_fasta}
        mafft {output.concat_fasta} > {output.alignment}
        """

rule analyze_consensus_vs_genbank:
    """Analyze consensus sequences from pileup versus Genbank."""
    output:
        csv=report('results/consensus_vs_genbank/stats.csv',
                   category='Deep sequencing vs Genbank',
                   caption='report/analyze_consensus_vs_genbank_csv.rst',
                   ),
        chart=report('results/consensus_vs_genbank/chart.html',
                     category='Deep sequencing vs Genbank',
                     caption='report/analyze_consensus_vs_genbank_chart.rst',
                     ),
        mismatches=report('results/consensus_vs_genbank/mismatches.csv',
                          category='Deep sequencing vs Genbank',
                          caption='report/analyze_consensus_vs_genbank_mismatches.rst',
                          )
    input:
        alignments=expand(rules.align_consensus_to_genbank.output.alignment,
                          aligner=config['aligners'],
                          genome=config['genomes'],
                          sample=[s for s in samples if 'genbank' in samples[s]],
                          ),
    params:
        descriptors=[{'aligner': aligner,
                      'genome': genome,
                      'sample': sample}
                     for aligner, genome, sample
                     in itertools.product(config['aligners'],
                                          config['genomes'],
                                          [s for s in samples if 'genbank' in samples[s]])
                     ],
        comparators=list(config['comparator_genomes'])
    log:
        notebook='results/logs/notebooks/analyze_consensus_vs_genbank.ipynb'
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_consensus_vs_genbank.py.ipynb'

rule analyze_pileups:
    """Analyze and plot BAM pileups per sample."""
    input:
        pileups=expand(rules.bam_pileup.output.pileup_csv,
                       genome=config['genomes'],
                       aligner=config['aligners'],
                       allow_missing=True)
    output:
        chart=(report("results/pileup/{sample}/interactive_pileup.html",
                      category='Viral deep sequencing analysis',
                      subcategory='Per-sample pileup files',
                      caption='report/analyze_pileups_interactive_pileup.rst',
                      )
               if config['per_sample_pileups_in_report'] else
               "results/pileup/{sample}/interactive_pileup.html"),
        diffs_from_ref="results/pileup/{sample}/diffs_from_ref.csv",
        frac_coverage="results/pileup/{sample}/frac_coverage.csv",
    params:
        consensus_min_frac=config['consensus_min_frac'],
        consensus_min_coverage=config['consensus_min_coverage'],
        report_frac_coverage=config['report_frac_coverage'],
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

rule aggregate_pileup_analysis:
    """Analyze viral deep sequencing aggregated across samples."""
    input:
        diffs_from_ref=expand(rules.analyze_pileups.output.diffs_from_ref,
                              sample=samples),
        frac_coverage=expand(rules.analyze_pileups.output.frac_coverage,
                             sample=samples),
        comparator_map=expand(rules.genome_comparator_map.output.site_map,
                              genome=config['genomes']),
        ivar_variants='results/ivar_variants/aggregated_ivar_variants.csv'
    output:
        frac_coverage_stats=report(
                'results/pileup/frac_coverage.csv',
                caption='report/aggregate_pileup_analysis_frac_coverage_stats.rst',
                category='Viral deep sequencing analysis'),
        frac_coverage_chart=report(
                'results/pileup/frac_coverage.html',
                caption='report/aggregate_pileup_analysis_frac_coverage_chart.rst',
                category='Viral deep sequencing analysis'),
        diffs_from_ref_stats=report(
                'results/pileup/diffs_from_ref.csv',
                caption='report/aggregate_pileup_analysis_diffs_from_ref_stats.rst',
                category='Viral deep sequencing analysis'),
        diffs_from_ref_chart=report(
                'results/pileup/diffs_from_ref.html',
                caption='report/aggregate_pileup_analysis_diffs_from_ref_chart.rst',
                category='Viral deep sequencing analysis'),
    params:
        samples=list(samples),
        genomes=list(config['genomes']),
    conda: 'environment.yml'
    log:
        notebook='results/logs/notebooks/aggregate_pileup_analysis.ipynb'
    notebook:
        'notebooks/aggregate_pileup_analysis.py.ipynb'

rule sex_chromosome_counts:
    """Count reads mapping perfectly to each host genome sex chromosome."""
    input:
        bam=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner],
        bai=lambda wc: {'bbmap': rules.align_bbmap.output.bam,
                        'bwa-mem2': rules.align_bwa_mem2.output.bam,
                        'minimap2': rules.align_minimap2.output.bam,
                        }[wc.aligner] + '.bai',
    output:
        counts="results/sex_chromosome/{aligner}/{genome}/{sample}.csv",
    params:
        male=lambda wc: (config['host_genomes'][wc.genome]
                               ['sex_chromosomes']['male']),
        female=lambda wc: (config['host_genomes'][wc.genome]
                                 ['sex_chromosomes']['female']),
    conda: 'environment.yml'
    shell:
        # get just primary alignments with no mismatches:
        # -F 256 is primary only: https://www.biostars.org/p/259963/#259965
        # NM tag should be 0: https://www.biostars.org/p/148858/#148879
        # the grep command is set to avoid error on empty:
        # https://unix.stackexchange.com/a/427598
        r"""
        printf "sex,chromosome,count\n" > {output.counts}
        printf "male,{params.male}," >> {output.counts}
        samtools view -F 256 {input.bam} {params.male} | \
            ( grep -P '\bNM:i:0\b' || [[ $? == 1 ]] ) | \
            wc -l >> {output.counts}
        printf "female,{params.female}," >> {output.counts}
        samtools view -F 256 {input.bam} {params.female} | \
            ( grep -P '\bNM:i:0\b' || [[ $? == 1 ]] ) | \
            wc -l >> {output.counts}
        """

rule analyze_sex:
    """Analyze the sex of the patients."""
    input:
        counts=expand(rules.sex_chromosome_counts.output.counts,
                      aligner=config['aligners'],
                      genome=config['host_genomes'],
                      sample=samples)
    output:
        stats=report('results/sex_chromosome/stats.csv',
                     category='Sex of patients',
                     caption='report/analyze_sex_stats.rst'),
        chart=report('results/sex_chromosome/chart.html',
                     category='Sex of patients',
                     caption='report/analyze_sex_chart.rst'),
    params:
        descriptors=[{'aligner': aligner,
                      'host_genome': host_genome,
                      'sample': sample}
                     for aligner, host_genome, sample
                     in itertools.product(config['aligners'],
                                          config['host_genomes'],
                                          samples)
                     ],
        reported_sex={sample: sample_d['sex'] if 'sex' in sample_d else 'unknown'
                      for sample, sample_d in samples.items()},
    log:
        notebook='results/logs/notebooks/analyze_sex.ipynb'
    conda: 'environment.yml'
    notebook:
        'notebooks/analyze_sex.py.ipynb'

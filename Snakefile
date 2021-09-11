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

# get all samples that are not plasmid
samples = {sample: sample_info for sample, sample_info in
           list(config['samples_fasterq_dump'].items()) +
           list(config['samples_wget'].items())
           if sample_info['patient_group'] != 'plasmid'
           }
# get all samples that are plasmid
plasmid_samples = {sample: sample_info for sample, sample_info in
                   list(config['samples_fasterq_dump'].items()) +
                   list(config['samples_wget'].items())
                   if sample_info['patient_group'] == 'plasmid'
                   }
# map all samples to accessions
samples_to_accessions = {sample: sample_info['accessions'] for
                         sample, sample_info in
                         list(samples.items()) + list(plasmid_samples.items())}

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
        expand("results/sra_downloads/{accession}.fastq.gz",
               accession=[acc for d in samples.values()
                          for acc in d['accessions']]),
        'results/pileup/coverage_all.html',
        'results/pileup/coverage_all.pdf',
        'results/pileup/coverage_region.html',
        'results/pileup/coverage_region.pdf',
        'results/pileup/samtools_pileup.csv',
        'results/early_sequences/annotated_filtered_substitutions.csv',
        'results/sra_file_info.csv',
        multiext('results/early_sequences/counts', '.html', '.pdf'),
        multiext('results/early_sequences/deltadist', '.html', '.pdf'),
        'results/early_sequences/deltadist.csv',
        multiext('results/early_sequences/deltadist_region', '.html', '.pdf'),
        multiext('results/deltadist_jitter', '.html', '.pdf'),
        'results/deleted_diffs.tex',
        'results/recovered_seqs.fa',
        'results/phylogenetics/tree_images/',
        'results/pileup/plasmid_mutations.csv',
        'results/deleted_seq_alignment_matches.csv',
        'results/phylogenetics/all_alignment_no_filter_rare.fa',
        'results/phylogenetics/all_alignment_no_filter_rare.csv',
        'results/consensus/consensus_seqs.csv',
        '_temp.txt',

checkpoint list_bigd_ftp:
    """List details from BIGD GSA FTP site."""
    output:
        files='results/bigd_files/files.txt',
        dirs='results/bigd_files/dirs.txt'
    params:
        host=config['bigd_host'],
        path=config['bigd_path'],
        user=config['bigd_user'],
        password=config['bigd_password'],
    conda: 'environment.yml'
    script: 'scripts/list_ftp.py'

def bigd_fastq_list(wildcards):
    """List of FASTQs to get from BIGD."""
    with checkpoints.list_bigd_ftp.get(**wildcards).output.dirs.open() as f:
        accessions = [os.path.basename(line.strip())
                      for line in f.readlines()]
    return [f"results/bigd_files/{acc}.fastq.gz" for acc in accessions]

rule get_bigd_fastq:
    """Download a file from BIGD (GSA)."""
    output: fastq="results/bigd_files/{bigd_acc}.fastq.gz"
    params:
        path=os.path.join('ftp://' + config['bigd_host'],
                          config['bigd_path'],
                          "{bigd_acc}",
                          "{bigd_acc}.fastq.gz")
    conda: 'environment.yml'
    shell: "wget {params.path} -O {output.fastq}"

rule aggregate_bigd_fastq:
    input: bigd_fastq_list
    output: '_temp.txt'
    conda: 'environment.yml'
    shell: 'echo "not implemented"'

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
                    'https://sra-pub-run-odp.s3.amazonaws.com/sra'],
        alternative_fastq_path=lambda wc: (config['alternative_fastq_locations'][wc.accession]
                                           if wc.accession in config['alternative_fastq_locations']
                                           else '')
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        if [[ "{params.alternative_fastq_path}" != "" ]]
        then
            echo "getting FASTQ for {wildcards.accession} from {params.alternative_fastq_path}"
            wget {params.alternative_fastq_path} -O {output.fastq_gz}
            mkdir -p {output.fastq_dir}
            mkdir -p {output.temp_dir}
            touch {output.sra_file}
        else
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
        fi
        """

rule sra_file_info:
    """Get info for ``*.sra`` files using ``vdb-dump``."""
    input:
        sra_files=expand(rules.download_sra.output.sra_file,
                         accession=[acc for d in samples.values()
                                    for acc in d['accessions']
                                    if acc not in config['alternative_fastq_locations']])
    output: csv='results/sra_file_info.csv'
    conda: 'environment.yml'
    script: 'scripts/sra_file_info.py'

rule preprocess_fastq:
    """Pre-process the FASTQ files by trimming adaptors etc with ``fastp``."""
    input:
        fastq_gz=rules.download_sra.output.fastq_gz,
    output:
        fastq_gz="results/preprocessed_fastqs/{accession}.fastq.gz",
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
            -l {params.min_read_length} \
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
                                 accession=samples_to_accessions[wc.sample]),
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
                                 accession=samples_to_accessions[wc.sample]),
        mmi=rules.minimap2_genome.output.mmi,
    output:
        sam=temp("results/alignments/minimap2/{sample}.sam"),
        unsorted_bam=temp("results/alignments/minimap2/{sample}.bam"),
        bam="results/alignments/minimap2/{sample}_sorted.bam",
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
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

rule adapter_sites:
    """Get all sites covered by adapters."""
    input:
        adapters=config['adapters_to_trim'],
        ref_genome=rules.get_ref_genome_fasta.output.fasta
    output: adapter_sites='results/adapter_sites/adapter_sites.csv'
    conda: 'environment.yml'
    script: 'scripts/adapter_sites.py'

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
        ref_fasta=rules.get_ref_genome_fasta.output.fasta,
        adapter_sites=rules.adapter_sites.output.adapter_sites,
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
            --exclude_sites {input.adapter_sites}
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
    params: gisaid_dirs=config['gisaid_comparator_dirs']
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
        site_map='results/genome_to_comparator/site_identity_map.csv'
    params:
        comparators=list(config['comparator_genomes'])
    conda: 'environment.yml'
    script:
        'scripts/genome_comparator_map.py'

rule analyze_plasmid_seqs:
    """Analyze pileups for plasmid samples."""
    input:
        pileups=expand(rules.bam_pileup.output.pileup_csv,
                       aligner=config['aligners'],
                       sample=plasmid_samples),
    output: plasmid_muts='results/pileup/plasmid_mutations.csv'
    params:
        consensus_min_frac=config['consensus_min_frac'],
        consensus_min_coverage=config['consensus_min_coverage'],
        descriptors=[{'aligner': aligner, 'sample': sample}
                     for aligner, sample in
                     itertools.product(config['aligners'], plasmid_samples)],
    conda: 'environment.yml'
    log: notebook='results/logs/notebooks/analyze_plasmid_seqs.ipynb'
    notebook: 'notebooks/analyze_plasmid_seqs.py.ipynb'

rule analyze_pileups:
    """Analyze and plot BAM pileups per sample."""
    input:
        pileups=expand(rules.bam_pileup.output.pileup_csv,
                       aligner=config['aligners'],
                       allow_missing=True)
    output:
        chart="results/pileup/{sample}/interactive_pileup.html",
        diffs_from_ref="results/pileup/{sample}/diffs_from_ref.csv",
        pileup_csv="results/pileup/{sample}/samtools_pileup.csv",
    params:
        consensus_min_frac=config['consensus_min_frac'],
        consensus_min_coverage=config['consensus_min_coverage'],
        descriptors=[{'aligner': aligner} for aligner in config['aligners']],
        chart_title="{sample}"
    conda: 'environment.yml'
    log: notebook="results/logs/notebooks/analyze_pileups_{sample}.ipynb"
    notebook: 'notebooks/analyze_pileups.py.ipynb'

rule plot_aggregate_pileup:
    """Plot pileups across all samples."""
    input:
        pileups=expand(rules.analyze_pileups.output.pileup_csv,
                       sample=samples),
    output:
        chart_all_html='results/pileup/coverage_all.html',
        chart_all_pdf='results/pileup/coverage_all.pdf',
        chart_region_html='results/pileup/coverage_region.html',
        chart_region_pdf='results/pileup/coverage_region.pdf',
        csv='results/pileup/samtools_pileup.csv',
    params:
        samples=list(samples),
        patient_groups=[d['patient_group'] for d in samples.values()],
        region_of_interest=config['region_of_interest'],
        aligners=config['aligners'],
        ref_name=config['ref_genome']['name'],
        consensus_min_coverage=config['consensus_min_coverage'],
    conda: 'environment.yml'
    log: notebook="results/logs/notebooks/plot_aggregate_pileup.ipynb"
    notebook: 'notebooks/plot_aggregate_pileup.py.ipynb'

rule diffs_from_ref:
    """Analyze differences from reference aggregated across deep sequencing samples."""
    input:
        diffs_from_ref=expand(rules.analyze_pileups.output.diffs_from_ref,
                              sample=samples),
        comparator_map=rules.genome_comparator_map.output.site_map,
    output:
        diffs_from_ref_stats='results/pileup/diffs_from_ref.csv',
        diffs_from_ref_chart='results/pileup/diffs_from_ref.html',
    params:
        samples=list(samples),
    conda: 'environment.yml'
    log: notebook='results/logs/notebooks/diffs_from_ref.ipynb'
    notebook: 'notebooks/diffs_from_ref.py.ipynb'

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

rule early_seqs_gisaid_fasta:
    """Get FASTA from GISAID augur download."""
    input: lambda wc: config['early_seqs']['gisaid'][wc.sequence_set]
    output: fasta="results/early_sequences/{sequence_set}.fa",
    params: props=config['early_seq_header_props']
    conda: 'environment.yml'
    script: 'scripts/gisaid_subdir_to_fasta.py'

rule align_early_seqs:
    """Align all the early sequences to the reference, stripping gaps."""
    input:
        ref_genome=rules.get_ref_genome_fasta.output.fasta,
        gisaid_fastas=expand(rules.early_seqs_gisaid_fasta.output.fasta,
                             sequence_set=config['early_seqs']['gisaid']),
    output:
        alignment='results/early_sequences/full_alignment.fa',
        concat_early_seqs=temp('results/early_sequences/_concat_early_seqs.fa'),
    conda: 'environment.yml'
    shell:
        # concatenate with newline at end of each file:
        # https://stackoverflow.com/a/25030513/4191652
        # then align to reference:
        # https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
        """
        awk 1 {input.gisaid_fastas} > {output.concat_early_seqs}
        mafft \
            --6merpair \
            --keeplength \
            --addfragments {output.concat_early_seqs} \
            {input.ref_genome} \
            > {output.alignment}
        """

rule early_seq_subs:
    """Get substitution mutations from early sequences."""
    input:
        alignment=rules.align_early_seqs.output.alignment,
        ref_genome=rules.get_ref_genome_fasta.output.fasta,
    output: csv='results/early_sequences/substitutions.csv'
    params:
        region_of_interest=config['region_of_interest'],
        props=config['early_seq_header_props'],
        ignore_muts_before=config['early_seqs_ignore_muts_before'],
        ignore_muts_after=config['early_seqs_ignore_muts_after'],
    conda: 'environment.yml'
    script: 'scripts/early_seq_subs.py'

rule annotate_early_seq_subs:
    """Annotate and filter early sequence substitutions."""
    input:
        subs_csv=rules.early_seq_subs.output.csv,
        comparator_map=rules.genome_comparator_map.output.site_map,
        who_china_report_cases_yaml=config['who_china_report_cases'],
        early_seqs_to_exclude_yaml=config['early_seqs_to_exclude'],
        wuhan_exports_yaml=config['wuhan_exports'],
    output: csv='results/early_sequences/annotated_filtered_substitutions.csv'
    params:
        comparators=list(config['comparator_genomes']),
        min_coverage=config['early_seqs_min_coverage'],
        max_subs=config['early_seqs_max_subs'],
        max_ambiguous=config['early_seqs_max_ambiguous'],
        max_date=config['early_seqs_max_date'],
        filter_runs=config['early_seqs_filter_runs'],
        who_china_report_last_date=config['who_china_report_last_date'],
    conda: 'environment.yml'
    log: notebook='results/logs/notebooks/annotate_early_seq_subs.ipynb'
    notebook: 'notebooks/annotate_early_seq_subs.py.ipynb'

rule outgroup_dist_analysis:
    """Analyze distances to comparator outgroups."""
    input:
        early_seq_subs=rules.annotate_early_seq_subs.output.csv,
        early_seq_alignment=rules.align_early_seqs.output.alignment,
        comparator_map=rules.genome_comparator_map.output.site_map,
        deleted_diffs=rules.diffs_from_ref.output.diffs_from_ref_stats,
        deleted_consensus=rules.aggregate_consensus_seqs.output.csv,
    output:
        early_seq_counts=multiext('results/early_sequences/counts', '.html', '.pdf'),
        early_seq_deltadist=multiext('results/early_sequences/deltadist', '.html', '.pdf'),
        early_seq_deltadist_region=multiext('results/early_sequences/deltadist_region', '.html', '.pdf'),
        deltadist_jitter=multiext('results/deltadist_jitter', '.html', '.pdf'),
        deltadist_csv='results/early_sequences/deltadist.csv',
        deleted_diffs_latex='results/deleted_diffs.tex',
        recovered_seqs='results/recovered_seqs.fa',
        alignment_all_fasta='results/phylogenetics/all_alignment.fa',
        alignment_all_csv='results/phylogenetics/all_alignment.csv',
        alignment_all_no_filter_rare_fasta='results/phylogenetics/all_alignment_no_filter_rare.fa',
        alignment_all_no_filter_rare_csv='results/phylogenetics/all_alignment_no_filter_rare.csv',
        deleted_csv='results/phylogenetics/deleted_seqs.csv',
        matches_in_gisaid='results/deleted_seq_alignment_matches.csv',
    params:
        region_of_interest=config['region_of_interest'],
        comparators=list(config['comparator_genomes']),
        min_frac_coverage=config['min_frac_coverage'],
        samples=samples,
        aligners=config['aligners'],
        ref_genome_name=config['ref_genome']['name'],
        ignore_muts_before=config['early_seqs_ignore_muts_before'],
        ignore_muts_after=config['early_seqs_ignore_muts_after'],
        phylo_last_date=config['phylo_last_date'],
        phylo_muts_to_ignore=config['phylo_muts_to_ignore'],
        phylo_collapse_rare_muts=config['phylo_collapse_rare_muts'],
        phylo_filter_rare_variants=config['phylo_filter_rare_variants'],
        phylo_min_frac_called=config['phylo_min_frac_called'],
        cat_colors=config['cat_colors'],
        subcat_colors=config['subcat_colors'],
    conda: 'environment.yml'
    log: notebook='results/logs/notebooks/outgroup_dist_analysis.ipynb'
    notebook: 'notebooks/outgroup_dist_analysis.py.ipynb'

checkpoint possible_progenitors:
    """Get possible progenitors with smallest distance to outgroup."""
    input: all_csv=rules.outgroup_dist_analysis.output.alignment_all_csv
    output: progenitors='results/phylogenetics/progenitors.txt'
    params: outgroups=list(config['comparator_genomes']),
    conda: 'environment.yml'
    script: 'scripts/possible_progenitors.py'

def progenitor_trees(_):
    with checkpoints.possible_progenitors.get().output.progenitors.open() as f:
        progenitors = [line.strip().replace('/', '%') for line in f]
        return [f"results/phylogenetics/all_{progenitor}.treefile"
                for progenitor in progenitors]

def progenitor_states(_):
    with checkpoints.possible_progenitors.get().output.progenitors.open() as f:
        progenitors = [line.strip().replace('/', '%') for line in f]
        return [f"results/phylogenetics/all_{progenitor}.state"
                for progenitor in progenitors]

rule iqtree:
    """Infer ``iqtree`` phylogenetic tree."""
    input:
        alignment="results/phylogenetics/all_alignment.fa",
    output:
        treefile="results/phylogenetics/all_{progenitor}.treefile",
        ancestral="results/phylogenetics/all_{progenitor}.state",
    params:
        pre=lambda wc, output: os.path.splitext(output.treefile)[0],
        progenitor=lambda wc: wc.progenitor.replace('%', '/')
    threads: config['max_cpus']
    conda: 'environment.yml'
    shell:
        """
        iqtree \
            -s {input.alignment} \
            -pre {params.pre} \
            -st DNA \
            -m GTR+F \
            -czb \
            --keep-ident \
            -nt {threads} \
            -redo \
            -seed 2 \
            -o {params.progenitor} \
            -asr
        """

rule visualize_trees:
    """Visualize the phylogenetic trees."""
    input:
        trees=progenitor_trees,
        states=progenitor_states,
        alignment="results/phylogenetics/all_alignment.fa",
        all_csv=rules.outgroup_dist_analysis.output.alignment_all_csv,
        deleted_csv=rules.outgroup_dist_analysis.output.deleted_csv,
        comparator_map=rules.genome_comparator_map.output.site_map,
    output: directory('results/phylogenetics/tree_images')
    params:
        site_offset=config['early_seqs_ignore_muts_before'] - 1,
        progenitors=lambda wc, input: [os.path.splitext(os.path.basename(f))[0]
                                       .replace('all_', '')
                                       .replace('%', '/')
                                       for f in input.trees],
        outgroups=list(config['comparator_genomes']),
        region_of_interest=config['region_of_interest'],
        cat_colors=config['cat_colors'],
        subcat_colors=config['subcat_colors'],
        wuhan_hu_1_add_muts=config['ref_genome']['add_mutations']
    conda: 'environment_ete3.yml'
    log: notebook='results/logs/notebooks/visualize_trees.ipynb'
    notebook: 'notebooks/visualize_trees.py.ipynb'

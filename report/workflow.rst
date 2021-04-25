This report summarizes an analysis of early Wuhan SARS-CoV-2 sequences by Jesse Bloom.
The goal of the analysis is to re-analyze the primary deep sequencing data of early viral samples to see if they shed any light on the origin and spread of the virus.

Briefly, the deep sequencing data for the following samples were downloaded from the `NCBI Sequence Read Archive <https://www.ncbi.nlm.nih.gov/sra>`_:
{% for sample, sample_d in snakemake.config['samples'].items() %}
 - `{{sample}} <{{sample_d['study_url']}}>`_
{% endfor %}

The reads were then aligned to the following reference genome(s):
{% for genome in snakemake.config['genomes'] %}
 - {{genome}}
{% endfor %}

The alignments were done using following short-read aligners.
{% for aligner in snakemake.config['aligners'] %}
 - {{aligner}}
{% endfor %}

Any base calls corresponding to read sites with a Q-score of at least {{snakemake.config['minq']}} were retained, lower-quality base calls were ignored for the rest of the analysis.
The pileup files from these alignments are detailed in Pileups_, which also reports any mismatches of the deep sequencing cnsensus relative to the reference genomes.
Consensus sequences from the deep sequencing were called by requiring a depth of {{ snakemake.config.min_consensus_coverage }} and more than {{ snakemake.config.min_consensus_coverage }} of the base calls to agree on the consensus identity.
Indels are ignored, so the analysis is only identifying nucleotide substitutions.

The consensus from the deep sequencing was then aligned to the Genbank sequence reported for each sample to see if any sites differed in the deep sequencing from what was reported as the Genbank sequence.
These differences are summarized in `Deep sequencing vs Genbank`_.

The full workflow is schematized below and the rest of this report contains more details about the configuration, etc.

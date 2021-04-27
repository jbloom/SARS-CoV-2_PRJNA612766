This report summarizes an analysis of early SARS-CoV-2 sequences by Jesse Bloom.
The goal of the analysis is to re-analyze the primary deep sequencing data of early viral samples to see if they shed any light on the origin and spread of the virus.

Briefly, the deep sequencing data for the following samples were downloaded from the `NCBI Sequence Read Archive <https://www.ncbi.nlm.nih.gov/sra>`_ (see the *Configuration* section for more details):
{% for sample, sample_d in snakemake.config['samples'].items() %}
 - `{{sample}} <{{sample_d['study_url']}}>`_
{% endfor %}

The reads were then aligned to the following reference viral genome(s):
{% for genome in snakemake.config['genomes'] %}
 - {{genome}}
{% endfor %}

The alignments were done using following short-read aligners.
{% for aligner in snakemake.config['aligners'] %}
 - {{aligner}}
{% endfor %}

For analysis of the viral deep sequencing, any base calls corresponding to read sites with a Q-score of at least {{snakemake.config['minq']}} were retained, lower-quality base calls were ignored for the rest of the analysis.
Consensus sequences from the deep sequencing were called by requiring a depth of {{ snakemake.config['min_consensus_coverage'] }} and more than {{ snakemake.config['min_consensus_coverage'] }} of the base calls to agree on the consensus identity.
Indels relative to the consensus are ignored, so the analysis is only identifying nucleotide substitutions.
The results are summarized in `Viral deep sequencing analysis`_.

The consensus from the deep sequencing was then aligned to the Genbank sequence reported for each sample to see if any sites differed in the deep sequencing from what was reported as the Genbank sequence.
At sites where there were differences, the identities in closely related comparator coronavirus genomes were also analyzed.
These deep-sequencing consensus versus Genbank comparisons are summarized in `Deep sequencing vs Genbank`_.

In addition, reads were mapped to the host genome (requiring exact primary matches) and relative mapping to the sex chromosomes was used to assess the sex (e.g., male or female) of the patients for which the viral reads were derived.
These results are summarized in `Sex of patients`_.

The full workflow is schematized below and the rest of this report contains more details about the configuration, etc.

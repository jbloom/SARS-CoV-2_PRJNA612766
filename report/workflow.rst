This report summarizes an analysis of early SARS-CoV-2 sequences by Jesse Bloom.
The goal of the analysis is to re-analyze the primary deep sequencing data of early viral samples to see if they shed any light on the origin and spread of the virus.

The workflow briefly works like this:

 - Deep sequencing data for viral samples are downloaded from the `NCBI Sequence Read Archive <https://www.ncbi.nlm.nih.gov/sra>`_

 - The deep sequencing data are aligned to the reference viral genome(s) using short-read aligner(s).

 - Pileup files are generated from the alignments ingoring any base calls with a Q score <{{snakemake.config['minq']}}.
   Indels are also ignored (so this workflow only analyzes substitutions).
   Consensus sequences from the deep sequencing are then called by requiring a depth of at least {{ snakemake.config['consensus_min_coverage'] }} and more than {{ snakemake.config['consensus_min_frac'] }} of the base calls to agree on the consensus identity.
   The results are summarized in `Viral deep sequencing analysis`_, which includes analyses of:

     - Coverage statistics over the viral genome.
     - Mutations relative to the reference viral genome.
     - Whether the mutations are to nucleotide identities found in other related comparator viral genomes.

The *Configuration* section has details about the samples, read aligners, reference genomes, and other settings.

The graph immediately below schematizes the workflow.

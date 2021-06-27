# Manually downloaded input data

## GISAID sequences
Due to [GISAID data sharing terms](https://www.gisaid.org/help/faq/), the actual FASTA sequences are not tracked in this repo.
However, the subdirectories with GISAID sequences contain `*.tsv.xz` files that list the accessions, as well as acknowledging the originating and submitting labs.
These are tracked.

- [./comparator_genomes_gisaid/](comparator_genomes_gisaid) contains bat virus comparator genomes as downloaded from GISAID on May-31-2021.

- [./gisaid_sequences_through_Feb2020/](gisaid_sequences_through_Feb2020) contains all sequences downloaded from GISAID on May-30-2021 from human SARS-CoV-2 isolated no later than Feb-29-2020.

## Annotations of early sequences
- [WHO_China_Report_Dec2019_cases.yaml](WHO_China_Report_Dec2019_cases.yaml) defines the cases from patients with onset date before Dec-31-2019 as defined in Tables 6 and 7 of the [joint China-WHO COVID-19 origins study](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part). Then looks at these sequences and decides which ones to "collapse" to a single sequence (ones that report says are from same patient). Also annotated which early sequences are from Huanan seafood market.

- [seqs_to_exclude.yaml](seqs_to_exclude.yaml) defines sequences manually for exclusion based on various criteria explained in the file.

- [Wuhan_exports.yaml](Wuhan_exports.yaml) defines sequences for which there is metadata indicating patient was infected in Wuhan even if sequenced elsewhere.

## Primers for adaptor trimming
- [Wang_et_al_primers.fasta](Wang_et_al_primers.fasta) gives the sequences of the primers used to amplify the sequenced regions, which should be trimmed from the data prior to alignment. These are taken from [Supplementary Table 1 of Wang et al (2020, medRxiv)](https://www.medrxiv.org/content/medrxiv/suppl/2020/03/06/2020.03.04.20029538.DC1/2020.03.04.20029538-1.pdf).

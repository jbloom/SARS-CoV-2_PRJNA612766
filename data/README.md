# Manually downloaded input data

## GISAID sequences
Due to [GISAID data sharing terms](https://www.gisaid.org/help/faq/), the actual FASTA sequences are not tracked in this repo.
However, the subdirectories with GISAID sequences contain `*.tsv.xz` files that list the accessions, as well as acknowledging the originating and submitting labs.
These are tracked.

- [./comparator_genomes_gisaid/](comparator_genomes_gisaid) contains bat virus comparator genomes as downloaded from GISAID on May-31-2021.

- [./gisaid_sequences_through_Feb2020/](gisaid_sequences_through_Feb2020) contains all sequences downloaded from GISAID on May-30-2021 from human SARS-CoV-2 isolated no later than Feb-29-2020.

## Annotations of early sequences
- [WHO_China_Report_Dec2019_cases.yaml](WHO_China_Report_Dec2019_cases.yaml) defines the cases from patients with onset date before Dec-31-2019 as defined in Tables 6 and 7 of the [joint China-WHO COVID-19 origins study](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part). Then looks at these sequences and decides which ones to "collapse" to a single sequence (ones that report says are from same patient).

- [seqs_to_exclude.yaml](seqs_to_exclude.yaml) defined sequences manually for exclusion based on various criteria explained in the file.

## Mutations extracted from GISAID data by Richard Neher

The [gisaid_mutations](gisaid_mutations) subdirectory contains data on mutations relative to Wuhan-Hu-1 for all full-length sequences in GISAID, compiled and provided to me by Richard Neher on May-17-2021 (not sure when he scraped GISAID for it).
I'm not sure if I'm allowed to track the mutations file on GitHub, so for now I just track the metadata file with information on the originating and submitting labs, accessions, etc.

# Recovery of deleted deep sequencing data sheds more light on the early Wuhan SARS-CoV-2 epidemic
This GitHub repository analyzes SARS-CoV-2 deep sequencing data recovered from the deleted BioProject PRJNA612766.
This analysis corresponds to the work described in [this paper](https://doi.org/10.1093/molbev/msab246).

Specifically:

 - [this tag](https://github.com/jbloom/SARS-CoV-2_PRJNA612766/tree/initial_bioRxiv_version) corresponds to the [initial bioRxiv pre-print](https://www.biorxiv.org/content/10.1101/2021.06.18.449051v1)
 - [this tag](https://github.com/jbloom/SARS-CoV-2_PRJNA612766/tree/second_bioRxiv_version) corresponds to the [revised bioRxiv pre-print](https://www.biorxiv.org/content/10.1101/2021.06.18.449051v2)
 - [this tag](https://github.com/jbloom/SARS-CoV-2_PRJNA612766/tree/published_MBE_version) corresponds to the [final version accepted by MBE](https://doi.org/10.1093/molbev/msab246)

## Running the analysis
The analysis is nearly fully automated by the `snakemake` pipeline included in [Snakefile](Snakefile).
The configuration for the analysis is in [config.yaml](config.yaml).
Note that the pipeline is somewhat convoluted and performs a variety of steps only tangentially related to the paper corresponding to this study.
The reason is that the study started simply as an effort to validate the analyses in the [joint WHO-China report on COVID-19 origins](https://www.who.int/publications/i/item/who-convened-global-study-of-origins-of-sars-cov-2-china-part), but then gradually shifted in goal upon the discovery of the deleted data set.
For this reason, there are still some vestigial parts of the code and analysis structure.

The only required manual step is to download existing coronavirus sequences from [GISAID](https://www.gisaid.org/), which must be done manually after creating a [GISAID](https://www.gisaid.org/) account since [GISAID data sharing terms](https://www.gisaid.org/help/faq/) prevent distribution of their sequences.
To get these sequences, download both the `*.metadata.tsv.xz` and `*.fasta.xz` files for the accessions in [data/gisaid_sequences_through_Feb2020/accessions.txt](data/gisaid_sequences_through_Feb2020/accessions.txt) to the subdirectory [data/gisaid_sequences_through_Feb2020/](data/gisaid_sequences_through_Feb2020/), and the same two files for the accessions in [data/comparator_genomes_gisaid/accessions.txt](data/comparator_genomes_gisaid/accessions.txt) to the subdirectory [data/comparator_genomes_gisaid/](data/comparator_genomes_gisaid/). 

After downloading these sequences and ensuring you have installed [conda](https://docs.conda.io/en/latest/), build the main [conda](https://docs.conda.io/en/latest/) environment for the pipeline with:

    conda env create -f environment.yml

Then activate the [conda](https://docs.conda.io/en/latest/) environment with:

    conda activate SARS-CoV-2_PRJNA612766

You can then run the entire analysis with:

    snakemake -j 1 --use-conda

Note that you need the `--use-conda` command because one of the rules in [Snakefile](Snakefile) uses a separate environment as specified in [environment_ete3.yml](environment_ete3.yml).

The above command will run the `snakemake` pipeline using just one computing core.
If you want to use more cores, adjust the value passed by `-j` appropriately.
If you have access to a computing cluster you can distribute the run across the cluster.
For the Fred Hutch computing cluster, that can be done using [cluster.yaml](cluster.yaml) by running the pipeline with the commands in [run_Hutch_cluster.bash](run_Hutch_cluster.bash).

## Input data, results, etc
The input data needed for the analysis are all available in the [./data/](data) subdirectory, which contains a README describing the files therein.

The results of running the pipeline are placed in the [./results/](results) subdirectory.
Most of these results are not tracked in this GitHub repo, but some key files are as described in the Methods of the paper associated with this study.

The code used to process the Excel supplementary table of accessions from project PRJNA612766 to generate the information found in [config.yaml](config.yaml) is in [./manual_analyses/PRJNA612766/](manual_analyses/PRJNA612766).

## Paper
The LaTex source for the paper and its figures are found in the [./paper/](paper) subdirectory.

## Updates to repo in late 2024 at request of GISAID compliance
In late 2024, I received a request from GISAID compliance to remove certain files in which they said they had been made aware of files in the repo that violated their data sharing.
I therefore made the following changes:

On-Nov-21-2024, I made commit `e9d972519` which removed these files from the current state of the repository.

On Dec-7-2024, I fully removed the files from the `git` history with the following commands:
```
git filter-repo --path data/gisaid_sequences_through_Feb2020/1622384383620.metadata.tsv.xz --invert-paths
git filter-repo --path data/comparator_genomes_gisaid/1622468911409.metadata.tsv.xz --invert-paths
git filter-repo --path results/early_sequences/deltadist.csv --invert-paths
```
I then made and committed new version of [results/early_sequences/deltadist.csv](results/early_sequences/deltadist.csv) where I had removed all columns that contained mutations, and committed it and this update to the README.
I then force-pushed these changes with `git push --force`.

## Updates to repo on May-2025 at request of GISAID compliance
On May-20-2025, I received a request from GISAID compliance to remove additoinal files that they said had GISAID metadata, and about which they had received a complaint from a contributor of the data.
Note that these were not files in the current state of the repo, but apparently were accessible in the history.

On May-21-2025, therefore made the following changes to fully remove the files from the `git` history:
```
git filter-repo --path data/comparator_genomes_gisaid/1621180534800.metadata.tsv.xz --invert-paths
git filter-repo --path data/gisaid_mutations/metadata.tsv.gz --invert-paths
```

I then force-pushed these changes with `git push --force origin main`.

## GISAID EpiSet identifier
At the time this project was performed in 2021, GISAID had not yet introduced EpiSet identifiers.

However, on May-29-2025, GISAID Compliance generated an EpiSet identifier of `EPI_SET_210531sk` (doi: [10.55876/gis8.210531sk](https://epicov.org/epi3/epi_set/210531sk?main=true)) that I have now added here to acknowledge the GISAID submitters.

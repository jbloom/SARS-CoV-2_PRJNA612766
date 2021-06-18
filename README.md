# Recovery of deleted deep sequencing data sheds more light on the early Wuhan SARS-CoV-2 epidemic
This GitHub repository analyzes SARS-CoV-2 deep sequencing data recovered from the deleted BioProject PRJNA612766.

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

# Accessions from BioProject PRJNA612766

[Farkas et al (2020)](https://peerj.com/articles/9255) include as Supplementary Table 1 of their paper all the Sequence Read Archive accessions relevant to SARS-CoV-2 that were present as of March-27-2020.
The Excel formatted [Supplementary_Table_1.xlsx](Supplementary_Table_1.xlsx) is included here.

Notably, this table lists a large number of accessions from BioProject PRJNA612766 by Wuhan University.
However, as of April-29-2021, this BioProject is not accessible on the SRA web interface.
However, the SRA files can still be obtained using commands like this:

    wget https://storage.googleapis.com/nih-sequence-read-archive/run/SRR11313395/SRR11313395

The Excel Supplementary Table 1 indicates the plasmid control samples are from Aisi Fu and the human samples from Renmin Hospital.
Googling these terms brings up [this pre-print](https://www.medrxiv.org/content/10.1101/2020.03.04.20029538v1) (posted March-6-2020) and [this paper](https://onlinelibrary.wiley.com/doi/full/10.1002/smll.202002169) (published June-24-2020), which clearly describe the study here.
Reading that paper makes the following clear:

 - The samples that are listed as plasmid controls from Aisi Fu are indeed controls, and so aren't used in any downstream analyses here.

 - They tested 45 swab samples from outpatients with suspected COVID-19 early in the outbreak. The pre-print says "45 nasopharyngeal swab samples from outpatients with suspected COVID-19 early in the epidemic." The final published paper changes this to read "45 nasopharyngeal swab samples from outpatients with suspected COVID-19 early in the epidemic (January 2020)". Both the pre-print and final paper agree that 34 of these 45 samples were positive using their method. They also say that they tested the method on 16 positive samples from hospitalized patients with confirmed COVID-19 that had been tested by qPCR on Feb-11-2020 and Feb-12-2020 (both the pre-print and final paper fully agree on details of these samples); all 16 of these confirmed COVID-19 samples were positive using their method.

 - The text (e.g., legend of Figure 4 in pre-print) makes clear that "NC" and "PC" refer to positive and negative controls, and can so be excluded.

 - The -4h and -10min samples are different reads of the same sample, and so can be combined.

 - The samples with names "R01" to "R16" are the hospitalized confirmed patients, and the samples with well names (C2, E5, etc) are the outpatients (see Figure 4 of pre-print).

 - The following six samples were negative, and so are excluded: D2, F12, D10, A5, H3, A10 (see Figure 4 of pre-print).

 - The following samples were inconclusive, and so are excluded: A3, A7, A8, G6, B5 (see Figure 4 of pre-print)

 - In the published paper, Table 1 gives all called mutations.

This directory consists of a Jupyter notebook [extract_accessions.ipynb](extract_accessions.ipynb) that extracts all of these accessions from [Supplementary_Table_1.xlsx](Supplementary_Table_1.xlsx), and writes them in YAML format to [for_config.yml](for_config.yml) in the format needed for the configuration file for the main analysis.
Note it only retains the swab samples from outpatients and confirmed hospitalized cases that were positive, grouping both timepoints of the same samples.

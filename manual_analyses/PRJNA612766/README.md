# Accessions from BioProject PRJNA612766

[Farkas et al (2020)](https://peerj.com/articles/9255) include as Supplementary Table 1 of their paper all the Sequence Read Archive accessions relevant to SARS-CoV-2 that were present as of March-27-2020.
The Excel formatted [Supplementary_Table_1.xlsx](Supplementary_Table_1.xlsx) is included here.

Notably, this table lists a large number of accessions from BioProject PRJNA612766 by Wuhan University.
However, as of April-29-2021, this BioProject is not accessible on the SRA.
However, the SRA files can still be obtained using commands like this:

    wget https://storage.googleapis.com/nih-sequence-read-archive/run/SRR11313395/SRR11313395

This directory consists of a Jupyter notebook [extract_accessions.ipynb](extract_accessions.ipynb) that extracts all of these accessions from [Supplementary_Table_1.xlsx](Supplementary_Table_1.xlsx).

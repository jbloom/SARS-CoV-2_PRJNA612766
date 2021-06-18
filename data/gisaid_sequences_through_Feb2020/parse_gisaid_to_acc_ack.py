"""Parse xzipped GISAID TSV to accessions and acknowledgments."""


import lzma

import pandas as pd


with lzma.open('1622384383620.metadata.tsv.xz') as f:
    df = (pd.read_csv(f, sep='\t')
          # next two lines handle entries with newlines
          .query('gisaid_epi_isl.notnull()')
          .query('gisaid_epi_isl.str.startswith("EPI_ISL_")')
          .reset_index(drop=True)
          )

accessions = df['gisaid_epi_isl'].tolist()
assert len(accessions) == len(set(accessions))
acc_file = 'accessions.txt'
print(f"Writing the {len(accessions)} to {acc_file}")
with open(acc_file, 'w') as f:
    f.write('\n'.join(accessions))

acknowledgments = (
    df
    [['originating_lab', 'submitting_lab', 'authors']]
    .drop_duplicates()
    )
ack_file = 'acknowledgments.txt'
print(f"Writing the {len(acknowledgments)} to {ack_file}")
acknowledgments.to_csv(ack_file, index=False)

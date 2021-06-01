"""Convert subdir with GISAID augur download into FASTA."""


import glob
import lzma
import os
import sys

import Bio.SeqIO

import pandas as pd


subdir = snakemake.input[0]
print(f"Reading GISAID data from {subdir}")

sequences_file = glob.glob(f"{subdir}/*.sequences.fasta.xz")
if len(sequences_file) != 1:
    raise IOError(f"Found {len(sequences_file)} `*.sequences.fasta.xz` files")
sequences_file = sequences_file[0]
print(f"The GISAID sequences are in {sequences_file}")
        
metadata_file = sequences_file.replace('.sequences.fasta.xz', '.metadata.tsv.xz')
if not os.path.isfile(metadata_file):
    raise IOError(f"Failed to find metadata file {metadata_file}")
print(f"The GISAID metadata are in {metadata_file}")

with lzma.open(metadata_file) as f:
    metadata = (pd.read_csv(f, sep='\t')
                # next two lines handle entries with newlines
                .query('gisaid_epi_isl.notnull()')
                .query('gisaid_epi_isl.str.startswith("EPI_ISL_")')
                .reset_index(drop=True)
                )
assert len(metadata) == metadata['gisaid_epi_isl'].nunique()

with lzma.open(sequences_file, mode='rt') as f:
    iseq = 0
    seqs = []
    for seq in Bio.SeqIO.parse(f, 'fasta'):
        strain = metadata.at[iseq, 'strain']
        assert seq.id == strain, f"{seq.id=}\n{strain=}"
        seq.description = ', '.join(f"{prop}={metadata.at[iseq, prop]}"
                                   for prop in snakemake.params.props)
        iseq += 1
        seqs.append(seq)
print(f"Read a total of {iseq} sequences")

print(f"Writing sequences to {snakemake.output.fasta}")
Bio.SeqIO.write(seqs, snakemake.output.fasta, 'fasta')

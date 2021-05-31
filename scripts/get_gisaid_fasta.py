"""Implements ``snakemake`` rule ``get_gisaid_fasta``."""


import glob
import lzma
import os
import sys

import Bio.SeqIO

import pandas as pd


gisaid_id = snakemake.wildcards.gisaid
print(f"Looking for GISAID ID {gisaid_id}")

for subdir in snakemake.config['gisaid_comparator_dirs']:
    for metadata in glob.glob(f"{subdir}/*.metadata.tsv.xz"):
        print(f"Looking in {metadata}")
        fasta_in = metadata.replace('.metadata.tsv.xz', '.sequences.fasta.xz')
        if not os.path.isfile(fasta_in):
            raise IOError(f"cannot find {fasta_in}")
        with lzma.open(metadata) as f:
            df = pd.read_csv(f, sep='\t')
        df = df.query('gisaid_epi_isl == @gisaid_id')
        if len(df) == 0:
            continue
        elif len(df) == 1:
            idx = df.index.values[0]
            strain = df['strain'].values[0]
            print(f"Looking for {gisaid_id} ({strain}) in {fasta_in}")
            with lzma.open(fasta_in, mode='rt') as f:
                for i, seq in enumerate(Bio.SeqIO.parse(f, 'fasta')):
                    if i == idx:
                        if seq.id != strain:
                            raise ValueError(f"expected id of {strain}, got {seq.id}")
                        print(f"Writing FASTA to {snakemake.output.fasta}")
                        Bio.SeqIO.write([seq], snakemake.output.fasta, 'fasta')
                        sys.exit(0)
                else:
                    raise ValueError(f"Failed to find {gisaid_id} in {fasta_in}")
        else:
            raise ValueError(f"Multiple matches for {gisaid_id} in {metadata}")
else:
    raise ValueError(f"Failed to find {gisaid_id} in {config['gisaid_comparator_dirs']}")

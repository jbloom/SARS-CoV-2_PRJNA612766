"""Implements ``snakemake`` rule `get_ref_genome_fasta`."""


import gzip
import io
import urllib.request

import Bio.SeqIO


with urllib.request.urlopen(snakemake.params.url) as f:
    contents = gzip.GzipFile(fileobj=f).read().decode('utf-8')
    
seq = str(Bio.SeqIO.read(io.StringIO(contents), 'fasta').seq)

for mut in snakemake.params.add_mutations:
    wt = mut[0]
    m = mut[-1]
    site = int(mut[1: -1]) - 1
    if seq[site] != wt:
        raise ValueError(f"Mismatch for {mut}, sequence has {seq[site]}")
    else:
        seq = seq[: site] + m + seq[site + 1:]

with open(snakemake.output.fasta, 'w') as f:
    f.write(f">{snakemake.params.name}\n{seq}\n")

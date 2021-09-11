"""Extract reads, alphabetize, write to FASTA with sequential headers."""


import gzip

import Bio.SeqIO


fastq = snakemake.input.fastq
fasta = snakemake.output.fasta

reads = []
with gzip.open(fastq, 'rt') as f:
    for seq in Bio.SeqIO.parse(f, 'fastq'):
        reads.append(str(seq.seq).upper())
reads = sorted(reads)

with open(fasta, 'w') as f:
    for i, read in enumerate(reads):
        f.write(f">seq_{i + 1}\n{read}\n")

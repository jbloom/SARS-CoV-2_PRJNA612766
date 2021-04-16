"""Script for rule ``scripts/trim3_polyA.py``."""


import Bio.SeqIO

seq = Bio.SeqIO.read(snakemake.input.fasta, 'fasta')
while seq.seq[-1] in ['A', 'a']:
    seq.seq = seq.seq[: -1]
Bio.SeqIO.write([seq], snakemake.output.fasta, 'fasta')

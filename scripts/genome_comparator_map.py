"""Build map of identities in reference to comparators from alignment."""

import Bio.SeqIO

import pandas as pd


# variables frome snakemake
alignment = snakemake.input.alignment
site_map = snakemake.output.site_map
comparators = snakemake.params.comparators

# read alignment
aligned_seqs = [str(s.seq).upper() for s in Bio.SeqIO.parse(alignment,
                                                            'fasta')]
assert len(aligned_seqs) == 1 + len(comparators)

# strip to reference
stripped_seqs = []
for s in aligned_seqs:
    stripped_s = []
    for nt, nt_ref in zip(s, aligned_seqs[0]):
        if nt_ref != '-':
            stripped_s.append(nt)
    stripped_seqs.append(''.join(stripped_s))

# make site map data frame and write to CSV
assert 'reference' not in comparators
assert 'site' not in comparators
(pd.DataFrame({name: list(s) for name, s in zip(['reference'] + comparators,
                                                stripped_seqs)
               })
 .assign(site=lambda x: x.index + 1)
 [['site', 'reference', *comparators]]
 .to_csv(site_map, index=False)
 )

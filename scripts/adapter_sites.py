"""Get adapter sites."""


import Bio.SeqIO

import pandas as pd


adapters = {s.id: str(s.seq) for s in Bio.SeqIO.parse(snakemake.input.adapters, 'fasta')}
adapters_rc = {s.id: str(s.seq.reverse_complement()) for s in Bio.SeqIO.parse(snakemake.input.adapters, 'fasta')}

ref_seq = str(Bio.SeqIO.read(snakemake.input.ref_genome, 'fasta').seq)

sites = {}

for a_name, a in list(adapters.items()) + list(adapters_rc.items()):
    if a in ref_seq:
        for r in range(ref_seq.index(a), ref_seq.index(a) + len(a)):
            assert r + 1 not in sites
            sites[r + 1] = a_name

(pd.Series(sites)
 .rename('primer')
 .rename_axis('site')
 .reset_index()
 .to_csv(snakemake.output.adapter_sites, index=False)
 )

"""Implements ``snakemake`` rule `early_seq_subs.py`."""


import Bio.SeqIO

import pandas as pd


print(f"Reading reference sequence from {snakemake.input.ref_genome}")
ref_genome = Bio.SeqIO.read(snakemake.input.ref_genome, 'fasta')
print(f"Reference genome is {ref_genome.id}")
length = len(ref_genome)
ref_seqstr = str(ref_genome.seq).upper()

print(f"Reading aligned sequences from {snakemake.input.alignment}")
iseq = 0
data = {'strain': [],
        'date': [],
        'gisaid_epi_isl': [],
        'genbank_accession': [],
        'n_gapped_to_ref': [],
        'n_ident_to_ref': [],
        'n_subs_to_ref': [],
        'n_ambiguous_to_ref': [],
        'substitutions': [],
        }
for seq in Bio.SeqIO.parse(snakemake.input.alignment, 'fasta'):
    if iseq == 0:
        if seq.id == ref_genome.id and seq.seq == ref_genome.seq:
            print('First sequence in alignment is reference, so skipping.')
    else:
        assert len(seq) == length
        seqstr = str(seq.seq).upper()
        data['n_gapped_to_ref'].append(seqstr.count('-'))
        data['n_ident_to_ref'].append(sum(ref_nt == seq_nt for ref_nt, seq_nt
                                          in zip(ref_seqstr, seqstr)))
        data['n_ambiguous_to_ref'].append(sum(nt not in {'-', 'A', 'C', 'G', 'T'}
                                              for nt in seqstr))
        subs = [f"{ref_nt}{site}{seq_nt}" for site, (ref_nt, seq_nt)
                in enumerate(zip(ref_seqstr, seqstr), start=1)
                if ref_nt != seq_nt and seq_nt in ['A', 'C', 'G', 'T']]
        data['n_subs_to_ref'].append(len(subs))
        assert length == sum(data[key][-1] for key in ['n_gapped_to_ref',
                                                       'n_ident_to_ref',
                                                       'n_subs_to_ref',
                                                       'n_ambiguous_to_ref'])
        data['substitutions'].append(','.join(subs))
        head_entries = seq.description.split()
        data['strain'].append(head_entries[0])
        for entry in head_entries[1:]:
            assert entry.count('=') == 1
            key, val = entry.split('=')
            data[key].append(val)
    iseq += 1
print(f"Read {iseq} sequences from the alignment including the reference")

assert all(len(val) == iseq - 1 for val in data.values())
df = pd.DataFrame(data)
assert len(df) == iseq - 1

print(f"Writing data to {snakemake.output.csv}")
(df
 .to_csv(snakemake.output.csv, index=False)
 )

"""Implements ``snakemake`` rule `add_outgroup_to_alignment`."""


import Bio.SeqIO

import pandas as pd


outgroup = snakemake.wildcards.outgroup
start = snakemake.params.region['start']
end = snakemake.params.region['end']
assert start <= end

print(f"Getting outgroup {outgroup} from {start} to {end} from {snakemake.input.comparator_map}")
comparator_map = (
    pd.read_csv(snakemake.input.comparator_map)
    .set_index('site')
    [outgroup]
    .to_dict()
    )
outgroup_seq = ''.join(comparator_map[site] for site in range(start, end + 1))
outgroup_seqrecord = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(outgroup_seq),
                                             id=outgroup,
                                             name='',
                                             description='',
                                             )

alignment = [outgroup_seqrecord]
print(f"Reading alignment without outgroup from {snakemake.input.alignment}")
for seq in Bio.SeqIO.parse(snakemake.input.alignment, 'fasta'):
    assert outgroup != seq.id
    alignment.append(seq)
    assert len(outgroup_seqrecord) == len(alignment[-1]) == end - start + 1, (
            f"{len(outgroup_seqrecord)=}\n{len(alignment[-1])=}\n{end - start + 1=}")

print(f"Writing the {len(alignment)} sequences including outgroup to {snakemake.output.alignment}")
Bio.SeqIO.write(alignment, snakemake.output.alignment, 'fasta')

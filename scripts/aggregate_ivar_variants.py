"""Aggregate variants called by ``ivar``."""


import pandas as pd


assert len(snakemake.input.tsvs) == len(snakemake.params.descriptors)

id_to_gene = {}
for genome, ref_gff in zip(snakemake.params.genomes, snakemake.input.ref_gffs):
    id_to_gene[genome] = {}
    with open(ref_gff) as f:
        for line in f:
            if line[0] != '#':
                entry = line.split('\t')[8]
                key = val = None
                for prop_val in entry.split(';'):
                    if prop_val.split('=')[0] == 'ID':
                        key = prop_val.split('=')[1]
                    elif prop_val.split('=')[0] == 'gene':
                        val = prop_val.split('=')[1]
                if key and val:
                    id_to_gene[genome][key] = val

df = pd.concat([pd.read_csv(f, sep='\t')
                  .assign(**d)
                  .assign(gene=lambda x: x['GFF_FEATURE'].map(id_to_gene[d['genome']]))
                for f, d in zip(snakemake.input.tsvs,
                                snakemake.params.descriptors)
                ],
               ignore_index=True,
               )

df.to_csv(snakemake.output.agg_csv,
          index=False)

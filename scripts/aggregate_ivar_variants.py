"""Aggregate variants called by ``ivar``."""


import pandas as pd


assert len(snakemake.input) == len(snakemake.params.descriptors)

df = pd.concat([pd.read_csv(f, sep='\t').assign(**d)
                for f, d in zip(snakemake.input,
                                snakemake.params.descriptors)
                ],
               ignore_index=True,
               )

df.to_csv(snakemake.output.agg_csv,
          index=False)

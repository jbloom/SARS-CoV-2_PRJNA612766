import pandas as pd

all_df = pd.read_csv(snakemake.input.all_csv, na_filter=False)

outgroups = snakemake.params.outgroups

outgroup_dist = (
    all_df
    .melt(id_vars=['substitutions', 'representative_strain'],
          value_vars=[f"{outgroup}_delta_dist" for outgroup in outgroups],
          var_name='outgroup',
          value_name='distance',
          )
    .assign(outgroup=lambda x: x['outgroup'].str.replace('_delta_dist', ''),
            min_distance=lambda x: x.groupby('outgroup')['distance'].transform('min'),
            )
    .query('distance == min_distance')
    )

print(outgroup_dist)

# **assume** (and then check) that all outgroups have same possible progenitors
possible_progenitors = outgroup_dist['representative_strain'].unique().tolist()
print(f"The possible progenitors are:\n" + '\n'.join(possible_progenitors))

assert all(set(possible_progenitors) == set(outgroup_dist
                                            .query('outgroup == @outgroup')
                                            ['representative_strain'])
           for outgroup in outgroups)

print(f"Writing to {snakemake.output.progenitors}")
with open(snakemake.output.progenitors, 'w') as f:
    f.write('\n'.join(possible_progenitors))

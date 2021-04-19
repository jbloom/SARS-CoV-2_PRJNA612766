import pandas as pd

df_list = []
cols = None
for f in snakemake.input:
    df = pd.read_csv(f)
    if cols is None:
        cols = set(df.columns)
    elif cols != set(df.columns):
        raise ValueError('data frames do not all have same columns')
    df_list.append(df)

merged_df = pd.concat(df_list, sort=False, ignore_index=True)
if len(merged_df) != len(merged_df.drop_duplicates()):
    raise ValueError('duplicate columns in merged data frame')

merged_df.to_csv(snakemake.output.merged_csv, index=False)

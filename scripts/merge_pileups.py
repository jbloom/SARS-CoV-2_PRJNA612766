import pandas as pd

pileup_cols = ['site', 'ref_nt', 'A', 'C', 'G', 'T']
codes = {'aligner': {}, 'genome': {}, 'sample': {}}

df = pd.DataFrame()
for csv_file, descriptor in zip(snakemake.input.csvs,
                                snakemake.params.descriptors):
    # get unique sequential numerical codes for each aligner/genome/sample
    csv_codes = {}
    for code_type, code_d in codes.items():
        if descriptor[code_type] not in code_d:
            if len(code_d):
                prior_code = max(code_d.values())
            else:
                prior_code = 0
            code_d[descriptor[code_type]] = prior_code + 1
        csv_codes[code_type] = code_d[descriptor[code_type]]

    # read data, keep columns of interest, and add codes
    df = df.append(pd.read_csv(csv_file)
                   [pileup_cols]
                   .assign(**csv_codes)
                   )

df.to_csv(snakemake.output.pileup, index=False)
for code_type, code_d in codes.items():
    outfile = getattr(snakemake.output, f"{code_type}_key")
    (pd.Series(code_d)
     .rename_axis(f"{code_type}_name")
     .rename(code_type)
     .to_csv(outfile)
     )
     

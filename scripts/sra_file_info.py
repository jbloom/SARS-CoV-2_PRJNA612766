"""Implements ``snakemake`` rule `sra_file_info`."""


import os
import re
import subprocess

import pandas as pd


sra_files = snakemake.input.sra_files
accessions = [os.path.splitext(os.path.basename(sra_file))[0]
              for sra_file in sra_files]

print(f"Getting information for {len(sra_files)} SRA files.")

records = []
for acc, f in zip(accessions, sra_files):

    obj_ts = subprocess.run(['vdb-dump', '--obj_timestamp', f],
                            capture_output=True,
                            text=True).stdout.strip()

    info = subprocess.run(['vdb-dump', '--info', f],
                          capture_output=True,
                          text=True).stdout
    info_m = re.search('TIME\s+\:\s\w+\s\((\d{2}/\d{2}/\d{4} \d{2}\:\d{2})\)',
                       info)
    if not info_m:
        raise ValueError(f"did not match info:\n{info}")
    info_date = info_m.group(1)

    records.append((acc, f, obj_ts, info_date))

print(f"Writing results to {snakemake.output.csv}")
(pd.DataFrame(records,
              columns=['accession', 'sra_file', 'object_timestamp',
                       'info_timestamp'])
 .assign(object_timestamp=lambda x: pd.to_datetime(x['object_timestamp']),
         info_timestamp=lambda x: pd.to_datetime(x['info_timestamp']),
         )
 .to_csv(snakemake.output.csv, index=False)
 )


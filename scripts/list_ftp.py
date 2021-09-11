"""Implements ``snakemake`` rule `list_ftp`."""

ftp_host = snakemake.params.host
ftp_path = snakemake.params.path
ftp_user = snakemake.params.user
ftp_password = snakemake.params.password
output_files = snakemake.output.files
output_dirs = snakemake.output.dirs


import ftputil


dir_list = []
file_list = []
print(f"Connecting to {ftp_host=} as {ftp_user=} with {ftp_password=}")
with ftputil.FTPHost(ftp_host, ftp_user, ftp_password, timeout=480) as ftp:
    ftp.encoding = 'utf-8'
    print(f"Listing contents of {ftp_path=}")
    path_list = ftp.listdir(ftp_path)
    print(f"Found {len(path_list)=} entries")
    for i, name in enumerate(path_list, 1):
        full_name = f"{ftp_path}/{name}"
        if ftp.path.isdir(full_name):
            print(f"{full_name=} is directory")
            dir_list.append(full_name)
        else:
            file_list.append(full_name)
            print(f"{full_name=} is file")

print(f"Writing {len(file_list)=} files to {output_files}")
with open(output_files, 'w') as f:
    f.write('\n'.join(file_list))

print(f"Writing {len(dir_list)=} files to {output_dirs}")
with open(output_dirs, 'w') as f:
    f.write('\n'.join(dir_list))

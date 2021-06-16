import glob
import os


for f in glob.glob('*%*'):
    new_f = f.replace('%', '-')
    print(f"Renaming {f} to {new_f}")
    os.rename(f, new_f)

import sys
import glob
import os

path  = "/data/local/savin/lwa/lwafiles/repeaters/DD002_049*"

dirlist = sorted(glob.glob(path))

for dir in dirlist:
    print(dir)
    os.chdir(dir)
    os.mkdir(f"{os.path.basename(dir)}_Images")
    os.system(f"mv *.png {os.path.basename(dir)}_Images")
    os.system(f"mv Images* {os.path.basename(dir)}_Images")
    os.system(f"scp -r {os.path.basename(dir)}_Images  savin@hercules.phys.unm.edu:/home/savin/frb_lwa_chime/repeaters/")

#!/usr/bin/env python3

import os
import glob


path = '/data/local/savin/lwa/lwafiles/DD002_90[3,4]*'

dirs = sorted(glob.glob(path))

print(dirs)

metadir = '/data/network/recent_data/savin/'

for dir_path in dirs:
    
    print(f"Going to {dir_path}")
    os.chdir(dir_path)
    metabase = os.path.basename(dir_path)
    metafile = metadir+metabase+'.tgz'
    files = sorted(glob.glob(dir_path+'/*.fits'))
    if len(files) == 3:
       os.system(f"python3  /data/local/savin/lwa/lwa_scripts/search_fits.py  -m {metafile} {files[0]} {files[1]} {files[2]}")

    elif len(files) == 2:
       os.system(f"python3  /data/local/savin/lwa/lwa_scripts/search_fits.py -m {metafile} {files[0]} {files[1]}")
     
    elif len(files) == 1:
       os.system(f"python3  /data/local/savin/lwa/lwa_scripts/search_fits.py -m {metafile} {files[0]} ") 
     
    else:
       print(f"No fits files to process in {dir_path}")
       continue

    


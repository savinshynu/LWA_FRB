#!/usr/bin/env python3

import os
import git
import sys
import argparse
import time
import glob

from astropy.coordinates import Angle as AstroAngle

from lsl.reader.ldp import DRSpecFile
from lsl.common import metabundle, metabundleADP
from lsl.misc import parser as aph

from detailed_info import load_repeaters


ANALYSIS_PATH = os.path.abspath(os.path.dirname(__file__))

OUT_PATH = "/data/local/savin/lwa/lwafiles/repeaters"

def main(filename):
    # Figure out our revision
    try:
        repo = git.Repo(os.path.dirname(os.path.abspath(__file__)))
        try:
            branch = repo.active_branch.name
            hexsha = repo.active_branch.commit.hexsha
        except TypeError:
            branch = '<detached>'
            hexsha = repo.head.commit.hexsha
        shortsha = hexsha[-7:]
        dirty = ' (dirty)' if repo.is_dirty() else ''
    except git.exc.GitError:
        branch = 'unknown'
        hexsha = 'unknown'
        shortsha = 'unknown'
        dirty = ''
        
    # Get the input filename
    #filename = args.filename
    
    try:
        ## Is it a raw data file?
        idf = DRSpecFile(filename)
        
        # Extract what is needed
        if args.source is None:
            args.source = "FRB"
        if args.ra is None:
            args.ra = "00:00:00.00"
        if args.dec is None:
            args.dec = "+00:00:00.0"
        ra = str(args.ra)
        dec = str(args.dec)
        if args.dm is None:
            raise RuntimeError("A DM is required when passing in a DRX file")
        name = f"{args.source}_DM{args.dm:.1f}_drx"
        mjd = int(idf.start_time.mjd)
        beam = int(idf.beam)
        drsfile = os.path.abspath(args.filename)
        
        idf.close()
        
    except:
        ## Treat it as metadata
        try:
            sdf = metabundle.get_sdf(filename)
            datafile = metabundle.get_session_metadata(filename)
        except RuntimeError:
            sdf = metabundleADP.get_sdf(filename)
            datafile = metabundleADP.get_session_metadata(filename)
            
        # Extract what we need
        name = sdf.sessions[0].observations[0].target
        mjd = sdf.sessions[0].observations[0].mjd
        ra = sdf.sessions[0].observations[0].ra
        dec = sdf.sessions[0].observations[0].dec
        beam = sdf.sessions[0].drx_beam
        assert(sdf.sessions[0].spectrometer_channels != 0)
        drsfile = os.path.abspath(os.path.dirname(filename))
        drsfile = os.path.join(drsfile, datafile[1]['tag'])
        
    # The DM is stored in the target list
    repeaters = load_repeaters()
    dm = repeaters[name]

    #dm_str = str(dm).zfill(3)

    # The "short name" is just FRB+the DM
    #sname = 'FRB%3.0f' % dm
    #sname = 'FRB'+dm_str

    sname = name

    # String versions of the coordinates
    ra_str = AstroAngle(ra, unit='hourangle')
    ra_str = str(ra_str).replace('h', ':').replace('m', ':').replace('s', '')
    dec_str = AstroAngle(dec, unit='deg')
    dec_str = str(dec_str).replace('d', ':').replace('m', ':').replace('s', '')
    
    print (f"DM: {dm} , RA: {ra_str}, Dec : {dec_str}, Source: {sname}, DRSFile: {drsfile}")

    DIR_PATH = os.path.join(OUT_PATH,os.path.splitext(os.path.basename(filename))[0])

    print (DIR_PATH)
    try:
       os.mkdir(DIR_PATH)
    except OSError:
       print("Directory alrealdy exists")

    os.chdir(DIR_PATH)
    t0 = time.time()

    print(drsfile)

    os.system(f"/usr/local/extensions/Pulsar/writePsrfits2FromDRSpec.py --source={sname} --ra={ra_str} --dec={dec_str} --nsblk=128 {drsfile} > writepsrfits.out")
    os.system("python /data/local/savin/lwa/lwa_scripts/check_writepsrfits_result.py writepsrfits.out")
    #os.system(f"combine_lwa2 -o drx_{mjd}_{sname}_b{beam} drx_{mjd}_{sname}_b{beam}t2_0001.fits drx_{mjd}_{sname}_b{beam}t1_0001.fits")

    t1 = time.time()

    print(f"Time taken for processing: {t1-t0} s")



if __name__ == '__main__':
    
    
    path = "/data/network/recent_data/savin/repeaters/DD002_0[5,6]*.tgz"
    files = sorted(glob.glob(path))
    print (files)
    for filename in files:
        main(filename)
        print(filename)


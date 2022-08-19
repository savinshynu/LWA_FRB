#!/usr/bin/env python3

import os
import glob
import git
import sys
import argparse
#import subprocess
import time
from astropy.coordinates import Angle as AstroAngle

from lsl.reader.ldp import DRXFile
from lsl.common import metabundle, metabundleADP
from lsl.misc import parser as aph


ANALYSIS_PATH = os.path.abspath(os.path.dirname(__file__))

OUT_PATH = "/data/local/savin/lwa/lwafiles/"

def main(filename):
        
    
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
    assert(sdf.sessions[0].spectrometer_channels == 0)
    drxfile = os.path.abspath(os.path.dirname(filename))
    drxfile = os.path.join(drxfile, datafile[1]['tag'])
        
    # The DM is currently encoded in the target name
    dm = name.split('_DM', 1)[1]
    dm = dm.split('_', 1)[0]
    dm = float(dm)
    
    # The "short name" is just FRB+the DM
    sname = 'FRB%3.0f' % dm

    # String versions of the coordinates
    ra_str = AstroAngle(ra, unit='hourangle')
    ra_str = str(ra_str).replace('h', ':').replace('m', ':').replace('s', '')
    dec_str = AstroAngle(dec, unit='deg')
    dec_str = str(dec_str).replace('d', ':').replace('m', ':').replace('s', '')
   
    print (f"DM: {dm} , RA: {ra_str}, Dec : {dec_str}, Source: {sname}, DRXFile: {drxfile}")
    
    DIR_PATH = os.path.join(OUT_PATH,os.path.splitext(os.path.basename(filename))[0])
    print (DIR_PATH)
    try:
       os.mkdir(DIR_PATH)
    except OSError:
       print("Directory alrealdy exists")
    
    os.chdir(DIR_PATH)
    t0 = time.time()
    
    
    os.system(f"/usr/local/extensions/Pulsar/writePsrfits2DMulti.py --source={sname} --ra={ra_str} --dec={dec_str} --nchan=4096 --nsblk=16384 --yes {dm} {drxfile} > writepsrfits.out")

    #os.system(f"/usr/local/extensions/Pulsar/writePsrfits2D.py --source={sname} --ra={ra_str} --dec={dec_str} --nchan=4096 --nsblk=16384  {dm} {drxfile} > writepsrfits.out")


    os.system("python /home/pulsar/bin/check_writepsrfits_result.py writepsrfits.out")


    os.system(f"combine_lwa2 -o drx_{mjd}_{sname}_b{beam} drx_{mjd}_{sname}_b{beam}t2_0001.fits drx_{mjd}_{sname}_b{beam}t1_0001.fits")

    t1 = time.time()

    print(f"Time taken for conversion : {t1-t0}")

if __name__ == '__main__':
    """
    parser = argparse.ArgumentParser(
        description='build a Makefile for analyzing FRB data from the LWA',
        epilog='NOTE: --source, --ra, --dec, and --dm are only used if a DRX file is given, otherwise the values are automatically extracted from the metadata.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filename', type=str, 
                        help='filename process - metadata or DRX')
    parser.add_argument('-s', '--source', type=str, 
                        help='source name')
    parser.add_argument('-r', '--ra', type=aph.hours, 
                        help='right ascension; HH:MM:SS.SS, J2000')
    parser.add_argument('-d', '--dec', type=aph.degrees, 
                        help='declination; sDD:MM:SS.S, J2000')
    parser.add_argument('-e', '--dm', type=float,
                        help='dispersion measure; pc/cm^3')
    args = parser.parse_args()
    """
    path = "/data/network/recent_data/savin/DD002_901[7,9]*.tgz"
    files = sorted(glob.glob(path))
    print (files)
    for filename in files:
        main(filename)
        print(filename)

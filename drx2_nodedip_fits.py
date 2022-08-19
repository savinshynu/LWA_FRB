#!/usr/bin/env python3

import os
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

OUT_PATH = "/data/local/savin/lwafiles/"

def main(args):
        
    # Get the input filename
    filename = args.filename
    
    try:
        ## Is it a raw data file?
        idf = DRXFile(args.filename)
        
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
        drxfile = os.path.abspath(args.filename)
        
        idf.close()
        
    except:
        ## Treat it as metadata
        try:
            sdf = metabundle.get_sdf(args.filename)
            datafile = metabundle.get_session_metadata(args.filename)
        except RuntimeError:
            sdf = metabundleADP.get_sdf(args.filename)
            datafile = metabundleADP.get_session_metadata(args.filename)
            
        # Extract what we need
        name = sdf.sessions[0].observations[0].target
        mjd = sdf.sessions[0].observations[0].mjd
        ra = sdf.sessions[0].observations[0].ra
        dec = sdf.sessions[0].observations[0].dec
        beam = sdf.sessions[0].drx_beam
        assert(sdf.sessions[0].spectrometer_channels == 0)
        drxfile = os.path.abspath(os.path.dirname(args.filename))
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
   
    print (dm, ra_str, dec_str, sname, drxfile)
    
    
    DIR_PATH = os.path.join(OUT_PATH,os.path.splitext(os.path.basename(filename))[0])
    print (DIR_PATH)
    try:
       os.mkdir(DIR_PATH)
    except OSError:
       print("Directory alrealdy exists")
    
    os.chdir(DIR_PATH)
    print (time.time())
    print (dm)
    
    os.system(f"/usr/local/extensions/Pulsar/writePsrfits2Multi.py --source={sname} --ra={ra_str} --dec={dec_str} --nchan=4096 --nsblk=16384 --yes  {drxfile} > writepsrfits.out")
    
    print (time.time()) 

    
    os.system("python /home/pulsar/bin/check_writepsrfits_result.py writepsrfits.out")


    os.system(f"combine_lwa2 -o drx_{mjd}_{sname}_b{beam} drx_{mjd}_{sname}_b{beam}t2_0001.fits drx_{mjd}_{sname}_b{beam}t1_0001.fits")

    """
    python /home/pulsar/bin/create_ignore_chan.py $< > $*.ignorechan

    rfifind -time 30 -o $* $< -ignorechan $*.ignorechan > $*_rfifind.out

    python2.7 /home/pulsar/bin/rfifind.py $^

    python /home/pulsar/bin/crab_plot.py $*_scatteredsearch/*.scatteredsearch



    prepsubband -lodm {dm-25} -numdms 51 -dmstep 1 -nsub 512 -ignorechan $*.ignorechan -o $*_DM{dm}_scatteredsearch/$* $*_0001.fits -mask $*_rfifind.mask
    python3 {ANALYSIS_PATH}/multi_crab_search.py -n 3 $*_DM{dm}_scatteredsearch/*.dat
    rm $*_DM{dm}_scatteredsearch/*.dat
    """


if __name__ == '__main__':
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
    main(args)

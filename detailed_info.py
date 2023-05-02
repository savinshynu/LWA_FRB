#!/usr/bin/env python3
"""
Script to give details about a metadata tarball from a CHIME triggered follow-up
observation.
"""

import os
import sys
import ephem
import numpy
from datetime import datetime

from astropy.coordinates import Angle as AstroAngle

from lsl.common.stations import lwa1, lwasv
from lsl.common.mcs import mjdmpm_to_datetime
from lsl.common import metabundle, metabundleADP
from lsl.reader.drx import FILTER_CODES
from lsl.misc.dedispersion import delay


def load_repeaters():
    """
    Load in the table of repeaters that we monitor in case this metadata is from
    one of those observations.
    """
    
    filename = os.path.abspath(os.path.dirname(__file__))
    filename = os.path.join(filename, 'DD002_Targets.txt')
    
    data = {}
    with open(filename, 'r') as fh:
        for line in fh:
            if len(line) < 3:
                continue
            if line[0] == '#':
                continue
                
            name, ra, dec, dm, _, _, _ = line.split(None, 6)
            data[name] = float(dm)
    return data


def main(args):
    for filename in args:
        # Load it in
        try:
            sdf = metabundle.get_sdf(filename)
            datafile = metabundle.get_session_metadata(filename)
            site = lwa1
        except RuntimeError:
            sdf = metabundleADP.get_sdf(filename)
            datafile = metabundleADP.get_session_metadata(filename)
            site = lwasv
            
        # Extract what we need
        name = sdf.sessions[0].observations[0].target
        mjd = sdf.sessions[0].observations[0].mjd
        mpm = sdf.sessions[0].observations[0].mpm
        dur = sdf.sessions[0].observations[0].dur
        ra = sdf.sessions[0].observations[0].ra
        dec = sdf.sessions[0].observations[0].dec
        freqL = sdf.sessions[0].observations[0].frequency1
        freqL = min([freqL, sdf.sessions[0].observations[0].frequency2])
        freqU = sdf.sessions[0].observations[0].frequency1
        freqU = max([freqU, sdf.sessions[0].observations[0].frequency2])
        bw = FILTER_CODES[sdf.sessions[0].observations[0].filter]
        datafile = datafile[1]['tag']
        
        # Parse and process
        try:
            ## Trigger
            event, dm, _ = name.split('_', 2)
            event = event.split('#', 1)[1]
            event = datetime.strptime(event, "%Y-%m-%d-%H:%M:%S.%fUTC")
            dm = float(dm.replace('DM', ''))
            is_trigger = True
        except ValueError:
            ## Repeater
            repeaters = load_repeaters()
            try:
                dm = repeaters[name]
                event = name
                is_trigger = False
            except KeyError:
                print(f"ERROR: cannot find information on {name}, skipping {os.path.basename(filename)}")
                continue
                
        start = mjdmpm_to_datetime(mjd, mpm)
        ra_str = AstroAngle(ra, unit='hourangle')
        ra_str = str(ra_str)
        dec_str = AstroAngle(dec, unit='deg')
        dec_str = str(dec_str)
        
        # Calculate the elevation
        bdy = ephem.FixedBody()
        bdy._ra = ra*numpy.pi/12
        bdy._dec = dec*numpy.pi/180
        bdy._epoch = ephem.J2000
        mid = mjdmpm_to_datetime(mjd, mpm+dur//2)
        site.date = mid.strftime("%Y-%m-%d %H:%M:%S")
        site.compute(bdy)
        mid_el = bdy.alt*180/numpy.pi
        
        # Calculate the lags
        if is_trigger:
            trigger_lag = (start-event).total_seconds()
            disp_lagU = delay([freqU+bw/2.0, 400e6], dm)[0]
            disp_lagL = delay([freqL-bw/2.0, 400e6], dm)[0]
            is_repeater = False
            
        # Report
        print(f"{os.path.basename(filename)}:")
        print(f"  Event with DM {dm:.1f} pc/cm^3 at {event}")
        print(f"    Data file is {datafile} from {site.name}")
        if is_trigger:
            print(f"    Observation started {trigger_lag:.3f} s after and went for {dur/1000:.1f} s")
            print(f"    Delay from 400 MHz down to {(freqU+bw/2.0)/1e6:.1f} MHz is {disp_lagU:.3f} s")
            print(f"      Pulse should reach the top of the band {(disp_lagU-trigger_lag):.3f} s after the start")
            print(f"      Pulse should leave the bottom of the band {(dur/1000-disp_lagL+trigger_lag):.3f} s before the end")
        print(f"  Position:")
        print(f"    RA: {ra_str}")
        print(f"    Dec: {dec_str}")
        print(f"    Alt: {mid_el:.1f} deg")

    
if __name__ == "__main__":
    main(sys.argv[1:])

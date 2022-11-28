"""
Modified version of writeHDF5fromPSRfits by J.Dowell for collecting data from PSRfits files
Further bandpass calibration, rfi removal, incoherent dedispersion, collecting candidates and plotting candidates
"""

import sys
import os
import re
import time
import argparse
import numpy as np
from astropy.time import Time as AstroTime
from astropy.io import fits as astrofits
from matplotlib import pyplot as plt
from utils import comp_bp, flag_rfi_time, boxcar_search, get_cand_list, plot_events, matchFilter_search, calc_triggerlag  
from dedispersion import incoherent, delay
import matplotlib
import tqdm as tq
matplotlib.use('Agg')


_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')


#Edit this list each time to get the slices


time_windows = [949290, 1258727, 1336750, 1632407, 1632408, 1664010, ]


def main(args):
    # Parse command line options
    filenames = args.filenames

    #opening a directory to save image candidates

    dirname = os.path.dirname(filenames[0])
    
    outdir = os.path.join(dirname, 'data_slices')
    try:
        os.mkdir(outdir)
    except OSError:
        print (f"{outdir} directory already exists")
       


    # Open the file and load in basic information about the observation's goal
    for c,filename in enumerate(filenames):
        print(filename)

        ## Ready the PSRFITS file
        hdulist = astrofits.open(filename, memmap=True)
        
        ## Try to find out the beam/tuning
        mtch = _fnRE.search(filename)

        if mtch is None:
            beam = 0
            tuning = 1
        else:
            beam = int(mtch.group('beam'))
            try:
                tuning = int(mtch.group('tuning'))
            except:
                tuning = 1

        ## File specifics
        sourceName = hdulist[0].header['SRC_NAME']
        ra = hdulist[0].header['RA']
        ra = ra.split(':', 2)
        ra = sum([float(v)/60**i for i,v in enumerate(ra)])*15.0
        dec = hdulist[0].header['DEC']
        decSign = -1.0 if dec.find('-') != -1 else 1.0
        dec = dec.replace('-', '')
        dec = dec.split(':', 2)
        dec = decSign*sum([float(v)/60**i for i,v in enumerate(dec)])
        epoch = float(hdulist[0].header['EQUINOX'])
        
        dm = hdulist[0].header['CHAN_DM']

        tStart = AstroTime(hdulist[0].header['STT_IMJD'], (hdulist[0].header['STT_SMJD'] + hdulist[0].header['STT_OFFS'])/86400.0,
                           format='mjd', scale='utc')
        cFreq = hdulist[0].header['OBSFREQ']*1e6        # MHz -> Hz
        srate = hdulist[0].header['OBSBW']*1e6          # MHz -> Hz
        LFFT = hdulist[1].header['NCHAN']
        tInt = hdulist[1].header['TBIN']
        nSubs = hdulist[1].header['NSBLK']
        tSubs = nSubs*tInt
        nPol = hdulist[1].header['NPOL']
        if nPol == 1:
            data_products = ['I',]
        elif nPol == 2:
            if hdulist[0].header['FD_POLN'] == 'CIRC':
                data_products = ['LL', 'RR']
            else:
                data_products = ['XX', 'YY']
        else:
            data_products = ['I', 'Q', 'U', 'V']
        nChunks = len(hdulist[1].data)
        


        ## Pre-process the start time
        tStartI = int(tStart.unix)
        tStartF = tStart.unix - tStartI


        ## Report
        print("Filename: %s (%i of %i)" % (filename, c+1, len(filenames)))
        print("Date of First Frame: %s" % tStart.datetime)
        print("Beam: %i" % beam)
        print("Tuning: %i" % tuning)
        print("Sample Rate: %i Hz" % srate)
        print("Tuning Frequency: %.3f Hz" % cFreq)
        print("---")
        print("Target: %s" % sourceName)
        print("RA: %.3f hours" % (ra/15.0,))
        print("Dec: %.3f degrees" % dec)
        print("DM: %0.3f " % dm)
        print("Data Products: %s" % ','.join(data_products))
        print("Integration Time: %.3f ms" % (tInt*1e3,))
        print("Number of sub blocks: %f " % (nSubs))
        print("Sub Integration time: %f" % (tSubs))
        print("Sub-integrations: %i (%.3f s)" % (nChunks, nChunks*tSubs))
        print("---")
        #print("Offset: %.3f s (%i subints.)" % (args.skip, skip))
        #print("Duration: %.3f s (%i subints.)" % (args.duration, dur))
        print("Transform Length: %i" % LFFT)
        

        ## Frequency information
        freq = hdulist[1].data[0][12]*1e6               # MHz -> Hz
        
        if  args.metafile:
            tdelay = delay([np.max(freq), 400e6], dm)[0]
            trigger_lag  = calc_triggerlag(args.metafile)

            exp_time = tdelay - trigger_lag
            exp_samp = int(exp_time/tInt)

            print(f"Expected time of pulse arrival  is {exp_time}")



        skip = 0 # Number of sub integrations to skip in starting 
        dur =  nChunks # Number of sub integrations to load the data
        

        
        #Initialize the big array to capture data from different sub integrations
        wfdata = np.zeros((dur*nSubs, LFFT, len(data_products)), dtype = np.float32)
        tdat = np.arange(0,dur*nSubs)*tInt


        print("Collecting data")
        t0 = time.time()
        ## Read in the data and apply what ever scaling is needed
        for i in tq.tqdm(range(skip, skip+dur)):
            ### Access the correct subintegration
            subint = hdulist[1].data[i]

            ### Pull out various bits that we need, including:
            ###  * the start time of the subint. - tOff
            ###  * the weight mask, converted to binary - msk
            ###  * the scale and offset values - bscl and bzero
            ###  * the actual data - data
            tOff = subint[1] - subint[0]/2
            msk = np.where(subint[13] >= 0.5, False, True)
            bzero = subint[14]
            bscl = subint[15]
            bzero.shape = (LFFT,nPol)
            bscl.shape = (LFFT,nPol)
            bzero = bzero.T
            bscl = bscl.T
            data = subint[16]
            data.shape = (nSubs,LFFT,nPol)
            data = data.T
            #print(data.shape)
            ### Apply the scaling/offset to the data and save the results 
            ### to an array
            for j in range(nSubs):
                k = (i-skip)*nSubs + j
                d = data[:,:,j]*bscl + bzero
                
                
                for l,p in enumerate(data_products):
                    wfdata[k,:,l] = d[l,:]


        t1 = time.time()
        print(f"Time taken to collect data: {t1-t0} s")

        #Assuming data_products == 'I':
        
        if len(data_products) == 1: 
            wfdata = wfdata[:,:,0] 
        else:
            del wfdata
            sys.exit("More than one polarization product: check the data")
            
        #computing the bandpass    
        bp = comp_bp(wfdata,nSubs)
        
        #Correcting the bandpass
        wfdata /= bp
        
        print("Finished bandpass correction")

        #Flagging bad rfi and time chunks
        wfdata = flag_rfi_time(wfdata,nSubs)

        print("Finished RFI Flagging")

        #Apply incoherent dedispersion
        wfdata = incoherent(freq, wfdata, tInt, dm, boundary='wrap', fill_value=np.nan)
        
        t2 = time.time()
        print(f"Bandpass, flagging, dedispersion  {t2-t1} s")
        
        maxind = wfdata.shape[0]

        print('Starting data slicing')
        
        for twin in time_windows:
            slicename = filename+'_'+str(twin)
            slicename = os.path.join(outdir, slicename)
            print(slicename)
            np.savez(slicename, wfdat = wfdata[max(0,twin-5000):min(twin+5000,maxind),:], tdat = tdat[max(0,twin-5000):min(twin+5000, maxind)], freq = freq, dm = dm, ra = ra, dec = dec, source = sourceName, cfreq = cFreq, srate = srate, lfft = LFFT, tint = tInt)


        t3 = time.time()
        print (f"Data slicing took {t3-t2} s")
        del wfdata





        hdulist.close()
        print(f"Finished processing {filename}")
    


if  __name__ == "__main__": 


    
    parser = argparse.ArgumentParser(
        description='Processing and looking for single pulses in FRB FITS files',
        epilog='NOTE: along with FITS file, add the metafile also',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filenames', type=str, 
                        help='FITS filename to process',  nargs = '+')
    parser.add_argument('-m', '--metafile', type=str, 
                        help='Name of the observation metafile')
    args = parser.parse_args()




    main(args)


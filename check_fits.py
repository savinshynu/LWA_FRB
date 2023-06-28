
"""
Reads a fits file and make some basic plots and prints out the metadata
"""

import sys
import os
import re
import time
import tqdm as tq
import argparse
import numpy as np
from astropy.time import Time as AstroTime
from astropy.io import fits as astrofits
from matplotlib import pyplot as plt
#from utils import comp_bp, flag_rfi_time, boxcar_search, get_cand_list, plot_events, matchFilter_search, calc_triggerlag  
#from dedispersion import incoherent, delay
import matplotlib
matplotlib.use('Tkagg') #'Agg'
from mad import median_absolute_deviation as mad

_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')

def main(args):
    # Parse command line options
    filenames = args.filenames


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
        
        if len(data_products) != 1:
            print(f"data products: {data_products}")
            sys.exit("More than one polarization product: check the data") 

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
        

        skip = 0 # Number of sub integrations to skip in starting 
        dur =  nChunks # Number of sub integrations to load the data
        

        
        #Initialize the big array to capture data from different sub integrations
        mad_ar = np.zeros(dur, dtype = np.float32)
        wfdata = np.zeros((nSubs, LFFT, len(data_products)), dtype = np.float32)

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
                    wfdata[j,:,l] = d[l,:]
            
            mad_ar[i] = mad(wfdata, axis = None)


 
        hdulist.close()
        print(f"Finished processing {filename}")
    
        plt.plot(mad_ar)
        plt.xlabel("No. of Sub-integrations")
        plt.ylabel("MAD")
        plt.show() 

if  __name__ == "__main__": 


    
    parser = argparse.ArgumentParser(
        description='Processing and looking for single pulses in FRB FITS files',
        epilog='NOTE: along with FITS file, add the metafile also',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('filenames', type=str, 
                        help='FITS filename to process',  nargs = '+')
    args = parser.parse_args()




    main(args)


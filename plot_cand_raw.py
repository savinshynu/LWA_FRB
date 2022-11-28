"""
Based on writeHDF5fromPSRfits by J.Dowell

"""

import sys
import time
import numpy as np
from astropy.time import Time as AstroTime
from astropy.io import fits as astrofits
from matplotlib import pyplot as plt
from utils import comp_bp, flag_rfi_time, boxcar_search, get_cand_list 
from dedispersion import incoherent

_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')
def main():
    # Parse command line options
    filenames = [sys.argv[1]]
    #sample
    skip = int(sys.argv[2])
    dur = int(sys.argv[3])
    bin = [int(sys.argv[4])]

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

        tStart = AstroTime(hdulist[0].header['STT_IMJD'], (hdulist[0].header['STT_SMJD'] + hdulist[0].header['STT_OFFS'])/86400.0,
                           format='mjd', scale='utc')
        cFreq = hdulist[0].header['OBSFREQ']*1e6        # MHz -> Hz
        srate = hdulist[0].header['OBSBW']*1e6          # MHz -> Hz
        LFFT = hdulist[1].header['NCHAN']
        tInt = hdulist[1].header['TBIN']
        nSubs = hdulist[1].header['NSBLK']
        tSubs = nSubs*tInt
        dm = hdulist[0].header['CHAN_DM']

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
        ## File cross-validation to make sure that everything is in order
        try:
            validationPass = True
            for keyword in ('sourceName', 'ra', 'dec', 'epoch', 'tStart', 'srate', 'LFFT', 'tInt', 'tSubs', 'nPol', 'nChunks'):
                keywordOld = keyword+"Old"
                valid = eval("%s == %s" % (keyword, keywordOld))
                if not valid:
                    print("ERROR:  Detail '%s' of %s does not match that of the first file" % (keyword, os.path.basename(filename)))
                    print("ERROR:  Aborting")
                    validationPass = False

            if not validationPass:
                continue

        except NameError as e:
            sourceNameOld = sourceName
            raOld = ra
            decOld = dec
            epochOld = epoch
            tStartOld = tStart
            srateOld = srate
            LFFTOld = LFFT
            tIntOld = tInt
            tSubsOld = tSubs
            nPolOld = nPol
            nChunksOld = nChunks
        ## Pre-process the start time
        tStartI = int(tStart.unix)
        tStartF = tStart.unix - tStartI

        """
        ## Convert the skip and duration values to subblocks
        skip = int(round(args.skip / tSubs))
        dur  = int(round(args.duration / tSubs))
        dur  = dur if dur else 1
        args.skip = skip * tSubs
        args.duration = dur * tSubs
        """

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
        

        #skip = 0
        #dur = nChunks
        #width_search = [20, 40, 80, 120, 200, 400, 800, 1600, 3200]
        width_search = bin
        wfdata_pol = np.zeros((dur*nSubs, LFFT, len(data_products)))
        #tdata = np.zeros((nSubs,2) )
         
        print("Collecting data")
        t0 = time.time()
        ## Read in the data and apply what ever scaling is needed
        for i in range(skip, skip+dur):
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
                t = subint[1] + tInt*(j-nSubs//2)
                d = data[:,:,j]*bscl + bzero
                
                #print(d.shape) 
                #if c == 0:
                #    tdata[j,:] = [tStartI, tStartF + t]
                for l,p in enumerate(data_products):
                    wfdata_pol[k,:,l] = d[l,:]


        #Assuming data_products == 'I':
        if  len(data_products) == 1:
            wfdata = wfdata_pol[:,:,0] 
        else:
            sys.exit("Check the polarization data products")

        #Computing and applying the bandpass
        bp = comp_bp(wfdata, nSubs)
        wfdata /= bp
        
        #Flagging bad rfi and time
        wfdata = flag_rfi_time(wfdata, nSubs)
        

        #wfdata = incoherent(freq, wfdata, tInt, dm, boundary='wrap', fill_value=np.nan) 

        tseries = np.mean(wfdata,axis = 1)
            

        #Looking for pulses in the box car

        print("starting boxcar searches")
        filt_dat =  boxcar_search(tseries, width_search)
        t1 = time.time()
        print(f"{dur*tSubs} s of psrfits data collected in {t1-t0}s")

        hdulist.close()

        wfavg = np.zeros((int(wfdata.shape[0]/10), int(wfdata.shape[1]/10)))
        xlist = np.arange(0,int(wfdata.shape[0]/10)*10,10)
        ylist = np.arange(0,int(wfdata.shape[1]/10)*10,10)
        for i,x in enumerate(xlist):
            for j,y in enumerate(ylist):
                wfavg[i,j] = np.mean(wfdata[x:x+10, y:y+10])
        
        tdat = np.arange(0,nSubs)*tInt*1e3

        freq_avg = np.zeros((int(freq.shape[0]/10)*10))
        flist = np.arange(0,len(freq_avg),10)
        for f,fq in enumerate(flist):
            freq_avg[f] = np.mean(freq[fq:fq+10])

        tavg = np.zeros((int(tdat.shape[0]/10)*10))
        tlist = np.arange(0,len(tavg), 10)

        for t,ts in enumerate(tlist):
            tavg[t] = np.mean(tdat[ts:ts+10])
        
        vmin = np.percentile(wfavg,0.1)
        vmax = np.percentile(wfavg,100)
        print (tavg)
        print (freq[0],freq[-1]) 
        plt.pcolormesh(wfavg.T, vmin = vmin, vmax = vmax)
        plt.xlabel("Time steps")
        plt.ylabel("Frequency channels")
        #plt.xticks(tavg)
        #plt.yticks(freq_avg)
        plt.show() 
        
        plt.plot(range(filt_dat.shape[1]),filt_dat[0,:])
        plt.show()

if __name__ == "__main__":  
    main()


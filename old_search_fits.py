"""
Based on writeHDF5fromPSRfits by J.Dowell

"""

import sys
import os
import re
import time
import numpy as np
from astropy.time import Time as AstroTime
from astropy.io import fits as astrofits
from matplotlib import pyplot as plt
from utils import comp_bp, flag_rfi_time, boxcar_search, get_cand_list, plot_events, matchFilter_search  
from dedispersion import incoherent



_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')
def main(args):
    # Parse command line options
    filenames = args


    # Open the file and load in basic information about the observation's goal
    for c,filename in enumerate(filenames):
        print(filename)

        ## Ready the PSRFITS file
        hdulist = astrofits.open(filename, memmap=True)

        #hdulist.info()
        #print(hdulist[0].header)
        #print (hdulist[1].header)
        #subint = hdulist[1].data[0]
        #data = subint[16]
        #print(data.shape)
        
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
        

        skip = 12 #0
        dur =  10 #nChunks
        #width_search = [20, 40, 80, 120, 200, 400, 800, 1600, 3200]
        width_search = [40]
        wfdata = np.zeros((dur*nSubs, LFFT, len(data_products)), dtype = np.float32)
        tdat = np.arange(0,dur*nSubs)*tInt


        print("Collecting data")
        t0 = time.time()
        ## Read in the data and apply what ever scaling is needed
        for i in range(skip, skip+dur):
            print(i)
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
            
        bp = comp_bp(wfdata,nSubs)
        
        
        #Correcting the bandpass
        wfdata /= bp

        #plt.pcolormesh(wfdata)
        #plt.show()

        #Flagging bad rfi and time chunks
        wfdata = flag_rfi_time(wfdata,nSubs)
        
        #plt.pcolormesh(wfdata)
        #plt.show()


        #Apply incoherent dedispersion
        #print (freq)
        #print (tInt, dm)

        wfdata = incoherent(freq, wfdata, tInt, dm, boundary='wrap', fill_value=np.nan)

        #plt.pcolormesh(wfdata)
        #plt.show()

        tseries = np.mean(wfdata,axis = 1)
        
        #plt.plot(range(tseries.shape[0]), tseries)
        #plt.show()

        #del wfdata


        #Looking for pulses in the box car

        print("starting boxcar searches")
        #filt_dat =  boxcar_search(tseries, width_search)
        filt_dat = matchFilter_search(tseries, width_search, tInt*1e+3)

        #plt.plot(range(filt_dat[1,:].shape[0]), filt_dat[1,:])
        plt.plot(range(filt_dat[0,:].shape[0]), filt_dat[0,:])
        plt.ylabel("Power (a.u)")
        plt.xlabel("Time integrations")
        plt.show()
        
        t2 = time.time()
        print(f"Bandpass, flagging, dedispersion and boxcar takes {t2-t1} s")

        
        dirname = os.path.dirname(filename)
        outdir = dirname+'/Images/'
        try:
            os.mkdir(outdir)
        except OSError:
            print ("Directory already exists")
       

        file_cand = open(dirname+"/candidates.txt","w")

        #Calling the plotting class
        pob = plot_events(wfdata, tdat, freq, os.path.basename(filename), outdir)


        for b,bin in enumerate(width_search):
            std = np.std(filt_dat[b,:])
            mean = np.mean(filt_dat[b,:])
            #print(filt_dat[b,:].max(), mean, std)
            snr = (filt_dat[b,:] - mean)/std
            good_sig = np.where((snr > 4.0))[0]
                
            if  len(good_sig) > 0:
                best_cand, best_snr  = get_cand_list(good_sig,snr[good_sig])
                print (f" Number of good sig: {len(good_sig)}, Number of best_sig: {len(best_cand)}, Bin size :{bin}" )

                for ind, cand in enumerate(best_cand):
                    # Getting each good candidate and making a plot for each of them
                    pob.plot(filt_dat[b,:], best_snr[ind], bin, max(0, cand-5000), min(cand+5000, wfdata.shape[0]))
                    file_cand.write(f" {cand}  {best_snr[ind]} {bin} \n")
        
        
        t3 = time.time()
        print (f"Collection of candidates and plotting them takes {t3-t2} s")
        del wfdata

        hdulist.close()
        file_cand.close()
        print(f"Finished processing {filename}")


if __name__ == "__main__":  
    args = sys.argv[1:]
    main(args)


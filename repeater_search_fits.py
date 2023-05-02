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
from utils import comp_bp, flag_rfi_time, boxcar_search, get_cand_list, plot_events, matchFilter_search, calc_triggerlag, rebin_data, smooth_time  
from dedispersion import incoherent, delay
import matplotlib
from detailed_info import load_repeaters
matplotlib.use('Agg')


_fnRE = re.compile('.*_b(?P<beam>[1-4])(t(?P<tuning>[12]))?_.*\.fits')

def main(args):
    # Parse command line options
    filenames = args.filenames

    #opening a directory to save image candidates

    dirname = os.path.dirname(filenames[0])
    outdir1 = os.path.join(dirname, 'Images_box')
    outdir2 = os.path.join(dirname, 'Images_match')
    #outdir3 = os.path.join(dirname, 'data_slices')
    try:
        os.mkdir(outdir1)
    except OSError:
        print (f"{outdir1} directory already exists")

    try:
        os.mkdir(outdir2)
    except OSError:
        print (f"{outdir2} directory already exists")
    
    #try:
    #    os.mkdir(outdir3)
    #except OSError:
    #    print (f"{outdir3} directory already exists")
        
    #text file to save candidates
    file_cand1 = open(outdir1+"/box_candidates.txt","w")
    file_cand2 = open(outdir2+"/match_candidates.txt","w")


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
        
        #dm = hdulist[0].header['CHAN_DM']
        # This is usually zero since the intra channel dedispersion is not corrected for the spectrometer data
        # and we need to get the dm from the repeater list

        # The DM is stored in the target list
        
        repeaters = load_repeaters()
        dm = repeaters[sourceName]



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
        
        """
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
        """


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
        
        # Only valid for metafile associated with chime triggers
        """
        if  args.metafile:
            tdelay = delay([np.max(freq), 400e6], dm)[0]
            trigger_lag  = calc_triggerlag(args.metafile)

            exp_time = tdelay - trigger_lag
            exp_samp = int(exp_time/tInt)

            print(f"Expected time of pulse arrival  is {exp_time}")
        """

        meta_list = [sourceName, round(dm,1)]
        skip = 0 # Number of sub integrations to skip in starting 
        dur =  nChunks # Number of sub integrations to load the data
        

        # Different bin sizes for boxcar filtering
        width_search = [20, 40, 80, 120, 200, 400, 800, 1600, 3200, 6400]
        
        #Initialize the big array to capture data from different sub integrations
        wfdata = np.zeros((dur*nSubs, LFFT, len(data_products)), dtype = np.float32)
        tdat = np.arange(0,dur*nSubs)*tInt


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


        #Flagging bad rfi and time chunks
        wfdata = flag_rfi_time(wfdata,nSubs)

        #Apply incoherent dedispersion
        wfdata_cor, tdelay = incoherent(freq, wfdata, tInt, dm, boundary='wrap', fill_value=np.nan)
        
        t2 = time.time()
        print(f"Bandpass, flagging, dedispersion  {t2-t1} s")

        #print(wfdata.shape)
        
        dmax = tdelay.max()
        print(dmax*tInt)

        #wfdata_cor = wfdata_cor[:-dmax, :]
        #wfdata = wfdata[:-dmax, :]

        #print(dmax, wfdata.shape)

        #Plotting the data

        wfdata_avg, freq_avg, tdat_avg = rebin_data(wfdata, freq, tdat,  avgt = 200, avgf = 10)
        wfdata_cor_avg, freq_avg, tdat_avg = rebin_data(wfdata_cor, freq, tdat, avgt = 200,  avgf = 10)
        
        freq_mhz = freq_avg/1e+6
        delf = freq_mhz[1] - freq_mhz[0]
        delt = tdat_avg[1] - tdat_avg[0]

        f_ar = np.linspace(freq_mhz[0] - delf/2.0, freq_mhz[-1] + delf/2.0, len(freq_mhz) + 1)
        t_ar = np.linspace(tdat_avg[0] - delt/2.0, tdat_avg[-1] + delt/2.0, len(tdat_avg) + 1)

        plt.pcolormesh(t_ar, f_ar, wfdata_avg.T)
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (MHz)")
        plt.savefig(f"{os.path.basename(filename)}_cmap_uncor.png")
        plt.close()

        plt.pcolormesh(t_ar, f_ar, wfdata_cor_avg.T)
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (MHz)")
        plt.savefig(f"{os.path.basename(filename)}_cmap_cor.png")
        plt.close()

        #Averging the data across channels for box filtering
        
        tseries_uncor = np.nanmean(wfdata,axis = 1)
        plt.plot(tdat,tseries_uncor)
        plt.ylabel("Power (a.u.)")
        plt.xlabel("Time (s)")
        plt.savefig(f"{os.path.basename(filename)}_tseries_uncor.png")
        plt.close()


        tseries_cor = np.nanmean(wfdata_cor,axis = 1)
        
        plt.plot(tdat, tseries_cor)
        plt.ylabel("Power (a.u.)")
        plt.xlabel("Time (s)")
        plt.savefig(f"{os.path.basename(filename)}_tseries_cor.png")
        plt.close()
        
        #Smoothened time series
        tseries_smth = smooth_time(tseries_cor, winSize = 2500)


        tseries_cor -= tseries_smth

        plt.plot(tdat, tseries_cor)
        plt.ylabel("Power (a.u.)")
        plt.xlabel("Time (s)")
        plt.savefig(f"{os.path.basename(filename)}_tseries_cor_smth.png")
        plt.close()

        
        
        
        #Looking for pulses in the box car
        # From box car
        filt_dat1 =  boxcar_search(tseries_cor, width_search)
        
        # From match filtering function
        filt_dat2 = matchFilter_search(tseries_cor, width_search, tInt*1e+3)

        t3 = time.time()
        print(f"Averaging and Smoothening takes {t3-t2} s")
        

        #Calling the plotting class
        #pob1 = plot_events(wfdata, tdat, freq, os.path.basename(filename), outdir1, outdir3)
        
        pob1 = plot_events(wfdata, wfdata_cor, tdat, freq, os.path.basename(filename), outdir1, meta_list, mode = 'box')
        
        print('Starting Box candidate search')
        #Candidate filtering and plotting 
        for b,bin in enumerate(width_search):
            
            std = np.std(filt_dat1[b,:])
            mean = np.mean(filt_dat1[b,:])
            #print(filt_dat1[b,:].max(), mean, std)
            snr = (filt_dat1[b,:] - mean)/std
            
            #if args.metafile and exp_time < (nChunks*tSubs):

            # Plotting the location where pulse is expected to show up irrespective of the sigma values
            #    pob1.plot(filt_dat1[b,:], snr[exp_samp], bin, max(0, exp_samp-5000), min(exp_samp+6000, wfdata.shape[0]), True)

            good_sig = np.where((snr > 4.0))[0]
                
            if  len(good_sig) > 0:
                best_cand, best_snr  = get_cand_list(good_sig,snr[good_sig])
                print (f" Number of good sig: {len(good_sig)}, Number of best_sig: {len(best_cand)}, Bin size :{bin}" )

                for ind, cand in enumerate(best_cand):
                    # Getting each good candidate and making a plot for each of them
                    pob1.plot(filt_dat1[b,:], best_snr[ind], bin, max(0, cand-5000), min(cand+6000, wfdata.shape[0]), False)
                    file_cand1.write(f" {os.path.basename(filename)} {cand}  {best_snr[ind]} {bin} {tdat[cand]} \n")
        
        
        del pob1
        

        #Calling the plotting class
        pob2 = plot_events(wfdata, wfdata_cor, tdat, freq, os.path.basename(filename), outdir2, meta_list, mode = 'filter')

        print('Starting Match filtering search')

        #Candidate filtering and plotting 
        for b,bin in enumerate(width_search):
            std = np.std(filt_dat2[b,:])
            mean = np.mean(filt_dat2[b,:])
            #print(filt_dat2[b,:].max(), mean, std)
            snr = (filt_dat2[b,:] - mean)/std
            
            #if args.metafile and exp_time < (nChunks*tSubs):

            #    # Plotting the location where pulse is expected to show up irrespective of the sigma values
            #    pob2.plot(filt_dat2[b,:], snr[exp_samp], bin, max(0, exp_samp-5000), min(exp_samp+6000, wfdata.shape[0]), True)


            good_sig = np.where((snr > 4.0))[0]
                
            if  len(good_sig) > 0:
                best_cand, best_snr  = get_cand_list(good_sig,snr[good_sig])
                print (f" Number of good sig: {len(good_sig)}, Number of best_sig: {len(best_cand)}, Bin size :{bin}" )

                for ind, cand in enumerate(best_cand):
                    # Getting each good candidate and making a plot for each of them
                    pob2.plot(filt_dat2[b,:], best_snr[ind], bin, max(0, cand-5000), min(cand+6000, wfdata.shape[0]), False)
                    file_cand2.write(f" {os.path.basename(filename)} {cand}  {best_snr[ind]} {bin} {tdat[cand]} \n")

        
        
        del pob2

        t4 = time.time()
        print (f"Collection of candidates and plotting them takes {t4-t3} s")
        del wfdata

        
        


        hdulist.close()
        print(f"Finished processing {filename}")
    
    file_cand1.close()
    file_cand2.close()


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


import numpy as np
from mad import median_absolute_deviation as mad
from matplotlib import pyplot as plt

#For reading the tarball metafile
import ephem
from datetime import datetime
from astropy.coordinates import Angle as AstroAngle
from lsl.common.stations import lwa1, lwasv
from lsl.common.mcs import mjdmpm_to_datetime
from lsl.common import metabundle, metabundleADP
from lsl.reader.drx import FILTER_CODES


def comp_bp(data,nsubs):

    """
    Function to compute a smooth bandpass model
    """
    med_spectra =  np.median(data[:nsubs,:], axis = 0)
    smth = med_spectra*0.0
    winSize = 10 
    
    #Compute the smoothed bandpass model
    for i in range(smth.size):
        mn = max([0, i-winSize/2])
        mx = min([i+(winSize/2)+1, smth.size])
        smth[i] = np.median(med_spectra[int(mn):int(mx)])      
           
    return smth





def flag_rfi_time(data, nsubs):
    """
    Flagging bad freq channels and time chunks
    """
    #dur = data.shape[0]/nsubs #calculating number of subintegration chunks
    chunks = np.arange(0,data.shape[0],nsubs)
    for i in chunks:
        avg_spec = np.mean(data[i:i+nsubs,:], axis = 0)
        medf = np.median(avg_spec)
        sigmdf = mad(avg_spec)
        flgf = np.where(((avg_spec-medf)/sigmdf > 3.5))[0]
        #print(flgf)        
        data[i:i+nsubs,flgf] = medf

    for j in chunks:
        avg_time = np.mean(data[j:j+nsubs,:], axis = 1)
        medt = np.median(avg_time)
        sigmdt = mad(avg_time)
        flgt = np.where(((avg_time-medt)/sigmdt > 3.5))[0]
        flgt += j
        #print(flgt)
        data[flgt,:] = medt


    return data 


def boxcar_search(tseries, width):

    """
    Returns a boxcar average series for different bin widths

    """
    datp = len(tseries)
    boxdat = np.zeros((len(width),datp))
    for k,bin in enumerate(width):
        #winsize = int(bin/2.0)
        for i in range(len(tseries)):
            #mn = max([0, i-winsize])
            mx = min([i+bin, datp])
            boxdat[k,i] = np.mean(tseries[i:mx])
    
    return boxdat

"""
def boxcar_wf(wf, bin):

    
    Apply a boxcar filter on a waterfall plot


    datp = wf.shape[0] #length of the time series
    boxdat_wf = wf*0.0  # A new waterfall to store afte filtering
    #Collecting the averaged waterfall into a new waterfall array
    for i in range(datp):
        mx = min([i+bin, datp])
        boxdat_wf[i,:] = np.mean(wf[i:mx,:], axis = 0)

    return boxdat_wf
"""



def get_cand_list(tm,snr):
    """
    Making a list of candidates which are not 
    repeatitions in the time series
    """
    
    best_tm = []
    best_snr = []
    sample_now = tm[0]
    snr_now = 0
    sample_last = tm[0]
    
    for s,sample in enumerate(tm):
        if abs(sample-sample_last) < 100: # Everything close to 100 time samples is considered as the same pulse
            if snr[s] > snr_now:
               sample_now = sample 
               snr_now = snr[s]
            if  (s == len(tm) -1) :
                best_tm.append(sample_now)
                best_snr.append(snr_now)
 
        else:
             
            best_tm.append(sample_now)
            best_snr.append(snr_now)
            snr_now = snr[s]
            sample_now = sample
            if  (s == len(tm) -1):
                best_tm.append(sample)
                best_snr.append(snr[s])


        sample_last = sample

     


    return best_tm, best_snr


class plot_events:

    """
    Class to plot the candidate files after the FRB search
    """
    def __init__(self, wfdata, tdat, freq, filename, outdir):
        
        self.wfdata = wfdata
        self.freq = freq/1e+6
        self.tdat = tdat
        self.filename = filename
        self.outdir = outdir
        #self.datadir = datadir

    def rebin_data(self, wfcut, tcut,  avgf = 10):

        avgf = int(avgf)
        #wfcut = self.wfdata[tstart:tstop,:] # Getting a section of data needed for plotting
        #tcut = self.tdat[tstart:tstop]

        # Averaging the waterfall data
        wfavg = np.zeros((int(wfcut.shape[0]/avgf), int(wfcut.shape[1]/avgf)))
        for i in range(wfavg.shape[0]):
            for j in range(wfavg.shape[1]):
                wfavg[i,j] = np.mean(wfcut[avgf*i:avgf*(i+1), avgf*j:avgf*(j+1)])

        #Averaging the frequency
        freq_avg = np.zeros(int(self.freq.shape[0]/avgf))
        for f in range(freq_avg.shape[0]):
            freq_avg[f] = np.mean(self.freq[avgf*f:avgf*(f+1)])
        
        #Averaging the time
        tavg = np.zeros(int(tcut.shape[0]/avgf))
        for t in range(tavg.shape[0]):
            tavg[t] = np.mean(tcut[avgf*t:avgf*(t+1)])

        return wfavg, freq_avg, tavg



    def plot(self, filt_time, snr, bin_width, tstart, tstop, exp_cand = False):
        
        wfcut = self.wfdata[tstart:tstop,:] # Getting a section of data needed for plotting
        tcut = self.tdat[tstart:tstop]

        #if  filter == 'box':                    
        #    wf_smth = boxcar_wf(wfcut, bin_width)
        #else:
        #    wf_smth = matchFilter_wf(wfcut, bin_width, tInt) 

        
        #outname_ar = self.datadir+'/%s_%0.1f_%0.1f_%0.1f' % (self.filename, tstart+5000, snr, bin_width)  
        #np.savez(outname_ar, wf =  wfcut.astype('float32'), time = tcut.astype('float32'), freq = self.freq.astype('float32'))

        wfplt, freqplt, tplt = self.rebin_data(wfcut, tcut, avgf = 20) #change avgf if needed to try another one
        
        filt_cut = filt_time[tstart:tstop]
        #tfilt_cut  = self.tdat[tstart:tstop]
        tInt = tcut[1] - tcut[0]

        # plotting 2 panels of wfdata and filtered time series
        plt.subplot(2,1,1)
        plt.pcolormesh(wfplt.T, cmap = 'viridis')
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (MHz)")
        if not exp_cand:
           plt.title(f"SNR: {round(snr,1)}, Bin width: {round(bin_width*tInt*1e+3)} ms")
        else:
           plt.title(f"Expected cand SNR: {round(snr,1)}, Bin width: {round(bin_width*tInt*1e+3)} ms")


        xr = np.arange(0, tplt.shape[0], 100)
        yr = np.arange(0, freqplt.shape[0], 100)
        tplt = np.round(tplt,1)
        freqplt = np.round(freqplt,1)
        
        plt.xticks(xr, tplt[xr])
        plt.yticks(yr, freqplt[yr])

        plt.subplot(2,1,2)
        plt.plot(tcut, filt_cut)
        plt.xlabel("Time (s)")
        plt.ylabel("Power (a.u.)")
        #plt.xticks(xr, tplt[xr])

        plt.tight_layout()
        # saving figure
        outname_img = self.outdir+'/%s_%0.1f_%0.1f_%0.1f.png' % (self.filename, tstart+5000, snr, bin_width)
        plt.savefig(outname_img, dpi = 150)
        plt.close()




def scatter_model(t, taud = 100.0, beta = 0.42, phase = 0, A = 1.0):
    """
    Function for matched filtering: t^beta exp(-t/taud)
    """
    beta = abs(beta)
    taud = abs(taud)
    A = abs(A)
    return np.piecewise(t, [t < phase, t >= phase], [0, lambda t: A*(t-phase)**beta*np.exp(-(t-phase)/taud)])




def matchFilter_search(tseries, width, tInt):

    """
    Returns a boxcar average series for different bin widths

    """
    datp = len(tseries)
    boxdat = np.zeros((len(width),datp))

    for k,bin in enumerate(width):
        tm = np.linspace(0, bin*tInt, bin)
        filt_func = scatter_model(tm, taud = (bin*tInt)/4.0, beta = 0.8)  
    
        for i in range(datp):
            
            # Matched filtering encountering an issue in the last set of bins which has less number of data 
            # After the datp -bin  in filtered data is filled with the last value had data in all bins
            
            tbox = tseries[i:min(i+bin,datp)]
            if i <=  datp-bin:
                boxdat[k,i] = np.sum(tbox*filt_func[:len(tbox)])
            else:
                boxdat[k,i] = boxdat[k,datp-bin]

    return boxdat

"""
def matchFilter_wf(wf, bin, tInt):

    
    #Applying a match filter on a waterfall plot

    
    wf = wf.T
    datp = wf.shape[1]
    match_wf = wf*0.0

    tm = np.linspace(0, bin*tInt, bin)
    filt_func = scatter_model(tm, taud = (bin*tInt)/4.0, beta = 0.8)

    for i in range(datp):

        # Matched filtering encountering an issue in the last set of bins which has less number of data 
        # After the datp -bin  in filtered, data in rest of the bins if filled with the last value (from datp-bin
        # Otherwise it would cause an increase in summed value towards the end

        tbox = wf[:,i:min(i+bin,datp)] 
        if i <=  datp-bin:
           match_wf[:,i] = np.sum((tbox*filt_func), axis = 1)
        else:
            match_wf[:,i] = wf[:,datp-bin]

    return match_wf.T
"""


def calc_triggerlag(filename):

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

    # Parse and process
    try:
        ## Trigger
        event, dm, _ = name.split('_', 2)
        event = event.split('#', 1)[1]
        event = datetime.strptime(event, "%Y-%m-%d-%H:%M:%S.%fUTC")
        dm = float(dm.replace('DM', ''))
        is_trigger = True
    except ValueError:
        
        print("Metafile not associated with trigger. Check the repeaters list")
        sys.exit(0)

    start = mjdmpm_to_datetime(mjd, mpm)
    ra_str = AstroAngle(ra, unit='hourangle')
    ra_str = str(ra_str)
    dec_str = AstroAngle(dec, unit='deg')
    dec_str = str(dec_str)


    # Calculate the lags
    if is_trigger:
       trigger_lag = (start-event).total_seconds()
       #disp_lagU = delay([freqU+bw/2.0, 400e6], dm)[0]
       #disp_lagL = delay([freqU-bw/2.0, 400e6], dm)[0]
       #disp_lagL = delay([freqL-bw/2.0, 400e6], dm)[0]
       #is_repeater = False

       return trigger_lag
    


import sys
import time
import numpy as np
from utils import  boxcar_plot, boxcar_wf 
from matplotlib import pyplot as plt

def main():
    # Parse command line options
    filename = sys.argv[1]
    bin_width = int(sys.argv[3])
    avgf = int(sys.argv[2])

    # Open the file and load in basic information about the observation's goal
    #loading the npz file and the all info
    dat_npz = np.load(filename)
    
    ## File specifics
    source = dat_npz['source']
    ra = dat_npz['ra']
    dec = dat_npz['dec']

    cfreq = dat_npz['cfreq']  #Hz
    srate = dat_npz['srate']  #Hz
    lfft = dat_npz['lfft']
    tint = dat_npz['tint']
    dm = dat_npz['dm']

    #Loading the data arrays
    wfdat = dat_npz['wfdat'] # Waterfall data
    tdat = dat_npz['tdat']   # Time info
    freq = dat_npz['freq']   # Frequency info


    ## Report
    print("Filename: %s" % (filename))
    print("Sample Rate: %i Hz" % srate)
    print("Tuning Frequency: %.3f Hz" % cfreq)
    print("---")
    print("Target: %s" % source)
    print("RA: %.3f hours" % (ra/15.0,))
    print("Dec: %.3f degrees" % dec)
    print("DM: %.3f pc/cm^3" % dm)
    print("Integration Time: %.3f ms" % (tint*1e3,))
    print("Transform Length: %i" % lfft)
        


    width_search = [bin_width,]
    tseries = np.mean(wfdat,axis = 1)
    
    #Box car averaging on both time series and waterfall
    filt_dat =  boxcar_plot(tseries, width_search)
    filt_dat = filt_dat[0,:]

    wfdat = (wfdat - np.mean(wfdat))/np.std(wfdat)
    
    wfdat_filt =  boxcar_wf(wfdat, bin_width)

    # Averaging the waterfall data
    wfavg = np.zeros((int(wfdat.shape[0]/avgf), int(wfdat.shape[1]/avgf)))
    for i in range(wfavg.shape[0]):
        for j in range(wfavg.shape[1]):
             wfavg[i,j] = np.mean(wfdat[avgf*i:avgf*(i+1), avgf*j:avgf*(j+1)])


    wfdat_fft = np.fft.fft2(wfavg)
    #wfdat_fft = np.fft.fftshift(wfdat_fft)
    amp = np.abs(wfdat_fft)**1.3
    phase = np.angle(wfdat_fft)
    wfdat_fft_win = amp*(np.exp(1.0j*phase))
    wfdat_ifft = np.fft.ifft2(wfdat_fft_win)

    
    wfavg_filt = np.zeros((int(wfdat_filt.shape[0]/avgf), int(wfdat_filt.shape[1]/avgf)))
    for i in range(wfavg_filt.shape[0]):
        for j in range(wfavg_filt.shape[1]):
             wfavg_filt[i,j] = np.mean(wfdat_filt[avgf*i:avgf*(i+1), avgf*j:avgf*(j+1)])
    

    #Averaging the frequency
    freq_avg = np.zeros(int(freq.shape[0]/avgf))
    for f in range(freq_avg.shape[0]):
        freq_avg[f] = np.mean(freq[avgf*f:avgf*(f+1)])

    #Averaging the time
    tavg = np.zeros(int(tdat.shape[0]/avgf))
    for t in range(tavg.shape[0]):
        tavg[t] = np.mean(tdat[avgf*t:avgf*(t+1)])
        
    vmin = np.percentile(wfavg, 0)
    vmax = np.percentile(wfavg, 100)

    # plotting 2 panels of wfdata and filtered time series
    fig, axs = plt.subplots(2, 3, constrained_layout=True, figsize = (12,10))
    
    
    tavg = np.round(tavg,1)
    freq_avg = np.round(freq_avg/1e+6,1)
    
    xr = np.arange(0, tavg.shape[0], 100)
    yr = np.arange(0, freq_avg.shape[0], 100)

    #tmesh, freq_mesh  = np.meshgrid(tavg, freq_avg)
    #print(freq_mesh.shape, tmesh.shape, wfavg.shape)   

    
    #axs[0,0].pcolormesh(10*np.log10(wfavg.T), cmap = 'viridis')
    axs[0,0].pcolormesh(np.abs(wfavg).T, cmap = 'viridis', shading = 'gouraud')
    axs[0,0].set_xlabel("Time (s)")
    axs[0,0].set_ylabel("Frequency (MHz)")
    axs[0,0].set_title(f"Without smoothing")    


    
    axs[0,0].set_xticks(xr)
    axs[0,0].set_xticklabels(tavg[xr])
    axs[0,0].set_yticks(yr)
    axs[0,0].set_yticklabels(freq_avg[yr])
    
    
    axs[0,1].pcolormesh(np.abs(wfavg_filt.T), cmap = 'viridis')
    axs[0,1].set_xlabel("Time (s)")
    axs[0,1].set_ylabel("Frequency (MHz)")
    axs[0,1].set_title(f"With Smoothing")


    axs[0,2].pcolormesh(np.abs(wfdat_ifft).T, cmap = 'viridis')
    axs[0,2].set_xlabel("Time (s)")
    axs[0,2].set_ylabel("Frequency (MHz)")
    axs[0,2].set_title(f"With Fourier enhancing")


    #xr = np.arange(0, tavg.shape[0], 100)
    #yr = np.arange(0, freq_avg.shape[0], 100)

    #tavg = np.round(tavg,1)
    #freq_avg = np.round(freq_avg/1e+6,1)

    #ax0.set_yticks(xr, tavg)
    #ax0.set_xticks(yr, freq_avg)

    axs[0,1].set_xticks(xr)
    axs[0,1].set_xticklabels(tavg[xr])
    axs[0,1].set_yticks(yr)
    axs[0,1].set_yticklabels(freq_avg[yr])




    axs[1,0].plot(tdat, tseries)
    axs[1,0].set_xlabel("Time (s)")
    axs[1,0].set_ylabel("Power (a.u.)")
    
    axs[1,1].plot(tdat, filt_dat)
    axs[1,1].set_xlabel("Time (s)")
    axs[1,1].set_ylabel("Smoothed Power (a.u.)")


    axs[1,2].plot(tavg, np.mean(np.abs(wfdat_ifft), axis = 1), '.')
    axs[1,2].set_xlabel("Time (s)")
    axs[1,2].set_ylabel("Enhanced Power (a.u.)")
    
    fig.suptitle(f"FRB: {source}")
    # saving figure
    #outname_img = self.outdir+'/%s_%0.1f_%0.1f_%0.1f.png' % (self.filename, tstart+5000, snr, bin_width)
    #plt.savefig(outname_img, dpi = 150)
    #plt.close()
    plt.show()

if __name__ == "__main__":  
    main()


import sys
import time
import numpy as np
from utils import  boxcar_plot, boxcar_wf 
from matplotlib import pyplot as plt
from scipy import signal
from skimage import exposure

def main():
    # Parse command line options
    filename = sys.argv[1]
    #bin_width = int(sys.argv[2])
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
        


    #width_search = [bin_width,]
    #tseries = np.mean(wfdat,axis = 1)
    
    ##Box car averaging on both time series and waterfall
    #filt_dat =  boxcar_plot(tseries, width_search)
    #filt_dat = filt_dat[0,:]
    
    # Averaging the waterfall data
    wfavg = np.zeros((int(wfdat.shape[0]/avgf), int(wfdat.shape[1]/avgf)))
    for i in range(wfavg.shape[0]):
        for j in range(wfavg.shape[1]):
             wfavg[i,j] = np.mean(wfdat[avgf*i:avgf*(i+1), avgf*j:avgf*(j+1)])
    
    wfavg = (wfavg - wfavg.min())/(wfavg.max() - wfavg.min())

    #Contrast stretching
    p2, p98 = np.percentile(wfavg, (2, 98))
    wfavg_rescale = exposure.rescale_intensity(wfavg, in_range=(p2, p98))

    # Equalization
    wfavg_eq = exposure.equalize_hist(wfavg)

    # Adaptive Equalization
    wfavg_adapteq = exposure.equalize_adapthist(wfavg, clip_limit=0.03)


    plt.subplot(2,3,1)
    plt.pcolormesh(wfavg.T, cmap = 'viridis')
    plt.title("Averaged data")

    plt.subplot(2,3,2)
    plt.pcolormesh(wfavg_rescale.T, cmap = 'viridis')
    plt.title("Contrast stretching")

    plt.subplot(2,3,3)
    plt.pcolormesh(wfavg_eq.T, cmap = 'viridis')
    plt.title("Histogram Equal..")

    plt.subplot(2,3,4)
    plt.pcolormesh(wfavg_adapteq.T, cmap = 'viridis')
    plt.title("Adaptive Histogram...")

    #plt.show()


    xf = np.arange(wfavg.shape[0])
    yf = np.arange(wfavg.shape[1])
   
    wind1 = 1/signal.windows.tukey(wfavg.shape[0], 0.3)
    wind2 = 1/signal.windows.tukey(wfavg.shape[1], 0.3)

    wind2d = np.outer(wind1, wind2)
    print(wfavg.shape)
    print(wind2d.shape)

    x, y  = np.meshgrid(xf, yf)
    z = wind2d
    
    #fig = plt.figure()
    #ax = plt.axes(projection = '3d')
    #ax.plot_surface(x,y,z.T)
    #plt.show()



    wfdat_fft = np.fft.fft2(wfavg)
    wfdat_fft = np.fft.fftshift(wfdat_fft)
    #wfdat_fft_win = np.zeros(wfdat_fft.shape)
    amp = np.abs(wfdat_fft)**1.5
    phase = np.angle(wfdat_fft)
    wfdat_fft_win = amp*(np.exp(1.0j*phase))

    #wfdat_fft_win[240:260,:] = 0

    wfdat_ifft = np.fft.ifft2(wfdat_fft_win)
    #wfdat_ifft = np.fft.ifftshift(wfdat_ifft)

    plt.subplot(2,3,5)
    plt.pcolormesh(np.abs(wfdat_ifft).T, cmap = 'viridis')
    plt.title("Fourier domain enhancing")
    plt.show()


     
    plt.subplot(2,2,1)
    plt.pcolormesh(np.abs(wfavg).T, cmap = 'viridis')
    plt.subplot(2,2,2)
    plt.pcolormesh(np.log(np.abs(wfdat_fft).T), cmap = 'viridis')
    
    plt.subplot(2,2,3)
    plt.pcolormesh(np.log(np.abs(wfdat_fft_win).T), cmap = 'viridis')

    plt.subplot(2,2,4)
    plt.pcolormesh(np.abs(wfdat_ifft).T, cmap = 'viridis')
    plt.show()

    #wfdat = (wfdat - np.mean(wfdat))/np.std(wfdat)
    
    """    
    wfdat_filt =  boxcar_wf(wfdat, bin_width)

    # Averaging the waterfall data
    wfavg = np.zeros((int(wfdat.shape[0]/avgf), int(wfdat.shape[1]/avgf)))
    for i in range(wfavg.shape[0]):
        for j in range(wfavg.shape[1]):
             wfavg[i,j] = np.mean(wfdat[avgf*i:avgf*(i+1), avgf*j:avgf*(j+1)])

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
        
    vmin = np.percentile(wfavg, 2)
    vmax = np.percentile(wfavg, 98)

    # plotting 2 panels of wfdata and filtered time series
    fig, axs = plt.subplots(2, 2, constrained_layout=True, figsize = (12,10))
    
    #axs[0,0].pcolormesh(10*np.log10(wfavg.T), cmap = 'viridis')
    axs[0,0].pcolormesh(wfavg.T, cmap = 'gray')
    axs[0,0].set_xlabel("Time (s)")
    axs[0,0].set_ylabel("Frequency (MHz)")
    axs[0,0].set_title(f"Without smoothing")    


    xr = np.arange(0, tavg.shape[0], 100)
    yr = np.arange(0, freq_avg.shape[0], 100)
        
    tavg = np.round(tavg,1)
    freq_avg = np.round(freq_avg/1e+6,1)

    #ax0.set_yticks(xr, tavg)
    #ax0.set_xticks(yr, freq_avg)
    
    axs[0,0].set_xticks(xr)
    axs[0,0].set_xticklabels(tavg[xr])
    axs[0,0].set_yticks(yr)
    axs[0,0].set_yticklabels(freq_avg[yr])
    
    #axs[0,1].pcolormesh(10*np.log10(wfavg_filt.T), cmap = 'viridis')
    axs[0,1].pcolormesh(wfavg_filt.T, cmap = 'gray')
    axs[0,1].set_xlabel("Time (s)")
    axs[0,1].set_ylabel("Frequency (MHz)")
    axs[0,1].set_title(f"With smoothing")


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
    axs[1,1].set_ylabel("Power (a.u.)")
    #plt.xticks(xr, tplt[xr])
    
    fig.suptitle(f"FRB: {source}, Bin width: {round(bin_width*tint*1e+3)} ms")
    # saving figure
    #outname_img = self.outdir+'/%s_%0.1f_%0.1f_%0.1f.png' % (self.filename, tstart+5000, snr, bin_width)
    #plt.savefig(outname_img, dpi = 150)
    #plt.close()
    plt.show()
    """

if __name__ == "__main__":  
    main()


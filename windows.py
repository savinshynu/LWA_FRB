#from __future__ import division, print_function, absolute_import


#import warnings

import numpy as np
"""
#from scipy import fftpack, linalg, special
#from scipy._lib.six import string_types

__all__ = ['boxcar', 'triang', 'parzen', 'bohman', 'blackman', 'nuttall',
           'blackmanharris', 'flattop', 'bartlett', 'hanning', 'barthann',
           'hamming', 'kaiser', 'gaussian', 'general_gaussian', 'chebwin',
           'slepian', 'cosine', 'hann', 'exponential', 'tukey', 'get_window']

"""
def _len_guards(M):
    """Handle small or incorrect window lengths"""
    if int(M) != M or M < 0:
        raise ValueError('Window length M must be a non-negative integer')
    return M <= 1


def _extend(M, sym):
    """Extend window by 1 sample if needed for DFT-even symmetry"""
    if not sym:
        return M + 1, True
    else:
        return M, False


def _truncate(w, needed):
    """Truncate window by 1 sample if needed for DFT-even symmetry"""
    if needed:
        return w[:-1]
    else:
        return w


def _cos_win(M, a, sym=True):
    r"""
    Generic weighted sum of cosine terms window
    Parameters
    ----------
    M : int
        Number of points in the output window
    a : array_like
        Sequence of weighting coefficients. This uses the convention of being
        centered on the origin, so these will typically all be positive
        numbers, not alternating sign.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.
    References
    ----------
    .. [1] A. Nuttall, "Some windows with very good sidelobe behavior," IEEE
           Transactions on Acoustics, Speech, and Signal Processing, vol. 29,
           no. 1, pp. 84-91, Feb 1981. :doi:`10.1109/TASSP.1981.1163506`.
    .. [2] Heinzel G. et al., "Spectrum and spectral density estimation by the
           Discrete Fourier transform (DFT), including a comprehensive list of
           window functions and some new flat-top windows", February 15, 2002
           https://holometer.fnal.gov/GH_FFT.pdf
    Examples
    --------
    Heinzel describes a flat-top window named "HFT90D" with formula: [2]_
    .. math::  w_j = 1 - 1.942604 \cos(z) + 1.340318 \cos(2z)
               - 0.440811 \cos(3z) + 0.043097 \cos(4z)
    where
    .. math::  z = \frac{2 \pi j}{N}, j = 0...N - 1
    Since this uses the convention of starting at the origin, to reproduce the
    window, we need to convert every other coefficient to a positive number:
    >>> HFT90D = [1, 1.942604, 1.340318, 0.440811, 0.043097]
    The paper states that the highest sidelobe is at -90.2 dB.  Reproduce
    Figure 42 by plotting the window and its frequency response, and confirm
    the sidelobe level in red:
    >>> from scipy import signal
    >>> from scipy.fftpack import fft, fftshift
    >>> import matplotlib.pyplot as plt
    >>> window = signal._cos_win(1000, HFT90D, sym=False)
    >>> plt.plot(window)
    >>> plt.title("HFT90D window")
    >>> plt.ylabel("Amplitude")
    >>> plt.xlabel("Sample")
    >>> plt.figure()
    >>> A = fft(window, 10000) / (len(window)/2.0)
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
    >>> plt.plot(freq, response)
    >>> plt.axis([-50/1000, 50/1000, -140, 0])
    >>> plt.title("Frequency response of the HFT90D window")
    >>> plt.ylabel("Normalized magnitude [dB]")
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    >>> plt.axhline(-90.2, color='red')
    """
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    fac = np.linspace(-np.pi, np.pi, M)
    w = np.zeros(M)
    for k in range(len(a)):
        w += a[k] * np.cos(k * fac)

    return _truncate(w, needs_trunc)

def hann(M, sym=True):
    r"""
    Return a Hann window.
    The Hann window is a taper formed by using a raised cosine or sine-squared
    with ends that touch zero.
    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an empty
        array is returned.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.
    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if `M` is even and `sym` is True).
    Notes
    -----
    The Hann window is defined as
    .. math::  w(n) = 0.5 - 0.5 \cos\left(\frac{2\pi{n}}{M-1}\right)
               \qquad 0 \leq n \leq M-1
    The window was named for Julius von Hann, an Austrian meteorologist. It is
    also known as the Cosine Bell. It is sometimes erroneously referred to as
    the "Hanning" window, from the use of "hann" as a verb in the original
    paper and confusion with the very similar Hamming window.
    Most references to the Hann window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.
    References
    ----------
    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
           spectra, Dover Publications, New York.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics",
           The University of Alberta Press, 1975, pp. 106-108.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 425.
    Examples
    --------
    Plot the window and its frequency response:
    >>> from scipy import signal
    >>> from scipy.fftpack import fft, fftshift
    >>> import matplotlib.pyplot as plt
    >>> window = signal.hann(51)
    >>> plt.plot(window)
    >>> plt.title("Hann window")
    >>> plt.ylabel("Amplitude")
    >>> plt.xlabel("Sample")
    >>> plt.figure()
    >>> A = fft(window, 2048) / (len(window)/2.0)
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
    >>> plt.plot(freq, response)
    >>> plt.axis([-0.5, 0.5, -120, 0])
    >>> plt.title("Frequency response of the Hann window")
    >>> plt.ylabel("Normalized magnitude [dB]")
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    """
    # Docstring adapted from NumPy's hanning function
    if _len_guards(M):
        return np.ones(M)
    M, needs_trunc = _extend(M, sym)

    w = _cos_win(M, [0.5, 0.5])

    return _truncate(w, needs_trunc)


hanning = hann


def tukey(M, alpha=0.5, sym=True):
    r"""Return a Tukey window, also known as a tapered cosine window.
    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an empty
        array is returned.
    alpha : float, optional
        Shape parameter of the Tukey window, representing the fraction of the
        window inside the cosine tapered region.
        If zero, the Tukey window is equivalent to a rectangular window.
        If one, the Tukey window is equivalent to a Hann window.
    sym : bool, optional
        When True (default), generates a symmetric window, for use in filter
        design.
        When False, generates a periodic window, for use in spectral analysis.
    Returns
    -------
    w : ndarray
        The window, with the maximum value normalized to 1 (though the value 1
        does not appear if `M` is even and `sym` is True).
    References
    ----------
    .. [1] Harris, Fredric J. (Jan 1978). "On the use of Windows for Harmonic
           Analysis with the Discrete Fourier Transform". Proceedings of the
           IEEE 66 (1): 51-83. :doi:`10.1109/PROC.1978.10837`
    .. [2] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function#Tukey_window
    Examples
    --------
    Plot the window and its frequency response:
    >>> from scipy import signal
    >>> from scipy.fftpack import fft, fftshift
    >>> import matplotlib.pyplot as plt
    >>> window = signal.tukey(51)
    >>> plt.plot(window)
    >>> plt.title("Tukey window")
    >>> plt.ylabel("Amplitude")
    >>> plt.xlabel("Sample")
    >>> plt.ylim([0, 1.1])
    >>> plt.figure()
    >>> A = fft(window, 2048) / (len(window)/2.0)
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
    >>> plt.plot(freq, response)
    >>> plt.axis([-0.5, 0.5, -120, 0])
    >>> plt.title("Frequency response of the Tukey window")
    >>> plt.ylabel("Normalized magnitude [dB]")
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    """
    if _len_guards(M):
        return np.ones(M)

    if alpha <= 0:
        return np.ones(M, 'd')
    elif alpha >= 1.0:
        return hann(M, sym=sym)

    M, needs_trunc = _extend(M, sym)

    n = np.arange(0, M)
    width = int(np.floor(alpha*(M-1)/2.0))
    n1 = n[0:width+1]
    n2 = n[width+1:M-width-1]
    n3 = n[M-width-1:]

    w1 = 0.5 * (1 + np.cos(np.pi * (-1 + 2.0*n1/alpha/(M-1))))
    w2 = np.ones(n2.shape)
    w3 = 0.5 * (1 + np.cos(np.pi * (-2.0/alpha + 1 + 2.0*n3/alpha/(M-1))))

    w = np.concatenate((w1, w2, w3))

    return _truncate(w, needs_trunc)


if __name__ == "__main__":
 import matplotlib.pyplot as plt

 x = tukey(51)
 plt.ylim([0, 1.1])
 plt.plot(x)
 plt.show()

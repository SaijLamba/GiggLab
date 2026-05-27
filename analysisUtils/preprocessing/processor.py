from scipy import signal, fft
import numpy as np
import warnings

class NeuropixelsProcessor:
    """
    processor class, processes chunks at a time to reduce computation time from doing it to the whole recording
    init defines self, n channels, parameters for recording, shifts for deskewing and the parametors for butterworth filter
    """
    def __init__(self, n_channels, params):
        self.n_ch = n_channels
        self.p = params
        #defines numerator (b) and denominator (a) for butterworth polynomial
        self.b, self.a = signal.butter(3, 300, btype='low', fs=params['fs_raw'])
        self.bad_mask = None
        #defines shifts for deskewing
        self.shifts = np.arange(n_channels) / n_channels

    def detect_bad_channels(self, raw_chunk):
        """
        detects bad channels of chunk
            1. scales to voltage from int16
            2. computes rms
            3. computes correlations against common median
            4. defines bad channels via threshold in params
            5. returns bad channels/mask as self.bad_mask, and counts of bad channels
        """
        #scale to volts
        data = raw_chunk.astype(np.float32) * self.p['voltage_scaling']
        #computer rms
        rms = np.sqrt(np.mean(data**2, axis=1))
        median_rms = np.nanmedian(rms)
        #correlate to common median
        global_sig = np.median(data, axis=0)
        corrs = np.zeros(self.n_ch)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i in range(self.n_ch):
                if rms[i] > 0.1: corrs[i] = np.corrcoef(data[i], global_sig)[0,1]
        #define threshold
        t = self.p['bad_ch_thresh']
        #assing bad channel mask
        self.bad_mask = (rms < t['rms_low']) | (rms > median_rms * t['rms_high_factor']) | (corrs < t['corr_min'])
        return {'bad_mask': self.bad_mask, 
                'counts': (np.sum(rms < t['rms_low']), np.sum(rms > median_rms * t['rms_high_factor']), np.sum(corrs < t['corr_min'])),
                'metrics': {
                    'rms': rms,
                    'corrs': corrs,
                    'median_rms': median_rms
                }
               }
 
    def process_chunk(self, raw_chunk, clean=True):
        """
        chunk processor
            1. scales to volatage
            2. clean = false then applies butterworth filter and downsamples, return causes function break i.e stops doing cleaning if already clean
            3. clean = True cleans data by:
                1. deskewing data
                2. setting bad channels to nan
                3. detrending and cmr
                4. applies butterworth filter
                5. returns downsampled data
        thus end return of chunk is same, resulting in a processed lfp for the chunk, but this method handles if the chunk has already been cleaned
        """
        data = raw_chunk.astype(np.float32) * self.p['voltage_scaling']
        if not clean:
            data = signal.filtfilt(self.b, self.a, data, axis=1)
            q = int(self.p['fs_raw'] / self.p['fs_lfp'])
            return data[:, ::q]

        # Deskew
        n_t = data.shape[1]
        f = fft.rfft(data, axis=1)
        freqs = fft.rfftfreq(n_t)
        shift = np.exp(-1j * 2 * np.pi * freqs[None, :] * self.shifts[:, None])
        data = fft.irfft(f * shift, n=n_t, axis=1)
        
        # Mask Bad
        if self.bad_mask is not None: data[self.bad_mask, :] = np.nan
            
        # Detrend & CMR
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = data - np.nanmedian(data, axis=1, keepdims=True)
            data = data - np.nanmedian(data, axis=0, keepdims=True)
            
        # Filter
        data = np.nan_to_num(data, nan=0.0) 
        data = signal.filtfilt(self.b, self.a, data, axis=1)
        if self.bad_mask is not None: data[self.bad_mask, :] = np.nan
            
        q = int(self.p['fs_raw'] / self.p['fs_lfp'])
        return data[:, ::q]



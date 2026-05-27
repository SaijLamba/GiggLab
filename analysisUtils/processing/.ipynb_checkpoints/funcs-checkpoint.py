from joblib import Parallel, delayed
import warnings
import numpy as np
from scipy import stats
from scipy.ndimage import gaussian_filter

def avg_lfp_by_trigger(data, triggers, params, processor):
    """
    averages out lfp by trigger using params
    designed to use memmap as input for data, in theory can be used with any binary/numpy data of shape [samples, channels]
    uses flattened triggers ie. [trigger 1, trigger 2] etc, can not yet handle trials, will implement in future (or write seperate func)
    in theory could define epochs outside of function and interate over trials to produce
    returns mean_cleaned lfp data over trigger window
    needs specific processor for this recording as part of input, in theory can be logiced out I recon, however as redundency (different recordings may use different processors)
    """
    
    def _process_single_trigger(trig):
            s = int(trig + params['window_lfp'][0] * params['fs_raw'])
            e = int(trig + params['window_lfp'][1] * params['fs_raw'])
            
            # Return None if out of bounds so we can filter it later
            if s < 0 or e > data.shape[0]: 
                return None
                
            raw_trace = data[s:e, :384].T
            return processor.process_chunk(raw_trace, clean=True)
    
    # 2. Execute the loop in parallel across all CPU cores
    epochs_clean_raw = Parallel(n_jobs=-1, prefer="threads")(
        delayed(_process_single_trigger)(trig) for trig in triggers
    )
    
    # 3. Filter out the skipped (out of bounds) triggers
    epochs_clean = [epoch for epoch in epochs_clean_raw if epoch is not None]
            
    # 4. Original mathematical summary remains identical
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean_clean = np.nanmean(epochs_clean, axis=0)
        
        sem_clean = stats.sem(epochs_clean, axis=0, nan_policy='omit')
        
        if np.ma.isMaskedArray(sem_clean):
            sem_clean = sem_clean.filled(np.nan)
            
    return mean_clean, sem_clean
    
def calculate_spiking_mua_depth(spike_times, spike_positions, triggers, fs=30000, 
                                pre_stim=100, post_stim=500, time_bin_ms=5, depth_bin_um=20):
    """
    Calculates trial-averaged, depth-resolved MUA from discrete spike times.
    
    Parameters:
    spike_times (np.array): Array of spike timestamps in samples.
    spike_positions (np.array): Nx2 array where col 1 is the y-coordinate (depth).
    triggers (np.array): Array of stimulus onset timestamps in samples.
    fs (int): Sampling rate in Hz.
    pre_stim (float): Milliseconds to extract before stimulus.
    post_stim (float): Milliseconds to extract after stimulus.
    time_bin_ms (float): Resolution of the time axis in ms.
    depth_bin_um (float): Resolution of the depth axis in um.
    
    Returns:
    mua_map (np.ndarray): 2D smoothed array (depth x time) of firing rates.
    time_ax (np.ndarray): 1D array of time bin centers.
    depth_ax (np.ndarray): 1D array of depth bin centers.
    """
    
    print("Calculating depth-resolved Spike Density MUA...")
    
    # Extract just the depths
    spike_depths = spike_positions[:, 1]
    
    # Convert window parameters from ms to samples
    pre_samples = int((pre_stim / 1000.0) * fs)
    post_samples = int((post_stim / 1000.0) * fs)
    
    # Define physical bins for the 2D histogram
    # We use min/max of the actual spike positions to cover the whole probe
    min_depth = np.floor(np.min(spike_depths) / depth_bin_um) * depth_bin_um
    max_depth = np.ceil(np.max(spike_depths) / depth_bin_um) * depth_bin_um
    
    depth_edges = np.arange(min_depth, max_depth + depth_bin_um, depth_bin_um)
    time_edges = np.arange(-pre_stim, post_stim + time_bin_ms, time_bin_ms)
    
    # Initialize an empty 2D histogram (Depth x Time)
    total_hist = np.zeros((len(depth_edges)-1, len(time_edges)-1))
    
    # Iterate through each checkerboard reversal trial
    for trigger in triggers:
        window_start = trigger - pre_samples
        window_end = trigger + post_samples
        
        # Find spikes that occurred within this specific trial window
        mask = (spike_times >= window_start) & (spike_times <= window_end)
        trial_spikes = spike_times[mask]
        trial_depths = spike_depths[mask]
        
        # Convert spike times to relative milliseconds from stimulus onset
        relative_times_ms = ((trial_spikes - trigger) / fs) * 1000.0
        
        # Bin the spikes for this trial in 2D space (time and depth)
        trial_hist, _, _ = np.histogram2d(
            trial_depths, relative_times_ms, 
            bins=[depth_edges, time_edges]
        )
        total_hist += trial_hist
        
    # Convert absolute spike counts to Average Firing Rate (Spikes/Second/Trial)
    num_trials = len(triggers)
    bin_width_seconds = time_bin_ms / 1000.0
    mua_rate = total_hist / (num_trials * bin_width_seconds)
    
    # Apply a 2D Gaussian filter to smooth the discrete bins into a continuous density map
    # Sigma is given in units of bins (e.g., sigma=1 smooths over adjacent bins)
    mua_smoothed = gaussian_filter(mua_rate, sigma=(1.5, 2.0))
    
    # Calculate bin centers for accurate plotting axes
    depth_ax = depth_edges[:-1] + (depth_bin_um / 2)
    time_ax = time_edges[:-1] + (time_bin_ms / 2)
    
    print("MUA Spiking map generated.")
    return mua_smoothed, time_ax, depth_ax



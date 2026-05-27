import numpy as np
import warnings
import pandas as pd
import os

def create_virtual_probe(data, channel_depths):
    """
    Creates a virtual 1D probe by averaging electrodes that share exact 
    vertical depth coordinates, natively ignoring NaN-masked bad channels.
    
    Parameters:
    data (np.ndarray): 2D array of shape (n_channels, n_samples) containing LFP data.
    channel_depths (np.ndarray): 1D array of shape (n_channels,) with integer depth values.
    
    Returns:
    virtual_data (np.ndarray): 1D collapsed probe data (n_unique_depths, n_samples).
    unique_depths (np.ndarray): The specific integer depth coordinates.
    """
    # 1. Enforce integer matching for precise spatial filtering
    depths_int = np.asarray(channel_depths, dtype=int)
    unique_depths = np.sort(np.unique(depths_int))
    n_depths = len(unique_depths)
    
    # 2. Initialize output
    virtual_data = np.zeros((n_depths, data.shape[1]))
    
    # 3. Collapse horizontally, letting np.nanmean handle missing data natively
    for i, depth in enumerate(unique_depths):
        mask = (depths_int == depth)
        depth_slice = data[mask, :]
        
        # Only average if the slice exists and isn't 100% NaN (all coplanar channels dead)
        if np.any(mask) and not np.all(np.isnan(depth_slice)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                virtual_data[i, :] = np.nanmean(depth_slice, axis=0)
        else:
            virtual_data[i, :] = np.nan
            
    # 4. Patch any depths where all channels were marked as bad
    df_virt = pd.DataFrame(virtual_data)
    virtual_data = df_virt.interpolate(method='linear', axis=0, limit_direction='both').values
    #why was this interpolation done with pandas? this may be bloating the process surely there is an easier np only way

    
    return virtual_data, unique_depths

def get_depth_map(ks_dir):
    """
    uses ks_dir to return depth_map of channels indexed by id such as depth_map[cid] = depth
    used to sort channels by depth
    returns depthmap
    """
    chan_map = np.load(os.path.join(ks_dir, 'channel_map.npy')).flatten()
    chan_pos = np.load(os.path.join(ks_dir, 'channel_positions.npy'))
    depth_map = np.full(384, np.nan)
    for i, cid in enumerate(chan_map):
        if cid < 384: 
            depth_map[cid] = chan_pos[i, 1]
    return depth_map
    
def get_depth_map_full(ks_dir):
    """
    uses ks_dir to return depth_map of channels indexed by id such as depth_map[cid] = depth
    used to sort channels by depth
    returns depthmap
    """
    chan_map = np.load(os.path.join(ks_dir, 'channel_map.npy')).flatten()
    chan_pos = np.load(os.path.join(ks_dir, 'channel_positions.npy'))
    depth_map = np.zeros((384, 2))
    for i, cid in enumerate(chan_map):
        if cid < 384: 
            depth_map[cid] = chan_pos[i]
    return depth_map



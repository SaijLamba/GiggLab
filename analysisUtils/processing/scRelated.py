import numpy as np
from bg_atlasapi import BrainGlobeAtlas
import pandas as pd

def find_sc_channel(data, params, depth, depth_axis):
    """
    function to find channel id by depth, uses depth from virtual probe
    takes calculated depth from main pipeline find channels within a 100 micron bin around that depth
    isolates data with said channels and time mask of window
    finds most negative
    returns channel id of most negative. Essentially unpacks the virtual probe so instead of finding best depth of sSc (as per the sSc_sink_depth) variable
    it finds best channel around this identified depth and reports the depth of the probe (i.e unpacking the laminarity)
    """
    #data not at "0" as lfp averaging takes data before trigger, need to offset window
    offset = params["window_lfp"][0]
    window = params["window_visual"]
    start = window[0] + abs(offset)
    end = window[1] + abs(offset)
    fs = params["fs_lfp"]
    
    start = int(start * fs)
    end = int(end * fs)
    
    d_max = depth + 100
    d_min = depth - 100
    depth_mask = (depth_axis >= d_min) & (depth_axis <= d_max)
    subset = data[depth_mask, start:end]
    
    index_min = np.nanargmin(subset)
    row, col = np.unravel_index(index_min, subset.shape)

    valid_ids = np.where(depth_mask)[0]

    best_channel = valid_ids[row]
    best_depth = depth_axis[best_channel]
    
    best_time = (start + col) / fs + offset
    volt = subset[row, col]
    print(f"best channel ID {best_channel} @ {best_depth} @ {best_time}ms")
    print(f"nVEP is {volt}")
    return best_channel, best_depth, best_time, volt

    
def map_probe_sc_centered(channel_positions, best_channel_idx, 
                          entry_coords, terminus_coords, 
                          atlas_name='allen_mouse_25um', resolution_step=10):
    
    # --- 1. SETUP ---
    channel_x = channel_positions[:, 0] - 758
    channel_depths = channel_positions[:, 1]
    
    atlas = BrainGlobeAtlas(atlas_name)
    
    # --- 2. TRAJECTORY SAMPLING ---
    p0 = np.array(entry_coords, dtype=float)
    p1 = np.array(terminus_coords, dtype=float)
    traj_vec = p1 - p0
    traj_len = np.linalg.norm(traj_vec)
    traj_unit = traj_vec / traj_len
    
    steps = np.arange(0, traj_len, resolution_step)
    sample_points = p0 + (steps[:, np.newaxis] * traj_unit)
    
    # --- 3. FIND SUPERFICIAL SC INTERSECTION ---
    scs_indices = []
    
    # Define the target ID for the superficial/sensory SC
    # In Allen CCF, 'SCs' (Superior colliculus, sensory related) is usually ID 302.
    # You can also pass a list of IDs if targeting specific layers like SCop, SCsg, SCzo
    superficial_target_ids = [302] 
    
    # Pre-fetch hierarchy for speed
    for i, pt in enumerate(sample_points):
        try:
            sid = atlas.structure_from_coords(pt, microns=True)
            if sid > 0:
                # Check ancestry: Is this point inside ANY of the superficial target structures?
                is_superficial = any(
                    sid == target_id or atlas.hierarchy.is_ancestor(target_id, sid) 
                    for target_id in superficial_target_ids
                )
                
                if is_superficial:
                    scs_indices.append(i)
        except IndexError:
            continue
    
    if not scs_indices:
        raise ValueError("Trajectory does not pass through the Superficial Superior Colliculus.")
    
    # Calculate Superficial SC Geometry
    scs_entry = sample_points[scs_indices[0]]
    scs_exit = sample_points[scs_indices[-1]]
    

    # Geometrical center of the superficial segment
    scs_center_coords = (scs_entry + scs_exit) / 2.0
    scs_span = np.linalg.norm(scs_exit - scs_entry)
    
    print(f"Superficial SC Intersection: Span {scs_span:.0f} um")
    print(f"Superficial SC Center Coords: {scs_center_coords}")
    
    # --- 4. MAP PROBE GEOMETRY ---

    facing_vec = np.array([-1, 0, 0]) 
    lateral_vec = np.cross(traj_unit, facing_vec)
    lateral_vec /= np.linalg.norm(lateral_vec)
    #lateral_vec = -lateral_vec
    
    anchor_depth = channel_depths[best_channel_idx]
    dist_from_anchor = channel_depths - anchor_depth
    
    channel_coords = (scs_center_coords 
                      - (dist_from_anchor[:, np.newaxis] * traj_unit)
                      + (channel_x[:, np.newaxis] * lateral_vec))
    """
    
    target_ml_axis = np.array([0.0, 0.0, -1.0]) 
    
    # 2. Gram-Schmidt Orthogonalization to build a mathematically perfect lateral vector
    dot_product = np.dot(target_ml_axis, traj_unit)
    lateral_vec = target_ml_axis - (dot_product * traj_unit)
    lateral_vec /= np.linalg.norm(lateral_vec)
    
    # 3. Calculate channel distances from the top of the track (or your chosen anchor)
    dist_from_top = channel_depths - channel_depths[0]
    
    # 4. Generate the final 3D coordinates for all channels
    channel_coords = (entry_coords
                      + (dist_from_top[:, np.newaxis] * traj_unit) 
                      + (channel_x[:, np.newaxis] * lateral_vec))

    """
    
    # --- 5. ATLAS LOOKUP FOR CHANNELS ---
    regions = []
    acronyms = []
    
    for coord in channel_coords:
        try:
            sid = atlas.structure_from_coords(coord, microns=True)
            if sid > 0:
                acr = atlas.structures[sid]['acronym']
            else:
                sid = 0
                acr = "Void"
        except IndexError:
            sid = 0
            acr = "Out"
        regions.append(sid)
        acronyms.append(acr)
        
    # --- 6. OUTPUT ---
    df = pd.DataFrame({
        'channel_id': np.arange(len(channel_positions)),
        'rel_x': channel_x,
        'rel_y': channel_depths,
        'ccf_ap': channel_coords[:, 0],
        'ccf_dv': channel_coords[:, 1],
        'ccf_ml': channel_coords[:, 2],
        'region_id': regions,
        'acronym': acronyms
    })
    
    sc_info = {
        'entry': scs_entry,
        'exit': scs_exit,
        'center': scs_center_coords,
        'span_um': scs_span
    }
    
    return df, sc_info

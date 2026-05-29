import os
import numpy as np
import pandas as pd


from .pipeline import run_sc_pipeline
from .helper import get_depth_map_full
from .processing.scRelated import find_sc_channel, map_probe_sc_centered
from .visualisation.funcs import plot_channel_sem, save_pub_fig




def run_localisation(mouse, probe_path, recording_path, ks_dir, triggers, checkerboard_output_dir, ccf_output_dir, results_path, params, force_compute, group):
    """
    runs the whole probe localisation pipeline for a particular mouse - needs the mouse number and 
    1. preprocesses checkerboard neural data -> saves raw + preprocessed plots and channel detection plots
    2. averages by checkerboard -> saves plots
    3. checks if probe trajectory exists -> skips if not
    4. if exits marries probe trajectory and ephys together to get ccf coordinates of all channels -> outputs as csv of channle for mouse
    5. updates and saves .csv of checkerboard results - just min, latency and depth features of the response for each mouse
    """

    results = run_sc_pipeline(recording_path, ks_dir, triggers, mouse, params, save_dir = checkerboard_output_dir, force_compute = force_compute)
    depths = get_depth_map_full(ks_dir)
    scID, scDepth, latency, volt = find_sc_channel(results["lfp"], params, results["anatomy"], results["depths"])

    fig, ax = plot_channel_sem(scID, results["lfp"], results["SEM"])
    save_pub_fig(fig, "10_Best_Channel", mouse, checkerboard_output_dir)

    has_trajectory = os.path.exists(probe_path)

    if has_trajectory:
        probe = np.load(probe_path)
        probe *= 25

        df, info = map_probe_sc_centered(depths, scID, probe[0], probe[1])
        df.to_csv(ccf_output_dir + rf"\Mouse {mouse}.csv", index = False)
    else:
        print("No probe trajectory found, skipping mapping to ccf")

    checkerboard_data = pd.DataFrame([{
        "Mouse": mouse,
        "Sc Channel": scID,
        "Depth (µm)": scDepth,
        "Latency (s)": latency,
        "nVEP (µV)": volt,
        "Treatment Group": group
    }])

    if os.path.exists(results_path):
        existing_df = pd.read_csv(results_path)
        existing_df = existing_df[existing_df["Mouse"] != mouse]
        updated_df = pd.concat([existing_df, checkerboard_data], ignore_index = True)
    else:
        updated_df = new_row

    updated_df = updated_df.sort_values(by="Mouse").reset_index(drop=True)
    updated_df.to_csv(results_path, index=False)
    print(f"Checkerboard results saved/updated to csv at {results_path}")
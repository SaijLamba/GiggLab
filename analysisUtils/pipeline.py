import os
import warnings
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable

from scipy import signal, stats
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde


from .helper import get_depth_map, create_virtual_probe
from .preprocessing.processor import NeuropixelsProcessor
from .storage.funcs import load_scaled_lfp_from_dat, save_scaled_lfp_to_dat
from .visualisation.funcs import save_pub_fig, plot_lfp_depth_profile, plot_channel_diagnostics
from .processing.funcs import calculate_spiking_mua_depth, avg_lfp_by_trigger

def run_sc_pipeline(file_path, ks_dir, triggers, mouse, PARAMS, save_dir=None, force_compute = False):
    print(f"--- PROCESSING: {os.path.basename(file_path)} ---")
    
    # 1. LOAD GEOMETRY
    print("Step 1: Loading Probe Geometry...")
    depth_map = get_depth_map(ks_dir)
    
    
    # 2. CALIBRATION
    print("Step 2: Calibrating (Detecting Bad Channels)...")
    processor = NeuropixelsProcessor(384, PARAMS)
    mm = np.memmap(file_path, dtype='int16', mode='r', pr
                   shape=(os.path.getsize(file_path)//(385*2), 385))
    
    calib_start = 30000 * 60 
    calib_chunk = mm[calib_start:calib_start+60000, :384].T
    calib_stats = processor.detect_bad_channels(calib_chunk)
    print(f"  > Flagged {np.sum(calib_stats['bad_mask'])} Bad Channels")

    
    # 3. SAMPLE VISUALIZATION (Defining sample_raw!)
    print("Step 3: Extracting Diagnostic Sample...")
    # Grab 1 second of data near the first trigger
    samp_start = int(triggers[0])
    samp_end = samp_start + 30000 
    samp_chunk = mm[samp_start:samp_end, :384].T
    
    # THIS IS WHERE sample_raw IS DEFINED
    sample_raw = processor.process_chunk(samp_chunk, clean=False)
    sample_clean = processor.process_chunk(samp_chunk, clean=True)
    
    sample_virtual, depth_axis_sample = create_virtual_probe(sample_clean, depth_map)

    
    # 4. TRIAL AVERAGING
    mean_path =save_dir + r"\data\meanLFP\meanLFP.dat"
    if os.path.isfile(mean_path) and not force_compute:
        print("Mean already computed, loading in as memmap")
        mean_clean, mean_meta = load_scaled_lfp_from_dat(mean_path)
        mean_sem = np.load(save_dir + r"\data\meanLFP\meanLFPSEM.np.npy")
    else:
        print(f"Step 4: Averaging {len(triggers)} Trials...")
        mean_clean, mean_sem = avg_lfp_by_trigger(mm, triggers, PARAMS, processor)
        mean_dir = save_dir + r"\data\meanLFP\meanLFP"
        save_scaled_lfp_to_dat(mean_clean, mean_dir, PARAMS['fs_lfp'])    
        np.save(save_dir + r"\data\meanLFP\meanLFPSEM.np", mean_sem)
        
    # 5. VIRTUAL PROBE (AVERAGE)
    print("Step 5: Creating Virtual Probe (Average)...")
    virt_lfp_final, depth_axis = create_virtual_probe(mean_clean, depth_map )
    mean_virtual_dir = save_dir + r"\data\mean_virtualLFP\mean_virtualLFP"
    save_scaled_lfp_to_dat(virt_lfp_final, mean_virtual_dir, PARAMS['fs_lfp'])

    
    # 6. ANATOMY & CSD
    print("Step 6: Calculating CSD...")
    time_ax = np.linspace(*PARAMS['window_lfp'], mean_clean.shape[1])
    t_mask = (time_ax >= PARAMS['window_visual'][0]) & (time_ax <= PARAMS['window_visual'][1])
    min_profile = np.min(virt_lfp_final[:, t_mask], axis=1)
    depth_mask = depth_axis > 2000 
    valid_depths = depth_axis[depth_mask]
    valid_profile = min_profile[depth_mask]
    ssc_sink_depth = valid_depths[np.argmin(valid_profile)]
    
    lfp_smooth = gaussian_filter(virt_lfp_final, sigma=[5, 0])
    csd = -1 * np.diff(lfp_smooth, n=2, axis=0)
    depth_csd = depth_axis[1:-1]
    print(f"  > sSC Center: {ssc_sink_depth:.1f} um")

    #csd_dir = save_dir + r"\csd"
    #save_scaled_lfp_to_dat(csd, csd_dir)
    

    #Calculate psd
    signal_length = virt_lfp_final.shape[1]  # Assuming axis 1 is time; should equal 600
    
    # Dividing by 2 means a 300-sample window, allowing for roughly 3 overlapping segments
    nperseg_adjusted = signal_length // 2 
    
    # Calculate PSD with the adjusted parameters
    fs = PARAMS['fs_lfp']
    frequencies, psd = signal.welch(
        virt_lfp_final, 
        fs=fs, 
        window='hann', 
        nperseg=nperseg_adjusted,      
        noverlap=nperseg_adjusted // 2,  
        axis=1           
    )

    
    # 7. SPIKES
    print("Step 7: Integrating Spikes...")
    st = np.load(os.path.join(ks_dir, 'spike_times.npy')).flatten()
    sc = np.load(os.path.join(ks_dir, 'spike_clusters.npy')).flatten()
    sp = np.load(os.path.join(ks_dir, 'spike_positions.npy'))
    
    df_units = pd.DataFrame({'cid': sc, 'depth': sp[:, 1]})
    unit_depths = df_units.groupby('cid')['depth'].mean()
    
    unit_stats = []
    for cid in unit_depths.index:
        u_spikes = st[sc==cid]
        #if len(u_spikes) < 20: continue 
        
        counts_pre, counts_post = [], []
        for trig in triggers:
            n_pre = np.searchsorted(u_spikes, trig + PARAMS['spike_win_pre'][1]*PARAMS['fs_raw']) - \
                    np.searchsorted(u_spikes, trig + PARAMS['spike_win_pre'][0]*PARAMS['fs_raw'])
            n_post = np.searchsorted(u_spikes, trig + PARAMS['spike_win_post'][1]*PARAMS['fs_raw']) - \
                     np.searchsorted(u_spikes, trig + PARAMS['spike_win_post'][0]*PARAMS['fs_raw'])
            counts_pre.append(n_pre)
            counts_post.append(n_post)
        
        counts_pre, counts_post = np.array(counts_pre), np.array(counts_post)
        
        if np.all(counts_pre == counts_post) or (np.sum(counts_pre)+np.sum(counts_post)==0):
            p_val = 1.0
        else:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    _, p_val = stats.wilcoxon(counts_pre, counts_post)
            except ValueError: p_val = 1.0
        
        m_pre = np.mean(counts_pre) / 0.2
        m_post = np.mean(counts_post) / 0.2
        mod = (m_post - m_pre)/(m_post + m_pre + 1e-9)
        
        unit_stats.append({'cluster_id': cid, 'depth': unit_depths[cid], 'p_value': p_val, 'modulation': mod, 'sig': p_val<0.05})
        
    df_res = pd.DataFrame(unit_stats)
    
    
    
    # --- VISUALIZATION & SAVING LOGIC ---
    print("Generating Standalone Plots...")

    
    # ==========================================
    # Plot 1: Continuous Raw LFP
    # ==========================================
    fig1, ax0 = plt.subplots(figsize=(8, 6))
    sorted_raw = sample_raw[np.argsort(depth_map), :] 
    vm_raw = np.percentile(np.abs(sorted_raw), 99)
    
    im0 = ax0.imshow(sorted_raw, aspect='auto', cmap='RdBu_r', origin='lower', 
                     vmin=-vm_raw, vmax=vm_raw,
                     extent=[0, 1000, np.min(depth_map), np.max(depth_map)])
    ax0.set_title("Continuous Raw Data (1 sec)")
    ax0.set_xlabel("Time (ms)")
    ax0.set_ylabel("Depth (µm)")
    fig1.colorbar(im0, ax=ax0, label="Voltage (µV)")
    save_pub_fig(fig1, "1_Raw_Continuous", mouse, save_dir)
    
    # ==========================================
    # Plot 2: Continuous Virtual Probe
    # ==========================================
    fig2, ax1 = plt.subplots(figsize=(8, 6))
    vm_virt = np.percentile(np.abs(sample_virtual), 99)
    
    im1 = ax1.imshow(sample_virtual, aspect='auto', cmap='RdBu_r', origin='lower', 
                     vmin=-vm_virt, vmax=vm_virt,
                     extent=[0, 1000, depth_axis_sample[0], depth_axis_sample[-1]])
    ax1.set_title("Continuous Virtual Probe Data")
    ax1.set_xlabel("Time (ms)")
    ax1.set_ylabel("Depth (µm)")
    fig2.colorbar(im1, ax=ax1, label="Voltage (µV)")
    save_pub_fig(fig2, "2_Virtual_Continuous", mouse, save_dir)
    
    # ==========================================
    # Plot 3: Averaged Virtual LFP (Evoked)
    # ==========================================
    fig3, ax2 = plt.subplots(figsize=(8, 6))
    vm_avg = np.percentile(np.abs(virt_lfp_final), 99)
    
    im2 = ax2.imshow(virt_lfp_final, aspect='auto', cmap='RdBu_r', origin='lower', 
                     vmin=-vm_avg, vmax=vm_avg,
                     extent=[time_ax[0], time_ax[-1], depth_axis[0], depth_axis[-1]])
    ax2.axhline(ssc_sink_depth, color='cyan', lw=2, linestyle='--', label=f'Sink @ {ssc_sink_depth:.0f}µm')
    ax2.axvline(0, color='Black', lw=2, linestyle='--', label=f'Stimulus onset')
    ax2.set_title("Averaged Evoked LFP")
    ax2.set_xlabel("Time from Checkerboard Reversal (s)")
    ax2.set_ylabel("Depth (µm)")
    #ax2.legend(loc='upper right') legend looks ugly af
    fig3.colorbar(im2, ax=ax2, label="Voltage (µV)")
    save_pub_fig(fig3, "3_Avg_LFP", mouse, save_dir)
    
    # ==========================================
    # Plot 4: Current Source Density (CSD)
    # ==========================================
    fig4, ax3 = plt.subplots(figsize=(8, 6))
    vm_csd = np.percentile(np.abs(csd), 99)
    
    im3 = ax3.imshow(csd, aspect='auto', cmap='viridis', origin='lower', 
                     vmin=-vm_csd, vmax=vm_csd,
                     extent=[time_ax[0], time_ax[-1], depth_csd[0], depth_csd[-1]])
    ax3.axhline(ssc_sink_depth, color='cyan', lw=2, linestyle='--', label='Primary Sink')
    ax3.axvline(0, color='black', lw=2, linestyle='--', label=f'Stimulus onset')
    ax3.set_title("Current Source Density (CSD)")
    ax3.set_xlabel("Time from Checkerboard Reversal (s)")
    ax3.set_ylabel("Depth (µm)")
    #ax3.legend(loc='upper right')
    fig4.colorbar(im3, ax=ax3, label="CSD (µV/mm²)")
    save_pub_fig(fig4, "4_CSD", mouse, save_dir)
    
    # ==========================================
    # Plot 5: Anatomy vs Physiology (Unit Spikes)
    # ==========================================
    fig5, ax4 = plt.subplots(figsize=(6, 10))
    vis_units = df_res[(df_res['sig']) & (df_res['modulation']>0)]
    motor_units = df_res[~df_res['sig']]

    # Scatter plots for units
    ax4.scatter(motor_units['modulation'], motor_units['depth'],  
                c='gray', alpha=0.3, label='Non-Responsive', edgecolors='none')
    ax4.scatter(vis_units['modulation'], vis_units['depth'], 
                c='red', alpha=0.8, label='Visual Units', edgecolors='white', linewidths=0.5)
    
    # Density plotting
    if len(vis_units) > 1:
        y_grid = np.linspace(PARAMS['probe_tip_offset'], np.max(depth_map), 200)
        kde = gaussian_kde(vis_units['depth'], bw_method=0.1)
        dens = kde(y_grid)
        
        # Create a secondary axis for the density plot to keep scales distinct
        ax4_dens = ax4.twiny()
        ax4_dens.plot(dens, y_grid, 'k--', lw=2, label='Visual Unit Density')
        ax4_dens.set_xlabel("Density", color='k')
        ax4_dens.spines['top'].set_visible(False)
    
    # Landmarks
    ax4.axhline(ssc_sink_depth, color='cyan', lw=3, label='LFP Sink (Superficial)')
    
    ax4.set_title("KDE of visually responsive units")
    ax4.set_ylabel("Probe Depth (µm)")
    ax4.set_xlabel("Modulation Index")
    
    # Combining legends from both axes if twinx was used
    lines_1, labels_1 = ax4.get_legend_handles_labels()
    if len(vis_units) > 1:
        lines_2, labels_2 = ax4_dens.get_legend_handles_labels()
        ax4_dens.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower left', facecolor = "lightgrey", framealpha = 1)
    else:
        ax4.legend(loc='lower left',  facecolor = "lightgrey", framealpha = 1, zorder = 100)
    
    save_pub_fig(fig5, "5_KDE_spikes", mouse, save_dir)



    # ==========================================
    # Plot 6: Spiking MUA (Spike Density Map)
    # ==========================================
    mua_map, t_ax, d_ax = calculate_spiking_mua_depth(st, sp, triggers)
    fig6, ax5 = plt.subplots(figsize=(8, 6))
    
    vm_spiking = np.percentile(mua_map, 99) # Cap the maximum color at the 99.5th percentile to prevent blowout from massive individual units
    
    im5 = ax5.imshow(mua_map, aspect='auto', cmap='magma', origin='lower', 
                     vmin=0, vmax=vm_spiking,
                     extent=[t_ax[0], t_ax[-1], d_ax[0], d_ax[-1]])
    
    # Overlay the previously identified CSD sink for anatomical correlation
    ax5.axhline(ssc_sink_depth, color='cyan', lw=2, linestyle='--', label=f'LFP Sink @ {ssc_sink_depth:.0f}µm')
    ax5.axvline(0, color='white', lw=2, linestyle='--', label=f'Stimulus Onset')
    
    ax5.set_title("Evoked Spiking Multi-Unit Activity (Spike Density)")
    ax5.set_xlabel("Time from Checkerboard Reversal (ms)")
    ax5.set_ylabel("Probe Depth (µm)")
    #ax5.legend(loc='upper right')
    fig6.colorbar(im5, ax=ax5, label="Firing Rate (Spikes/s)")
    
    save_pub_fig(fig6, "6_Spiking_MUA_Map", mouse, save_dir)
    # ==========================================
    # Plot 7: psd
    # ==========================================
    
    freq_mask = (frequencies >= 1) & (frequencies <= 100)
    freqs_plot = frequencies[freq_mask]
    psd_plot = psd[:, freq_mask]

    
    fig7, ax6 = plt.subplots(figsize=(8, 6))
    mesh = ax6.pcolormesh(
        freqs_plot, 
        depth_axis, 
        psd_plot, 
        shading='auto', 
        cmap='viridis',
        norm=LogNorm() # Log-scale normalization for power
    )
    

    ax6.set_xlabel('Frequency (Hz)')
    ax6.set_ylabel('Depth ($\mu m$)')
    ax6.set_title('Local Field Potential Power Spectral Density')
    ax6.axhline(ssc_sink_depth, color='cyan', lw=2, linestyle='--', label=f'LFP Sink @ {ssc_sink_depth:.0f}µm')
    #ax6.legend(loc='upper right')
    cbar = fig7.colorbar(mesh, ax=ax6, label='Power ($V^2/Hz$)')
    save_pub_fig(fig7, "7_PSD", mouse, save_dir)


    # ==========================================
    # Waveforms all channels
    # ==========================================
    
    #Waveforms all channels
    num_channels, num_samples = mean_clean.shape
    time_vector = np.linspace(-0.1, 0.5, num_samples)
    depths = depth_map
    # 4. Set up the colormap based on depth
    # We normalize the depth values to a range of 0 to 1 for the colormap
    norm = Normalize(vmin=np.min(depths), vmax=np.max(depth_axis))
    cmap = plt.get_cmap('coolwarm') # 'viridis', 'plasma', or 'coolwarm' work well
    
    # 5. Create the plot
    fig8, ax7 = plt.subplots(figsize=(10, 8))

    # 6. Loop through each channel and plot
    for i in range(num_channels):
        # Get the depth for this specific channel
        channel_depth = depths[i]
        
        # Get the specific color for this depth
        color = cmap(norm(channel_depth))
        
        # Plot the waveform
        ax7.plot(time_vector, mean_clean[i, :], color=color, alpha=0.7, linewidth=0.8)
    
    # 7. Add a vertical line at t=0 (stimulus onset)
    ax7.axvline(x=0, color='red', linestyle='--', label='Stimulus Onset')
    
    # 8. Formatting the axes and labels
    ax7.set_title('Mean LFP Waveforms Across Channels', fontsize=14)
    ax7.set_xlabel('Time from Checkerboard Reversal (s)', fontsize=12)
    ax7.set_ylabel('Voltage (\u00B5V)', fontsize=12)
    ax7.set_xlim([-0.1, 0.5])
    
    # 9. Add a colorbar to act as a depth legend
    # We create a ScalarMappable to link the colormap and the normalization
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([]) # Required for older versions of matplotlib
    cbar = fig8.colorbar(sm, ax=ax7)
    cbar.set_label('Probe Depth (\u00B5m)', rotation=270, labelpad=15)
    save_pub_fig(fig8, "8_Waveforms", mouse, save_dir)

    fig11, ax10 = plot_lfp_depth_profile(mean_clean, time_vector, depths, np.max(depth_axis), scale_x_len=0.1, scale_y_len=2)
    save_pub_fig(fig11, "11_Waveforms_bar", mouse, save_dir)


    # ==========================================
    # virtual Waveforms all channels
    # ==========================================
    #Virtual LFP waveforms by depth
    num_channels, num_samples = virt_lfp_final.shape
    time_vector = np.linspace(-0.1, 0.5, num_samples)
    depths = depth_axis
    # 4. Set up the colormap based on depth
    # We normalize the depth values to a range of 0 to 1 for the colormap
    norm = Normalize(vmin=np.min(depths), vmax=np.max(depth_axis))
    cmap = plt.get_cmap('coolwarm') # 'viridis', 'plasma', or 'coolwarm' work well
    
    # 5. Create the plot
    fig9, ax8 = plt.subplots(figsize=(10, 8))
    
    # 6. Loop through each channel and plot
    for i in range(num_channels):
        # Get the depth for this specific channel
        channel_depth = depths[i]
        
        # Get the specific color for this depth
        color = cmap(norm(channel_depth))
        
        # Plot the waveform
        ax8.plot(time_vector, virt_lfp_final[i, :], color=color, alpha=0.7, linewidth=0.8)
    
    # 7. Add a vertical line at t=0 (stimulus onset)
    ax8.axvline(x=0, color='red', linestyle='--', label='Stimulus Onset')
    
    # 8. Formatting the axes and labels
    ax8.set_title('Mean LFP Waveforms Across Channels')
    ax8.set_xlabel('Time from Checkerboard Reversal (s)')
    ax8.set_ylabel('Voltage (\u00B5V)')
    ax8.set_xlim([-0.1, 0.5])
    
    # 9. Add a colorbar to act as a depth legend
    # We create a ScalarMappable to link the colormap and the normalization
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([]) # Required for older versions of matplotlib
    cbar = fig9.colorbar(sm, ax=ax8)
    cbar.set_label('Probe Depth (\u00B5m)', rotation=270, labelpad=15)
    save_pub_fig(fig9, "9_Waveforms_virtual", mouse, save_dir)



    #############################################
    #### BAD CHANNELS PLOT######################
    ############################################
    
    
    
    # Extract the necessary arrays and parameters
    bad_mask = calib_stats['bad_mask']
    metrics = calib_stats['metrics']
    thresholds = processor.p['bad_ch_thresh']
    
    print("2. Generating Channel Diagnostics Plot...")
    # Call the standalone plotting function
    fig10, ax9 = plot_channel_diagnostics(
        rms=metrics['rms'], 
        corrs=metrics['corrs'], 
        median_rms=metrics['median_rms'], 
        bad_mask=bad_mask, 
        thresholds=thresholds
    )

    save_pub_fig(fig10, "10_Bad_Channel_Detection", mouse, save_dir)


    
    
    return {'lfp_virt': virt_lfp_final, 
            'csd': csd, 
            'units': df_res, 
            'anatomy': ssc_sink_depth, 
            "depth_axis": depth_axis, 
            "lfp": mean_clean, 
            "depths": depth_map, 
            "SEM": mean_sem}

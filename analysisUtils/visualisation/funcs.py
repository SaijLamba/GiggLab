import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from matplotlib.cm import ScalarMappable

def save_pub_fig(fig, suffix, mouse, save_dir):
    """
    Saves a figure in both high-res PNG and vector PDF formats.
    """
    save_dir = save_dir + r"\plots"
    if not save_dir:
        return
        
    os.makedirs(save_dir, exist_ok=True)
    base_name = f"Mouse{mouse}"
    fname = f"{base_name}_{suffix}"
    
    png_path = os.path.join(save_dir, f"{fname}.png")
    pdf_path = os.path.join(save_dir, f"{fname}.pdf")
    
    # bbox_inches='tight' prevents axis labels from being cut off
    fig.savefig(png_path, dpi=300, bbox_inches='tight', transparent=False)
    fig.savefig(pdf_path, format='pdf', bbox_inches='tight', transparent=True)
    print(f"Saved: {fname}")
    plt.close(fig) # Frees memory
    

def plot_channel_diagnostics(rms, corrs, median_rms, bad_mask, thresholds):
    """
    Generates a scatter plot of channel Correlation vs. RMS for quality control.
    
    Parameters:
        rms (np.ndarray): 1D array of RMS values for each channel.
        corrs (np.ndarray): 1D array of Pearson correlation coefficients to the global median.
        median_rms (float): The median RMS value across all channels.
        bad_mask (np.ndarray): 1D boolean array where True indicates a bad channel.
        thresholds (dict): Dictionary containing 'rms_low', 'rms_high_factor', and 'corr_min'.
        
    Returns:
        fig (matplotlib.figure.Figure): The figure object.
        ax (matplotlib.axes.Axes): The axes object containing the plot.
    """
    # 1. Initialize the figure
    fig, ax = plt.subplots(figsize=(10, 6))

    # 2. Extract threshold boundaries from the passed dictionary
    rms_low = thresholds['rms_low']
    rms_high = median_rms * thresholds['rms_high_factor']
    corr_min = thresholds['corr_min']

    # 3. Create the inverse mask for valid channels
    good_mask = ~bad_mask

    # 4. Plot Good and Bad Channels as distinct scatter series
    ax.scatter(corrs[good_mask], rms[good_mask], 
               color='#2ca02c', label='Good Channels', alpha=0.7, edgecolors='w', linewidth=0.5)
    ax.scatter(corrs[bad_mask], rms[bad_mask], 
               color='#d62728', label='Bad Channels', alpha=0.7, edgecolors='w', linewidth=0.5)

    # 5. Render Threshold Lines
    ax.axhline(rms_low, color='black', linestyle='--', linewidth=1.5, 
               label=f'RMS Low ({rms_low:.2f})')
    ax.axhline(rms_high, color='black', linestyle='-.', linewidth=1.5, 
               label=f'RMS High Factor ({thresholds["rms_high_factor"]}x Med)')
    ax.axvline(corr_min, color='black', linestyle=':', linewidth=1.5, 
               label=f'Corr Min ({corr_min:.2f})')

    # 6. Formatting and Labeling
    ax.set_xlabel('Pearson Correlation to Global Median')
    ax.set_ylabel('RMS Voltage')
    ax.set_title('Channel Quality: RMS vs. Common Median Correlation')
    
    # 7. Dynamic Shading of the Acceptance Region
    # We calculate the normalized y-coordinates (0.0 to 1.0) relative to the current axes limits
    y_min_norm = (rms_low - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    y_max_norm = (rms_high - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    
    ax.axvspan(corr_min, 1.0, ymin=y_min_norm, ymax=y_max_norm, 
               color='green', alpha=0.05, label='Acceptance Region')

    # 8. Final UI Adjustments
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    return fig, ax


def plot_lfp_depth_profile(mean_clean, time_vector, depths, max_depth, scale_x_len=0.1, scale_y_len=50):
    """
    Plots mean LFP waveforms colored by probe depth, using a clean scale bar 
    instead of standard axes.
    
    Parameters:
        mean_clean (np.ndarray): 2D array of waveforms (channels x samples).
        time_vector (np.ndarray): 1D array of time points in seconds.
        depths (np.ndarray): 1D array of depths for each channel.
        max_depth (float): The maximum depth value for colormap normalization.
        scale_x_len (float): Length of the horizontal scale bar in seconds.
        scale_y_len (float): Length of the vertical scale bar in µV.
        
    Returns:
        fig (matplotlib.figure.Figure): The figure object.
        ax (matplotlib.axes.Axes): The axes object.
    """
    num_channels = mean_clean.shape[0]
    
    # 1. Initialize the figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 2. Set up the colormap based on depth
    norm = Normalize(vmin=np.min(depths), vmax=max_depth)
    cmap = plt.get_cmap('coolwarm')
    
    # 3. Plot each channel's waveform
    for i in range(num_channels):
        color = cmap(norm(depths[i]))
        ax.plot(time_vector, mean_clean[i, :], color=color, alpha=0.7, linewidth=0.8)
        
    # 4. Add stimulus onset marker
    ax.axvline(x=0, color='red', linestyle='--', linewidth=1.5, label='Stimulus Onset')
    
    # 5. Strip the standard axes (spines and ticks)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # 6. Construct the L-shaped Scale Bar
    # Anchor the scale bar to the bottom left of the data extent
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # Calculate a dynamic 5% margin for clean padding
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    
    # Offset slightly inwards from the bottom-left corner
    anchor_x = x_min + x_margin
    anchor_y = y_min + y_margin
    
    # Draw Horizontal Bar (Time) - draws left to right
    ax.plot([anchor_x, anchor_x + scale_x_len], [anchor_y, anchor_y], 
            color='black', linewidth=5, solid_capstyle='butt')
            
    # Text for Horizontal Bar - centered below the bar
    ax.text(anchor_x + (scale_x_len / 2), anchor_y - (y_margin * 0.4), 
            f"{scale_x_len} s", ha='center', va='top')
            
    # Draw Vertical Bar (Voltage) - draws bottom to top
    ax.plot([anchor_x, anchor_x], [anchor_y, anchor_y + scale_y_len], 
            color='black', linewidth=5, solid_capstyle='butt')
            
    # Text for Vertical Bar - placed to the left of the bar
    ax.text(anchor_x - (x_margin * 0.4), anchor_y + (scale_y_len / 2), 
            f"{scale_y_len} \u00B5V", ha='right', va='center', rotation=90)
                

    # 7. Add Depth Colorbar Legend
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([]) 
    cbar = fig.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label('Probe Depth (µm)', rotation=270, labelpad=20)
    
    # 8. Final touches
    ax.set_title('Mean LFP Waveforms Across Channels', pad=15)
    
    return fig, ax

    
def plot_channel_sem(target_channel, meanlfp, semlfp):
    
    num_channels, num_samples = meanlfp.shape
    
    # 2. Generate the time vector
    # 600 samples at 1000 Hz spans exactly 0.6 seconds.
    time_vector = np.linspace(-0.1, 0.5, num_samples)
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot the central mean waveform
    ax.plot(time_vector, meanlfp[target_channel, :], color='black', alpha=1, linewidth=1.5, label='Mean LFP')
    
    # Plot the SEM as a shaded region using the loaded semlfp array
    ax.fill_between(time_vector, 
                     meanlfp[target_channel, :] - semlfp[target_channel, :], 
                     meanlfp[target_channel, :] + semlfp[target_channel, :], 
                     color='gray', alpha=0.3, label='$\pm$ 1 SEM')
    
    # Add stimulus onset marker
    ax.axvline(x=0, color='red', linestyle='--', label='Checkerboard Onset')
    
    # Formatting the axes and labels
    ax.set_title(f'Mean LFP Waveform (Channel {target_channel})', fontsize=14)
    ax.set_xlabel('Time relative to stimulus (s)', fontsize=12)
    ax.set_ylabel('Voltage (\u00B5V)', fontsize=12)
    ax.set_xlim([-0.1, 0.5])
    
    ax.legend(loc='upper right')
    return fig, ax



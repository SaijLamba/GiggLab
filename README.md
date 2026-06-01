# Read Me
This repo is part of a Msci project looking at probe localisation with in vivo electrophysiology with Neuropixels 2.0.
This is done in the context of a broader experiment looking at the anxiolytic properties of Psilocybin by Micheal Hogan

Current state of the repo is very messy and is going through a cleanup phase

The functions were built within the jupyter notebooks themseles for ease of protoyping and using, however as the code and pipeline is essentially near its "finished" form the functions will become a package to use
The hope is to use this as a basis for many e-phys functions, and as such ton generalise these functions to use outside of ephys
- This may occur as a seperate repo, to keep the (cleaned) code used for this experiment and code turned into package distinct

Saij Lamba

# Installation
It is recommended to install this using a virtual environment to avoid any dependency clashes, the requirements.txt file contains the dependencies needed and can be installed via your favourite package/virtual environment manager

Currently the package is not available on pip (hopefully it will be) so to use the packages you need to
1. Clone the repo
2. Put the analysisUtils folder within your working directory to utilise it in your python scripts

Eventually the package will contain setup parameters so you may pip install it from the repo even if it isn't available on pip, but hopefully it will be on pip


# Usage
## Histology/reference atlas registration

This project utilises [HERBS](https://github.com/Whitlock-Group/HERBS) to register brain slices to the allen mouse atlas and track probe trajectories in the common ccf
Please refer to the HERBS cookbook for a far better tutorial on how to do so than I can write
The end goal is to find the vector trajectory of the probe within CCF, this can be done with HERBS or whatever method you wish. The rest of the pipeline utilises the np array of the probe vector in the format. We store the vector as voxels and load and scale to microns as needed, downstream analysis functions take the probe vectors as microns

``` python
probe = [insertion_coords, terminus_coords]
```

In theory any two points on the line which contains the probe can be utilised, however this is untested with the code base and is discouraged. To go from HERBS pkl file to a np file containing insertion and terminus coords used by the rest of the code:

```python
from analysisUtils.histology.probetrajec import pkl2traj

probe_dir = "path to pkl file from HERBS output"
output = "output folder for np array"

pkl2traj(probe_dir, output)

```

## Ephys (pre)processing
This describes the core driver of the analysis uses a custome written processor class that can be initalised with various recording parameters
This class (and relevant methods) is utlised to take a chunk of ephys data takes a np.array of shape (n_channels, timepoints) to output a cleaned LFP of said data in format (n_channels, time_point)

This is natively implemented in the pipeline functions as such an understanding of this can be skipped, but it is here for the curious reader or if you'd like to reverse engineer the code

The analysis funcs and this processor need a set of parameters to be initialised, recommended paramaters are as follows

 ```python
 PARAMS = {
    # Hardware Specs
    'fs_raw': 30000,          # Acquisition Rate (Hz)
    'fs_lfp': 1000,           # Downsampled Rate (Hz)
    'voltage_scaling': 0.195, # uV per bit scaling, the .dat file saved as 16-bit binary needs to be scaled to voltages
    'probe_tip_offset': 720,  # Microns from tip with no channels

    # Analysis Windows
    'window_lfp': (-0.1, 0.5),    # Time window around trigger to average data around(s)
    'window_visual': (0.075, 0.150),# Window to find sSC sink (75 - 150ms post-stim)
    'bin_size_um': 15,            # Spatial bin for virtual probe, I think this is redundant now but still written to be safe, func does it dynmaically i believe
    
    # Bad Channel Thresholds
    'bad_ch_thresh': {
        'rms_low': 1.0,        # Dead (<1uV)
        'rms_high_factor': 15, # Noisy (>15x median)
        'corr_min': 0.6        # Artifact (<0.6 correlation to brain)
    },
    
    # Spike Statistics
    'spike_win_pre': (-0.2, 0.0), # Baseline window (s)
    'spike_win_post': (0, 0.2)# Visual response window (s)
}
```

With these paramaters defined the processor can be initalised for your specific animal, recording or unique set of paramaters as follows

```python 
from analysisUtils.preprocessing.processor import NeuropixelsProcessor

processor = NeuropixelsProcessor(384, PARAMS) #processor takes argument of n_channels (0 indexed) and params
```
The initalisation of the processor creates attributes such as the numerators and denomitors the butterworth polynomial needs in order to bandpass for LFP (3-300hz) and the shifting needed to deskew the data in the time domain (to account for the multiplexing)

### Bad channel detection 
Next the processor can be used to determine bad channels from the defined params with the following defintions:
1. Any channel that has a correlation less than 0.6 with the probe median is bad
2. Any channel with a RMS of less than 1 is bad (dead channel)
3. Any channel with a RMS of 15 x probe median RMS is bad (inherently noisy channel)

the processor returns the bad channel metrics (which is further utilised downstream to plot) and stores the bad channels as self.bad_mask

Bad channel detection is recommended for a chunk of data pre any stimulus presentation or experimental manipulation i.e the 1st second of recording

```python

bad_channel_stats = processor.detect_bad_channels(calibration_chunk) #calibration chunk is a chunk of ephys data passed into this method to detect bad channels
```

The bad channel detection is not dynamic, i.e it is done once and the mask is then used for the rest of the recording or processing of chunks, with a mask of said channels stored as the bad_mask attribute for the processor

### Data preprocessing
With bad channel detection complete the proccessor can now preprocess the data 
This processor:
1. Scales the chunk to voltage
2. If clean = false then just returns the data bandpass filtered and downsampled
3. if clean = true then:
    1. deskews the data 
    2. sets bad channels to nan (not 0, as technically 0 is still a datapoint) - used bad mask worked out above
    3. detrending and common mode referencing by subtracting probe median across time and channels
    4. apply the butterworth filter to filter between 3 and 300hz and then downsamples as specified by the parameters


```python 
clean_chunk = processor.process_chunk(raw_chunk, clean = True)
```


## Sc pipeline
The processor above is called and utlised in functions related to averaging LFP over checkerboard stimulus and running probe localisation and is neatly wrappep up in a probe localisation function that needs to take arguments:

1. mouse number
2. path to processed probe np array file (see histology) (assumes in voxels, loads and scales to microns assuming a 25micron atlas)
3. path to recording .dat file (assumes unscaled, 16 bit file in shape (timepoints, channels) and transposes to (channels, timepoints) as needed)
4. ks_dir containing relevant spike files and channel maps 
5. triggers -> trigger times (in samples) for checkerboard reversal - flattened numpy array
6. checkerboard output directory, save dir for each mouse to store meanLFP data and plots
7. ccf_output_dir -> output dir to store the channel localisaiton data for each animal
8. results_path -> path to .csv file to store extracted features (min, latency, depth) of checkerboard response for each animal
9. params (see above)
10. force compute - boolean to specify whether to force recomputation of mean instead of checking if it already exists and loading it in
11. group, treatment group information for particular animal

The pipeline: 
1. preprocesses checkerboard neural data -> saves raw + preprocessed plots and channel detection plots
2. averages by checkerboard -> saves plots
3. checks if probe trajectory exists -> skips if not
4. if exits marries probe trajectory and ephys together to get ccf coordinates of all channels -> outputs as csv of channle for mouse
5. updates and saves .csv of checkerboard results - just min, latency and depth features of the response for each mouse

```python
#example code, lets assume mouse 2
from analysisUtils.probelocalisation import run_localisation #the actual pipeline function

mouse = 2
recording_path = "path to mouse 2 .dat file"
ks_dir = "path to ks_dir for mouse 2"
triggers = [numpy array of trigger times]
checkerboard_output_dir = "output directory path for storing plots and saved lfp
ccf_output = "path to store channel information"
results_path = "path to .csv to store extracted features"
probe_path = "path to numpy array of insertion and terminus coords"
treamtent_group = "Vehicle" #or "Psilocybin"

run_localisation(
    mouse,
    probe_path,
    recording_path,
    ks_dir,
    triggers,
    checkerboard_output_dir,
    ccf_output,
    results_path,
    params,
    False,
    treatment_group
)

```

The pipeline utilises [BrainGlobe Atlas Api] (https://brainglobe.info/documentation/brainglobe-atlasapi/index.html) to sample within the allen ccf to find region information for each channel
import os
import json
import numpy as np



def save_scaled_lfp_to_dat(lfp_array, file_path, sampling_rate):
    """
    Saves a voltage-scaled LFP numpy array to a flat binary .dat file,
    safely creating parent directories if they do not exist, and creates
    a companion JSON file.
    """
    
    # 1. Verify data type
    if lfp_array.dtype != np.float32:
        lfp_array = lfp_array.astype(np.float32)
        
    # 2. Extract properties to conserve
    metadata = {
        "dtype": str(lfp_array.dtype),
        "shape": lfp_array.shape,
        "sampling_rate_hz": sampling_rate,
        "is_voltage_scaled": True
    }
    
    # 3. Ensure the parent directory exists (THE FIX)
    parent_directory = os.path.dirname(file_path)
    if parent_directory:
        os.makedirs(parent_directory, exist_ok=True)
    
    # 4. Define file paths (adding .dat extension if missing for clarity)
    if not file_path.endswith('.dat'):
        file_path = f"{file_path}.dat"
        
    base_name = os.path.splitext(file_path)[0]
    json_path = f"{base_name}.json"
    
    # 5. Write the binary data directly to disk
    lfp_array.tofile(file_path)
    
    # 6. Write the metadata companion file
    with open(json_path, 'w') as f:
        json.dump(metadata, f, indent=4)
        
    print(f"Saved LFP binary to {file_path}")
    print(f"Saved LFP metadata to {json_path}")

def load_scaled_lfp_from_dat(file_path):
    # 1. Load metadata to get shape and dtype
    base_name = os.path.splitext(file_path)[0]
    with open(f"{base_name}.json", 'r') as f:
        metadata = json.load(f)
        
    shape = tuple(metadata["shape"])
    dtype = np.dtype(metadata["dtype"])
    
    # 2. Memory-map the flat binary file back into a NumPy array
    lfp_array = np.memmap(file_path, dtype=dtype, mode='r', shape=shape)
    
    return lfp_array, metadata["sampling_rate_hz"]

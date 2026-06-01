import numpy as np
import pickle as pkl

def load_probe(probe_dir):
    """
    helper function to load HERBS probe pickle file
    """
    with open(probe_dir + r"\probe.pkl", "rb") as f:
        probe = pkl.load(f)
    return probe

def calculate_ccf_trajectory(herbs_insertion, herbs_terminus):
    """
    Takes insertion and terminus voxels from HERBS output and converts them to allen CCF voxels
    """
    Y_MAX = 528
    Z_MAX = 320
    

    ccf_ins_ap = Y_MAX - herbs_insertion[1]
    ccf_ins_dv = Z_MAX - herbs_insertion[2]
    ccf_ins_ml = herbs_insertion[0]
    ccf_insertion = np.array([ccf_ins_ap, ccf_ins_dv, ccf_ins_ml])

    ccf_term_ap = Y_MAX - herbs_terminus[1]
    ccf_term_dv = Z_MAX - herbs_terminus[2]
    ccf_term_ml = herbs_terminus[0]
    ccf_terminus = np.array([ccf_term_ap, ccf_term_dv, ccf_term_ml])


    
    return ccf_insertion, ccf_terminus


def pkl2traj(probe_dir, output):
    """
    Wraps above functions into one, taking in input directory of HERBS probe and output directory to save the CCF compatible insertion and terminus coordinate
    """
    probe = load_probe(probe_dir)
    insertion_vox = probe["data"]["insertion_vox"]
    terminus_vox = probe["data"]["terminus_vox"]

    ccf_insert, ccf_terminus = calculate_ccf_trajectory(insertion_vox, terminus_vox)

    probe_vox  = [ccf_insert, ccf_terminus]

    np.save(output + r"\probe_ccf_vox.npy", probe_vox)

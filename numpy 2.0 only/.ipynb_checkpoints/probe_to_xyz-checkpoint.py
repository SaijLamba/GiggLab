import pickle as pkl
import numpy as np
import json



def herbs_to_ibl(probe_dir, output_dir, bins):
    with open(probe_dir, "rb") as f:
        probe = pkl.load(f)

    probeInsertionCoords = np.array(probe["data"]["insertion_coords"])
    probeTerminusCoords = np.array(probe["data"]["terminus_coords"])

    #data already in mlapdv
    probeVector = probeTerminusCoords - probeInsertionCoords
    probeUnitVec = probeVector / bins

    xyz = [probeInsertionCoords]

    i = 1
    while i <= bins:
        coord = [probeInsertionCoords + (probeUnitVec * i)]
        xyz = np.concatenate((xyz, coord), axis = 0)
        i += 1


    xyz_picks = {"xyz_picks": xyz.tolist()}
    with open(output_dir + r'\xyz_picks.json', 'w') as fp:
        json.dump(xyz_picks, fp, indent= 2)


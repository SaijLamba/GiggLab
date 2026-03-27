import numpy as np
from brainrender import Scene
from brainrender.actors import Point, Line
from bg_atlasapi import BrainGlobeAtlas
import bg_atlasapi
atlas_name = "allen_mouse_25um" 
bg_atlas = BrainGlobeAtlas(atlas_name)



probe_path = r"Z:\Saij\Herbs probes\Mouse 7\probe_ccf_vox.npy"
probe = np.load(probe_path) * 25
print(probe)
# 2. Print the spatial metadata
print(f"Atlas Name: {bg_atlas.metadata['name']}")
print(f"Orientation String: {bg_atlas.metadata['orientation']}")
print(f"Resolution (um): {bg_atlas.metadata['resolution']}")
print(f"Shape (voxels): {bg_atlas.metadata['shape']}")


scene = Scene(title="Probe Trajectory Check", atlas_name="allen_mouse_25um")

scene.add_brain_region("SCm", alpha=0.5, color="green")
scene.add_brain_region("PAG", alpha=0.5, color="blue")

insertion_coord = probe[0]
terminus_coord = probe[1]

# Replaced 'Sphere' with 'Point'
insertion_marker = Point(insertion_coord, radius=80, color="red")
scene.add(insertion_marker)

# Replaced 'Cylinder' with 'Line' to guarantee a point-to-point connection
track_coordinates = np.array([insertion_coord, terminus_coord])
track_line = Line(track_coordinates, linewidth=10, color="yellow")
scene.add(track_line)

scene.slice("sagittal")
scene.render()
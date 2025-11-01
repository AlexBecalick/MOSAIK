#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from resegmentation import Resegmentation

fov = 303
channel_names = ['B', 'G', 'Y', 'R', 'U']
# zarr_dir = '/Volumes/SSD'
zarr_dir = '/Users/alvincruiziat/Downloads/London/KCL'
zarr_name = 'G16_morphology_2D.zarr'
seg_path = f'FOV_{fov}_segmentation_ready.png'

pipe = Resegmentation(zarr_dir, zarr_name, fov)
pipe.preprocess_image(channel_names)
# pipe.run_cellpose("path/to/preprocessed.png")
pipe.run_cellpose(seg_path)
new_sdata = pipe.update_spatialdata()
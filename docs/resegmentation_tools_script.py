import os
# path = '/Volumes/SSD/'
path = '/Users/alvincruiziat/Downloads/London/KCL'
os.chdir(path)
from regmentations_tools import check_gpu, preprocess_spatial_image
from regmentations_tools import run_cellpose_SAM, update_spatialdata
from regmentations_tools import process_masks_to_shapes, process_labels, process_points_and_tables

# Confirming system information, TensorFlow version & GPU access
check_gpu()

""" PRE SEGMENTATION: LOADING IMG """
FOV = 303
seg_path = f'FOV_{FOV}_segmentation_ready.png'
channels = ['B', 'G', 'Y', 'R', 'U']
downscale = (768, 768)
# base_dir = ''
base_dir = '/Users/alvincruiziat/Downloads/London/KCL'
zarr_name = 'G16_morphology_2D.zarr'
sdata = preprocess_spatial_image(zarr_dir=base_dir, zarr_name=zarr_name, fov=FOV, 
                                   channel_names=channels, channels_to_use=["Y", "U"], 
                                   output_path=seg_path, thumbnail_size=downscale) 

""" SEGMENTATION: CELLPOSE-SAM """
full_path = os.path.join(base_dir, seg_path)
masks, flows, styles = run_cellpose_SAM(img_path=full_path)

""" POST SEGMENTATION: UPDATE SDATA object"""
# Masks → Shapes
masks, gdf_polygons, shapes_model = process_masks_to_shapes(masks, FOV)

# Labels
masks, labels_model = process_labels(masks, FOV)

# Points + Tables
gdf_points, vdata = process_points_and_tables(sdata, gdf_polygons, FOV)

""" SAVING ZARR FILE"""
update_spatialdata(sdata, vdata, labels_model, shapes_model, base_dir, FOV)
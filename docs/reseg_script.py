import os
path = '/Users/k2481276/Documents/'
os.chdir(path)
from reseg_tools import stitching, vis_segmentation, qc_metric, qc_plot, filter_norm
from reseg import Resegmentation
import spatialdata as sd
import pandas as pd


sdata = sd.read_zarr("/Volumes/SSD/resegmentation/G16_morphology_2D.zarr")
fov_position = pd.read_csv("/Volumes/SSD/resegmentation/G16_fov_positions_file.csv")
template = ['_image', '_labels', '_shapes', '_points']

sample1 = list(range(1, 72 + 1))
S1 = [str(a) + b for a in sample1 for b in template]
sdata_S1 = sdata.subset(element_names=S1, include_orphan_tables=False)
stitching(sdata_S1, fov_position, name="S1")


channel_names = ['B', 'G', 'Y', 'R', 'U']
channels_to_use=['U']
factor_rescale = 8
image_name = 'Stitched'
label_name = 'labels'
shape_name = 'shapes'
point_name = 'points'
output = 'Stiched_segmentation_ready.png'
name = 'S1'
zarr_name = 'stitched_' + name + '.zarr'


pipe = Resegmentation(path, 
                      zarr_name, 
                      output,
                      factor_rescale,
                      image_name,
                      label_name,
                      shape_name,
                      point_name)

pipe.preprocess_image(channel_names, channels_to_use)
pipe.run_cellpose(flow_threshold=1.2,
                  cellprob_threshold=-3,
                  tile_overlap=0.15)
pipe.update_spatialdata()
    

### Start analysis    
sdata_S1 = sd.read_zarr(path + "/updated_stitched_S1.zarr")
vis_segmentation(sdata_S1, "S1_6")

adata_S1 = sdata_S1.tables["table"]   
qc_metric(adata_S1)  
qc_plot(adata_S1, sdata_S1)
filter_norm(adata_S1)

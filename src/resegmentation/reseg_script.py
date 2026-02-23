import os
import sys
from pathlib import Path

# Load parameters from src/params.yaml (section: resegmentation_cosmx)
_src_dir = Path(__file__).resolve().parent.parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))
from params_loader import get_params

p = get_params("resegmentation_cosmx")
path = p["base_path"]
os.chdir(path)

from reseg_tools import stitching, vis_segmentation, qc_metric, qc_plot, filter_norm
from reseg import Resegmentation
import spatialdata as sd
import pandas as pd

sdata = sd.read_zarr(p["zarr_path"])
fov_position = pd.read_csv(p["fov_positions_path"])
template = p["template"]

sample1 = list(range(p["sample_range_start"], p["sample_range_end"] + 1))
name = p["sample_name"]
S1 = [str(a) + b for a in sample1 for b in template]
sdata_S1 = sdata.subset(element_names=S1, include_orphan_tables=False)
stitching(sdata_S1, fov_position, name=name)

channel_names = p["channel_names"]
channels_to_use = p["channels_to_use"]
factor_rescale = p["factor_rescale"]
image_name = p["image_name"]
label_name = p["label_name"]
shape_name = p["shape_name"]
point_name = p["point_name"]
output = p["preprocessed_image_output"]
zarr_name = p["stitched_zarr_name"]

pipe = Resegmentation(
    path,
    zarr_name,
    output,
    factor_rescale,
    image_name,
    label_name,
    shape_name,
    point_name,
)

pipe.preprocess_image(channel_names, channels_to_use)
pipe.run_cellpose(
    flow_threshold=p["flow_threshold"],
    cellprob_threshold=p["cellprob_threshold"],
    tile_overlap=p["tile_overlap"],
)
pipe.update_spatialdata()

### Start analysis
updated_zarr = p["updated_zarr_name"]
updated_path = os.path.join(path, updated_zarr) if not os.path.isabs(updated_zarr) else updated_zarr
sdata_S1 = sd.read_zarr(updated_path)
vis_segmentation(sdata_S1, p["vis_element_name"])

adata_S1 = sdata_S1.tables["table"]   
qc_metric(adata_S1)  
qc_plot(adata_S1, sdata_S1)
filter_norm(adata_S1)

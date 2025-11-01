Extracting Files
================

Before initializing a ``SpatialData`` object, both CosMx and Xenium datasets require locating and preparing their input files.

- **CosMx** requires reading multiple CSVs (metadata, FOV positions, and polygons)

- **Xenium** only requires pointing to the dataset folder, as preprocessing is handled internally

.. code-block:: python

   # CosMx file handling
   zarr_path = "QuarterBrain.zarr"
   slide = "/flatFiles"

   flat_file_dir_slide = path + slide

   metafile = [item for item in os.listdir(flat_file_dir_slide) if 'metadata_file' in item][0]
   metafile_df = pd.read_csv(flat_file_dir_slide+ '/' + metafile)

   fovfile = [item for item in os.listdir(flat_file_dir_slide) if 'fov_positions_file' in item][0]
   fovfile_df = pd.read_csv(flat_file_dir_slide + '/' + fovfile)

   polygon = [item for item in os.listdir(flat_file_dir_slide) if 'polygons' in item][0]
   polygon_df = pd.read_csv(flat_file_dir_slide + '/' + polygon)

   # Xenium file handling
   # No CSV preprocessing needed for Xenium
   path = r"/Volumes/Extreme SSD/Xenium_Skin"
   zarr_path = "Xenium_Skin.zarr"
   os.chdir(path)


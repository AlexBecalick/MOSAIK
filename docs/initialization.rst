Initialization
==============

After file handling (see :doc:`extracting_files`), the datasets must be initialized and written into a Zarr store.  

- **CosMx** requires extra preprocessing (renaming columns and generating cell identifiers)
- **Xenium** can be initialized directly  


CosMx Integration
-----------------

.. code-block:: python

   first_run = input("Is it the first run (0: False, 1: True): ")

   if first_run == '1':
       # Reload CSVs to apply preprocessing
       fovfile_df = pd.read_csv(flat_file_dir_slide + '/' + fovfile)
       fovfile_df = fovfile_df.rename({'FOV': 'fov'}, axis='columns')
       fovfile_df.to_csv(flat_file_dir_slide + '/' + fovfile, index=False)

       polygon_df = pd.read_csv(flat_file_dir_slide + '/' + polygon)
       polygon_df = polygon_df.rename({'cellID': 'cell_ID'}, axis='columns')
       if 'cell' not in polygon_df.keys():
           polygon_df['cell'] = "c_" + polygon_df['fov'].astype(str) + '_' + polygon_df['cell_ID'].astype(str)
       polygon_df.to_csv(flat_file_dir_slide + '/' + polygon, index=False)

       # Load CosMx data (flip_image may be required depending on dataset version)
       sdata = cosmx.cosmx(cosmx_path + slide, flip_image=True)
       sdata.write(zarr_path)

Xenium Integration
------------------

.. code-block:: python

   first_run = input("Is it the first run (0: False, 1: True): ")

   if first_run == '1':
       # Load Xenium data and save as Zarr
       sdata = xenium.xenium(xenium_path)
       sdata.write(zarr_path)

   # Reload data for downstream analysis
   sdata = sd.read_zarr(zarr_path)
   adata = sdata.tables["table"]
   print(adata.obs.keys())

   # Example: Spatial plot
   xy = adata.obsm['spatial']
   plt.scatter(xy[:, 0], xy[:, 1], s=0.0001)

   # Example: Quality Control metrics
   sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
   cprobes = (adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
   cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
   print(f"Negative DNA probe count % : {cprobes}")
   print(f"Negative decoding count % : {cwords}")


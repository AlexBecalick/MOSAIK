Outputs
=======

Using the polygons generated during post-segmentation (see :doc:`post_segmentation`),
we can now create the cells x genes  mapping table and save the results into a new SpatialData Zarr file.

Cell x Gene Matrix
------------------

Assign transcript points to segmented cells and build an AnnData table:

- **gdf_points** – a ``GeoDataFrame`` of transcript points with a ``cell_ID`` column
  linking each point to a segmented cell.
- **vdata** – an ``AnnData`` object representing gene counts per cell.

.. code-block:: python

       def process_points_and_tables(self):
        """Assign transcript points to segmented cells and build AnnData table."""
        points = self.sdata.points[f"{self.fov}_points"]
        points_df = points.compute()

        self.gdf_points = gpd.GeoDataFrame(
            points_df,
            geometry=gpd.points_from_xy(points_df["x"], points_df["y"]),
            crs=self.gdf_polygons.crs,
        )

        joined = gpd.sjoin(
            self.gdf_points,
            self.gdf_polygons[["cell_id", "geometry"]],
            how="left",
            predicate="within",
        )
        self.gdf_points["cell_ID"] = joined["cell_id"].values

        all_genes = self.gdf_points["target"].unique().tolist()
        cell_gene_counts = {}

        for cell in self.gdf_points["cell_ID"].dropna().unique():
            cell_data = self.gdf_points[self.gdf_points["cell_ID"] == cell]
            gene_count = Counter(cell_data["target"])
            ordered_counts = OrderedDict((g, gene_count.get(g, 0)) for g in all_genes)
            cell_gene_counts[str(cell)] = ordered_counts

        df = pd.DataFrame.from_dict(cell_gene_counts, orient="index", columns=all_genes)
        df = df.reindex(columns=all_genes, fill_value=0)

        self.vdata = ad.AnnData(X=df.values)
        self.vdata.obs_names = df.index
        self.vdata.var_names = df.columns
        self.vdata.obs["cells"] = self.vdata.obs_names
        self.vdata.var["genes"] = self.vdata.var_names
        return self.gdf_points, self.vdata

Update SpatialData
------------------

Processed results are saved into a new SpatialData Zarr file. This includes:

- Original images (`images`)  
- Labels (`labels`)  
- Shapes (`shapes`)  
- Original transcript points (`points`)  
- Tables (`vdata`)

.. code-block:: python

       def update_spatialdata(self):
        """Save processed SpatialData into a new Zarr file."""
        self.process_masks_to_shapes()
        self.process_labels()
        self.process_points_and_tables()
        
        new_sdata = sd.SpatialData()
        dst_path = os.path.join(self.zarr_dir, f"{self.fov}.zarr")

        new_sdata.images[f"{self.fov}_image"] = self.sdata.images[f"{self.fov}_image"]
        new_sdata.labels[f"{self.fov}_labels"] = self.labels_model
        new_sdata.shapes[f"{self.fov}_shapes"] = self.shapes_model
        new_sdata.points[f"{self.fov}_points"] = self.sdata.points[f"{self.fov}_points"]
        new_sdata.tables = {'vdata': self.vdata}

        try:
            client = Client.current()
            client.close()
        except ValueError:
            pass

        new_sdata.write(dst_path, overwrite=True)

        for var in (self.masks, self.flows, self.styles):
            try:
                del var
            except Exception:
                pass
        gc.collect()

        return new_sdata


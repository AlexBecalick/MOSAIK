Post-Segmentation
=================

After segmentation, the raw masks are refined and converted into polygons to create labels and enable mapping of gene transcripts in newly segmented cells (see :doc:`outputs` for cells and gene transcripts table). This is done in three main steps:

Masks
-----

1. **Filter masks** using ``filter_cell_by_regionprops`` (removes overly small/cluster cells).

   .. code-block:: python

      def filter_cell_by_regionprops(seg_masks, max_eccentricity=0.95):
          """Filter segmented cell masks by size and eccentricity."""
          labeled = label(seg_masks)
          regions = regionprops(labeled)
          if not regions:
              return np.zeros_like(seg_masks, dtype=np.int32)

          areas = [r.area for r in regions]
          min_area = np.median(areas)
          cleaned_mask = np.zeros_like(seg_masks, dtype=np.int32)
          current_label = 1

          for region in regions:
              if region.area < min_area:
                  continue
              if region.eccentricity > max_eccentricity:
                  continue
              cleaned_mask[labeled == region.label] = current_label
              current_label += 1

          return cleaned_mask

Polygons
--------

2. **Convert masks into polygons** using ``masks_to_polygons`` and upscale back
   to 4256×4256. Then wrap the polygons into a ``GeoDataFrame`` and a ``ShapesModel``.

   .. code-block:: python

      @staticmethod
      def masks_to_polygons(seg_masks):
          """Convert segmentation masks into scaled polygons."""
          polygons = []
          for region in regionprops(seg_masks):
              mask = (seg_masks == region.label).astype(int)
              contours = find_contours(np.array(mask))
              if contours:
                  contour = max(contours, key=lambda x: len(x))
                  if not ((contour[0] == contour[-1]).all()):
                      contour = np.vstack([contour, contour[0]])
                  poly = Polygon(contour[:, [1, 0]]).buffer(0)
                  if poly.is_valid and not poly.is_empty:
                      polygons.append(poly)

          h, w = seg_masks.shape
          if h < 4256 and w < 4256:
              scale_factor = 4256 / h
              polygons = [
                  scale(poly, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))
                  for poly in polygons
              ]
          return polygons

   .. code-block:: python

      def process_masks_to_shapes(self, max_eccentricity: float = 0.95):
          """Filter masks, convert to polygons, and create shapes model."""
          self.masks = self.filter_cell_by_regionprops(self.masks, max_eccentricity)
          polygons = self.masks_to_polygons(self.masks)

          if not polygons:
              raise RuntimeError("No polygons extracted from mask — check segmentation.")

          self.gdf_polygons = gpd.GeoDataFrame(
              {"cell_id": [f"cell_{i+1}" for i in range(len(polygons))], "geometry": polygons}
          )
          shapes_df = gpd.GeoDataFrame(
              {"geometry": polygons, "region": [f"cell_{i+1}" for i in range(len(polygons))]}
          )
          shapes_df.set_geometry("geometry", inplace=True)
          self.shapes_model = ShapesModel.parse(shapes_df)
          return self.masks, self.gdf_polygons, self.shapes_model

Labels
------

3. **Create labels** by resizing masks to the full original resolution and
   converting to a ``Labels2DModel`` for downstream analysis.

   .. code-block:: python

      def process_labels(self, target_shape: tuple = (4256, 4256)):
          """Upscale masks and create Labels2DModel."""
          self.masks = resize(
              self.masks,
              target_shape,
              order=0,
              preserve_range=True,
              anti_aliasing=False,
          ).astype(np.int32)

          self.labels_model = Labels2DModel.parse(data=np.squeeze(self.masks), dims=("y", "x"))
          self.labels_model.name = f"{self.fov}_labels"
          return self.masks, self.labels_model


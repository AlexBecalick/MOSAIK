Post-Segmentation
=================

After segmentation, the raw masks are refined and converted into polygons to create labels and enable mapping of gene transcripts in newly segmented cells (see :doc:`outputs` for cells and gene transcripts table). This is done in three main steps:

Masks
-----

1. **Filter masks** using ``filter_cell_by_regionprops`` (removes artifacts and irregular cells based on multiple morphological criteria).

   .. code-block:: python

      def filter_cell_by_regionprops(
          seg_masks, 
          max_eccentricity=None,
          min_area='median_div2',
          min_absolute_area=50,
          max_area=None,
          min_solidity=None,
          max_solidity=None,
          min_extent=None,
          max_extent=None,
          min_compactness=None,
          max_convexity_deficit=0.20,
          max_perimeter_area_ratio=None,
          area_std_filter=None,
          save_debug=False,
          verbose=True):
          """Filter segmented cell masks by multiple morphological criteria.
          
          Analyzes region properties of segmented masks and filters cells based on
          shape metrics including area, eccentricity, solidity, extent, compactness,
          and convexity deficit. Provides detailed statistics and rejection reasons 
          when verbose mode is enabled.
          
          Args:
              seg_masks (np.ndarray): Input segmentation masks with integer labels.
              max_eccentricity (float, optional): Maximum eccentricity (elongation)
                  threshold. Range: 0-1, where 0=circle, 1=line. Typical: 0.95.
              min_area (str or float, optional): Minimum area filter. Options:
                  - 'median': Use median(areas) as threshold
                  - 'median_div2': Use median(areas)/2 as threshold (default)
                  - 0-100: Use percentile of area distribution
                  - Absolute value: Use as pixel count threshold
              min_absolute_area (int, optional): Hard minimum area to remove 
                  tiny artifacts regardless of distribution. Defaults to 50.
              max_area (float, optional): Maximum area threshold in pixels.
              min_solidity (float, optional): Minimum solidity (area/convex_area).
                  Range: 0-1. Filters highly concave cells. Typical: 0.7-0.9.
              max_solidity (float, optional): Maximum solidity. Filters perfectly
                  convex cells (unrealistic). Typical: 0.98-0.995.
              min_extent (float, optional): Minimum extent (area/bounding_box_area).
                  Range: 0-1. Filters sparse cells. Typical: 0.3-0.5.
              max_extent (float, optional): Maximum extent. Filters perfectly
                  rectangular cells. Typical: 0.9-0.95.
              min_compactness (float, optional): Minimum compactness (Polsby-Popper).
                  Compactness = 4π * area / perimeter². Range: 0-1, where 1=perfect circle.
                  Filters irregular shapes. Typical: 0.5-0.7 for cells.
              max_convexity_deficit (float, optional): Maximum convexity deficit.
                  Deficit = 1 - solidity = (convex_area - area) / convex_area.
                  Range: 0-1, where 0=perfectly convex. Filters concave cells like half-moons.
                  Typical: 0.10-0.20. Defaults to 0.20.
              max_perimeter_area_ratio (float, optional): Maximum perimeter/sqrt(area) ratio.
                  Normalized perimeter measure. Filters cells with excessive perimeter.
                  Typical cells: 3.5-5.0, irregular cells: >6.0.
              area_std_filter (float, optional): Remove area outliers beyond N
                  standard deviations from mean. E.g., 3.0 removes outliers beyond 3σ.
              save_debug (bool, optional): Whether to save debug image 
                  "cleaned_mask.png". Defaults to False.
              verbose (bool, optional): Whether to print filtering statistics
                  and rejection breakdown. Defaults to True.
          
          Returns:
              np.ndarray: Cleaned and relabeled mask array with consecutive 
                  integer labels starting from 1. Zero indicates background.
          """

**Shape Metrics Explained:**

- **Solidity** = area / convex_area
  
  - ~1.0: Very smooth, convex (circle, ellipse) - may be unrealistic
  - 0.7-0.9: Slightly irregular - typical realistic cells
  - <0.7: Highly concave/irregular - likely artifacts

- **Extent** = area / bounding_box_area
  
  - ~1.0: Fills bounding box (square, rectangle)
  - 0.5-0.8: Typical cell shapes
  - <0.3: Very sparse/thin - likely artifacts

- **Eccentricity** = elongation measure
  
  - 0: Perfect circle
  - 0.95: Very elongated
  - >0.95: Extremely elongated - may be artifacts

- **Compactness** = 4π * area / perimeter² (Polsby-Popper)
  
  - 1.0: Perfect circle
  - 0.6-0.8: Typical cells
  - <0.5: Irregular/elongated shapes (including half-moons)

- **Convexity Deficit** = 1 - solidity
  
  - 0: Perfectly convex
  - 0.05-0.15: Slight irregularity (normal cells)
  - >0.20: Significant concavity (half-moon cells, artifacts)

- **Perimeter/Area Ratio** = perimeter / sqrt(area)
  
  - 3.5-5.0: Typical cells
  - >6.0: Irregular perimeter (fragmented, concave cells)

**Note:** The new metrics (compactness, convexity deficit) are especially effective for removing half-moon/crescent-shaped cells.

Polygons
--------

2. **Convert masks into polygons** using ``masks_to_polygons`` and upscale back
   to original resolution. Then wrap the polygons into a ``GeoDataFrame`` and a ``ShapesModel``.

   .. code-block:: python

      @staticmethod
      def masks_to_polygons(seg_masks, factor_rescale):
          """Convert segmentation masks into scaled polygons.
          
          Extracts contours from each labeled region in the segmentation mask,
          converts them to Shapely polygons, and rescales coordinates to match
          the original image resolution.
          
          Args:
              seg_masks (np.ndarray): Labeled segmentation masks where each 
                  unique integer represents a different cell.
              factor_rescale (int): Factor by which to upscale polygon coordinates.
                  E.g., if masks are 4x smaller than original, use factor_rescale=4.
          
          Returns:
              list[shapely.geometry.Polygon]: List of upscaled polygon geometries,
                  one per valid segmented region. Invalid or empty polygons are
                  excluded.
          
          Note:
              - Finds the longest contour for each region
              - Ensures contours are closed (first point = last point)
              - Buffers polygons by 0 to fix any self-intersections
              - If factor_rescale=0, no scaling is applied
          """
          upscale_polygons = []
      
          for region in regionprops(seg_masks):
              mask = (seg_masks == region.label).astype(int)
              contours = find_contours(np.array(mask))
              if contours:
                  contour = max(contours, key=lambda x: len(x))
      
                  # ensure the contour is closed (first point = last point)
                  if not ((contour[0] == contour[-1]).all()):
                      contour = np.vstack([contour, contour[0]])
                  poly = Polygon(contour[:, [1, 0]]).buffer(0)
      
                  # only keep polygons that are valid and not empty
                  if poly.is_valid and not poly.is_empty:
                      upscale_polygons.append(poly)
      
          if factor_rescale != 0:
              upscale_polygons = [
                  scale(poly, xfact=factor_rescale, yfact=factor_rescale, origin=(0, 0))
                  for poly in upscale_polygons
              ]
      
          return upscale_polygons

   .. code-block:: python

      def process_masks_to_shapes(self):
          """Filter masks, convert to polygons, and create shapes model.
          
          Applies morphological filtering to segmentation masks, converts filtered
          masks to polygon geometries, mirrors coordinates, and creates a 
          SpatialData-compatible shapes model with cell IDs.
          
          Returns:
              tuple: Contains three elements:
                  - masks (np.ndarray): Filtered and relabeled segmentation masks.
                  - gdf_polygons (gpd.GeoDataFrame): GeoDataFrame with cell_id and
                    geometry columns.
                  - shapes_model (gpd.GeoDataFrame): SpatialData shapes model with
                    region names and geometries.
          
          Raises:
              RuntimeError: If no valid polygons are extracted from masks after
                  filtering.
          
          Note:
              - Filters masks using filter_cell_by_regionprops()
              - Converts to polygons and upscales by factor_rescale
              - Mirrors y-coordinates and shifts to ensure all y >= 0
              - Sets self.masks, self.gdf_polygons, and self.shapes_model
          """
          # Filter masks
          self.masks = self.filter_cell_by_regionprops(self.masks)

          # Convert masks → polygons
          polygons = self.masks_to_polygons(self.masks, self.factor_rescale)

          if not polygons:
              raise RuntimeError(
                  "No polygons extracted from mask — check your segmentation.")

          # Create shapes model
          shapes_df = gpd.GeoDataFrame({
              "geometry": polygons,
              "region": [f"cell_{i+1}" for i in range(len(polygons))]
          })
          shapes_df.set_geometry("geometry", inplace=True)
          self.shapes_model = ShapesModel.parse(shapes_df)
          
          # Apply coordinate transformations
          self.shapes_model["geometry"] = self.shapes_model["geometry"].apply(
              self.mirror_y0)
          
          miny = self.shapes_model.total_bounds[1]
          if miny < 0:
              self.shapes_model["geometry"] = self.shapes_model["geometry"].apply(
                  lambda g: translate(g, xoff=0, yoff=-miny))
                      
          self.gdf_polygons = self.shapes_model[["region", "geometry"]].rename(
              columns={"region": "cell_id"})

          return self.masks, self.gdf_polygons, self.shapes_model

Labels
------

3. **Create labels** by resizing masks to the full original resolution and
   converting to a ``Labels2DModel`` for downstream analysis.

   .. code-block:: python

      def process_labels(self):
          """Upscale masks and create Labels2DModel.
          
          Upscales the filtered segmentation masks to the original image resolution
          and creates a SpatialData-compatible Labels2DModel for visualization and
          further analysis.
          
          Returns:
              tuple: Contains two elements:
                  - masks (np.ndarray): Original filtered masks (not upscaled).
                  - labels_model (sd.models.Labels2DModel): Upscaled labels model
                    with name set to self.label_name.
          
          Note:
              - Uses nearest-neighbor interpolation (order=0) to preserve integer labels
              - Upscales by factor_rescale to match original image dimensions
              - Sets self.labels_model attribute with specified label name
          """
          rescale_size = tuple(map(lambda x: x * self.factor_rescale, self.masks.shape))
          masks_upscale = resize(
              self.masks,
              rescale_size,
              order=0,
              preserve_range=True,
              anti_aliasing=False,
          ).astype(np.int32)

          self.labels_model = Labels2DModel.parse(
              data=np.squeeze(masks_upscale),
              dims=("y", "x"),
          )
          self.labels_model.name = self.label_name

          return self.masks, self.labels_model

Coordinate Transformation
--------------------------

The pipeline includes a coordinate mirroring step to align the segmentation coordinates with the spatial data coordinate system:

.. code-block:: python

   @staticmethod
   def mirror_y0(geom):
       """Mirror a polygon geometry across the x-axis (y=0).
       
       Flips polygon coordinates vertically by negating y-values. Used to
       correct coordinate system orientation differences between image and
       spatial coordinate systems.
       
       Args:
           geom (shapely.geometry.Polygon or other): Input geometry to mirror.
       
       Returns:
           shapely.geometry: Mirrored geometry of the same type as input.
               If input is not a Polygon, returns input unchanged.
       
       Example:
           Point (x, y) becomes (x, -y)
       """
       return type(geom)([
           (x, -y) for x, y in geom.exterior.coords
       ]) if geom.geom_type == "Polygon" else geom

Advanced Filtering Examples
----------------------------

Basic Filtering
~~~~~~~~~~~~~~~

.. code-block:: python

   # Use default settings (median_div2 area, convexity deficit < 0.20)
   cleaned_masks = filter_cell_by_regionprops(masks)

Strict Filtering (Remove More Artifacts)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Remove elongated cells, irregular shapes, and half-moons
   cleaned_masks = filter_cell_by_regionprops(
       masks,
       max_eccentricity=0.90,           # Remove very elongated cells
       min_compactness=0.55,             # Remove irregular shapes
       max_convexity_deficit=0.15,      # Remove concave/half-moon cells
       max_perimeter_area_ratio=5.0,    # Remove cells with irregular perimeter
       verbose=True
   )

Permissive Filtering (Keep More Cells)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Keep more cells, only remove extreme artifacts
   cleaned_masks = filter_cell_by_regionprops(
       masks,
       min_area='median_div2',
       min_absolute_area=30,             # Lower minimum size
       max_convexity_deficit=0.30,      # Allow more concavity
       verbose=True
   )

Custom Area Filtering
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Use percentile-based area filtering
   cleaned_masks = filter_cell_by_regionprops(
       masks,
       min_area=10,                      # Keep cells above 10th percentile
       max_area=2000,                    # Remove very large cells
       area_std_filter=2.5,              # Remove area outliers beyond 2.5σ
       verbose=True
   )

Debugging
~~~~~~~~~

.. code-block:: python

   # Save debug image and see detailed statistics
   cleaned_masks = filter_cell_by_regionprops(
       masks,
       save_debug=True,
       verbose=True
   )
   # Creates "cleaned_mask.png" for visual inspection

Understanding Filter Statistics
--------------------------------

When ``verbose=True``, the function prints detailed statistics:

**Shape Statistics:**

- Area distribution (range, mean, median, std)
- Eccentricity distribution
- Compactness distribution (key for identifying irregular shapes)
- Solidity and convexity deficit (key for identifying concave cells)
- Perimeter/area ratio (key for identifying fragmented cells)
- Extent distribution (bounding box filling)

**Filtering Results:**

- Number and percentage of regions that passed
- Statistics for valid regions only
- Rejection breakdown by category (e.g., "High Convexity Deficit: 45 (12.3%)")

This helps you understand:

1. What types of cells are being filtered out
2. Whether your thresholds are too strict or too lenient
3. The shape characteristics of your dataset

Performance Notes
-----------------

- **Vectorized operations**: The filter uses vectorized operations and lookup tables for efficient relabeling
- **Memory efficient**: Processes regions without creating intermediate copies
- **Optimized for large datasets**: Can handle thousands of cells efficiently
- **Debug mode**: Use ``save_debug=True`` only when needed, as it creates additional file I/O


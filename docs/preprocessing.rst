Preprocessing
=============

The preprocessing step prepares a raw SpatialData image for segmentation by performing the following operations:

1. Load morphology image from a SpatialData Zarr file.

   .. code-block:: python

      zarr_path = os.path.join(self.zarr_dir, self.zarr_name)
      self.sdata = sd.read_zarr(zarr_path)

      img_xr = self.sdata.images[f"{self.fov}_image"]
      img = img_xr.transpose("y", "x", "c").data.compute()

2. Extract channels: ``Membrane (Y)``, ``DNA (U)``.

   .. code-block:: python

      channel_indices = [channel_names.index(ch) for ch in channels_to_use]
      img = img[..., channel_indices]

3. Normalize contrast (2–98% percentile stretch).

   .. code-block:: python

      p2, p98 = np.percentile(img, (2, 98))
      img_stretched = np.clip((img - p2) / (p98 - p2), 0, 1).astype(np.float32)

4. Convert to 8-bit RGB (required for Cellpose).

   .. code-block:: python

      img_8bit = (img_stretched * 255).astype(np.uint8)

      if img_8bit.ndim == 3 and img_8bit.shape[-1] == 2:
          zero_channel = np.zeros_like(img_8bit[..., :1])
          img_8bit = np.concatenate([img_8bit, zero_channel], axis=-1)
      elif img_8bit.ndim == 2:
          img_8bit = np.stack([img_8bit] * 3, axis=-1)

5. Downsize to ``768×768`` (ideal for Cellpose ViT cyto3 on low RAM).

   .. code-block:: python

      seg_pil = Image.fromarray(img_8bit, mode="RGB")
      seg_pil.thumbnail(thumbnail_size=(768, 768), Image.Resampling.LANCZOS)

6. Save a segmentation-ready ``.png`` file.

   .. code-block:: python

       if output_path is not None:
           seg_pil.save(output_path)

       self.image = seg_pil


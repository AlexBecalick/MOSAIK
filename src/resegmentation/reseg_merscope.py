#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MERSCOPE-specific resegmentation pipeline (Cellpose -> ProSeg)."""

from __future__ import annotations

import gc
import logging
import os
import re
import shutil
import sys
from pathlib import Path

import anndata as ad
import dask.dataframe as dd
import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata as sd
from PIL import Image
from spatialdata.models import PointsModel, TableModel
from tqdm import tqdm

from proseg_wrapper import (
    create_integrated_spatialdata,
    create_visualizations,
    export_data,
    fix_anndata_table,
    fix_zarr_metadata,
    generate_comparison_statistics,
    load_proseg_components,
    run_proseg_refinement,
    validate_and_parse_components,
)
from reseg import Resegmentation_xenium


def _get_reseg_params() -> dict:
    """Load shared resegmentation params from src/params.yaml."""
    src_dir = Path(__file__).resolve().parent.parent
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))
    from params_loader import get_params

    return get_params("resegmentation")


class Resegmentation_merscope(Resegmentation_xenium):
    """MERSCOPE adaptation of Xenium-style resegmentation."""

    def __init__(
        self,
        zarr_dir: str,
        zarr_name: str,
        output: str,
        factor_rescale: int,
        image_name: str,
        label_name: str,
        shape_name: str,
        point_name: str,
        z_start: int,
        z_end: int,
        image_key_prefix: str | None = None,
        transcript_gene_col: str = "gene",
        transcript_x_col: str = "x",
        transcript_y_col: str = "y",
        transcript_z_col: str = "z",
        transcript_cell_id_col: str | None = None,
        min_transcript_qv: float | None = None,
    ):
        super().__init__(
            zarr_dir=zarr_dir,
            zarr_name=zarr_name,
            output=output,
            factor_rescale=factor_rescale,
            image_name=image_name,
            label_name=label_name,
            shape_name=shape_name,
            point_name=point_name,
        )
        self.z_start = int(z_start)
        self.z_end = int(z_end)
        self.image_key_prefix = image_key_prefix
        self.transcript_gene_col = transcript_gene_col
        self.transcript_x_col = transcript_x_col
        self.transcript_y_col = transcript_y_col
        self.transcript_z_col = transcript_z_col
        self.transcript_cell_id_col = transcript_cell_id_col
        self.min_transcript_qv = min_transcript_qv
        self.selected_plane_keys: list[str] = []

    @staticmethod
    def _extract_image_array(image_obj):
        """Extract (H, W, C) numpy array + channel names from a SpatialData image."""
        img_xr = image_obj
        if hasattr(image_obj, "__contains__") and "scale0" in image_obj:
            img_xr = image_obj["scale0"].ds["image"]
        elif hasattr(image_obj, "__getitem__"):
            try:
                img_xr = image_obj["scale0"].ds["image"]
            except Exception:
                img_xr = image_obj

        if hasattr(img_xr, "transpose") and hasattr(img_xr, "dims"):
            dims = list(img_xr.dims)
            if all(d in dims for d in ("y", "x", "c")):
                arr = img_xr.transpose("y", "x", "c").data.compute()
                channels = [str(c) for c in img_xr.coords["c"].values] if "c" in img_xr.coords else []
                return arr, channels

            if all(d in dims for d in ("c", "y", "x")):
                arr = np.moveaxis(img_xr.data.compute(), 0, -1)
                channels = [str(c) for c in img_xr.coords["c"].values] if "c" in img_xr.coords else []
                return arr, channels

        arr = img_xr.data.compute() if hasattr(img_xr, "data") else np.asarray(img_xr)
        if arr.ndim == 2:
            arr = arr[..., np.newaxis]
        elif arr.ndim == 3 and arr.shape[0] <= 8 and arr.shape[-1] > 8:
            arr = np.moveaxis(arr, 0, -1)

        channels = [f"c{i}" for i in range(arr.shape[-1])]
        return arr, channels

    def _resolve_merscope_plane_keys(self) -> list[str]:
        """Resolve image keys for requested inclusive z-range."""
        if self.sdata is None:
            raise RuntimeError("SpatialData is not loaded.")

        pattern = re.compile(r"^(?P<prefix>.+)_z(?P<z>\d+)$")
        keys_by_z: dict[int, str] = {}
        for key in self.sdata.images.keys():
            match = pattern.match(str(key))
            if not match:
                continue
            if self.image_key_prefix and not str(key).startswith(self.image_key_prefix):
                continue
            z_val = int(match.group("z"))
            keys_by_z[z_val] = str(key)

        if not keys_by_z:
            available = list(self.sdata.images.keys())
            raise ValueError(
                "No MERSCOPE z-plane image keys matched pattern '<prefix>_zN'. "
                f"Available images: {available}"
            )

        z0, z1 = sorted((self.z_start, self.z_end))
        missing = [z for z in range(z0, z1 + 1) if z not in keys_by_z]
        if missing:
            raise ValueError(
                f"Requested z-range [{z0}, {z1}] is missing planes: {missing}. "
                f"Available z planes: {sorted(keys_by_z)}"
            )

        return [keys_by_z[z] for z in range(z0, z1 + 1)]

    def preprocess_image(
        self,
        channel_names: list[str] | None = None,
        channels_to_use: list[str] | None = None,
    ):
        """Create max-projection across z-range, then preprocess for Cellpose."""
        zarr_path = os.path.join(self.zarr_dir, self.zarr_name)
        self.sdata = sd.read_zarr(zarr_path)

        self.selected_plane_keys = self._resolve_merscope_plane_keys()
        projection_stack = []
        channels_from_data = None

        for key in self.selected_plane_keys:
            img_array, chs = self._extract_image_array(self.sdata.images[key])
            projection_stack.append(img_array)
            if channels_from_data is None:
                channels_from_data = chs

        img = np.max(np.stack(projection_stack, axis=0), axis=0)
        channels_from_data = channels_from_data or []

        if channels_to_use:
            ch_idx = []
            for ch in channels_to_use:
                if ch in channels_from_data:
                    ch_idx.append(channels_from_data.index(ch))
            if not ch_idx:
                raise ValueError(
                    f"None of channels_to_use={channels_to_use} found in {channels_from_data}"
                )
            img = img[..., ch_idx]

        # Keep Cellpose input as 1-3 channels.
        if img.shape[-1] > 3:
            img = img[..., :3]

        p2, p98 = np.percentile(img, (2, 98))
        img_stretched = np.clip((img - p2) / (p98 - p2 + 1e-8), 0, 1).astype(np.float32)
        img_8bit = (img_stretched * 255).astype(np.uint8)

        if img_8bit.ndim == 2:
            img_8bit = np.stack([img_8bit] * 3, axis=-1)
        elif img_8bit.shape[-1] == 1:
            img_8bit = np.repeat(img_8bit, 3, axis=-1)
        elif img_8bit.shape[-1] == 2:
            img_8bit = np.concatenate([img_8bit, np.zeros_like(img_8bit[..., :1])], axis=-1)

        seg_pil = Image.fromarray(img_8bit, mode="RGB")
        y_downscale = int(img_8bit.shape[0] / self.factor_rescale)
        x_downscale = int(img_8bit.shape[1] / self.factor_rescale)
        seg_pil.thumbnail((x_downscale, y_downscale), Image.Resampling.LANCZOS)
        if self.output is not None:
            seg_pil.save(self.output)

        # Keep first plane key for update_spatialdata image copy.
        self.image_name = self.selected_plane_keys[0]
        self.image = seg_pil
        self.clear_gpu_memory()
        return seg_pil

    def process_points_and_tables(self):
        """MERSCOPE transcript assignment and table construction."""
        point_frames = []
        for key in self.sdata.points.keys():
            points = self.sdata.points[key]
            frame = points.compute() if hasattr(points, "compute") else points.copy()
            frame = pd.DataFrame(frame)
            if not frame.empty:
                point_frames.append(frame)

        if not point_frames:
            raise RuntimeError("No transcript point layers found in SpatialData.")

        points_df = pd.concat(point_frames, ignore_index=True)
        columns = set(points_df.columns)

        def resolve_col(preferred: str | None, fallbacks: list[str], required: bool = True):
            if preferred and preferred in columns:
                return preferred
            for col in fallbacks:
                if col in columns:
                    return col
            if required:
                raise KeyError(f"Could not resolve required column. Tried: {preferred}, {fallbacks}")
            return None

        x_col = resolve_col(self.transcript_x_col, ["x", "global_x", "x_location"])
        y_col = resolve_col(self.transcript_y_col, ["y", "global_y", "y_location"])
        z_col = resolve_col(self.transcript_z_col, ["z", "global_z", "z_location"], required=False)
        gene_col = resolve_col(self.transcript_gene_col, ["feature_name", "gene", "target"])
        qv_col = "qv" if "qv" in columns else None

        if self.min_transcript_qv is not None and qv_col is not None:
            points_df = points_df[points_df[qv_col] >= self.min_transcript_qv]

        normalized = pd.DataFrame(
            {
                "x": points_df[x_col],
                "y": points_df[y_col],
                "z": points_df[z_col] if z_col is not None else 0,
                "feature_name": points_df[gene_col].astype(str),
            }
        )
        normalized = normalized.dropna(subset=["x", "y", "feature_name"])

        self.gdf_points = gpd.GeoDataFrame(
            normalized.copy(),
            geometry=gpd.points_from_xy(normalized["x"], normalized["y"]),
            crs=self.gdf_polygons.crs,
        )
        self.gdf_points["x_updated"] = self.gdf_points["x"].astype(float)
        self.gdf_points["y_updated"] = self.gdf_points["y"].astype(float)

        joined = gpd.sjoin(
            self.gdf_points,
            self.gdf_polygons[["cell_id", "geometry"]],
            how="left",
            predicate="within",
        )
        self.gdf_points["cell_id"] = joined["cell_id"].values

        assigned = self.gdf_points.dropna(subset=["cell_id"]).copy()
        assigned["cell_id"] = assigned["cell_id"].astype(str)
        all_genes = sorted(assigned["feature_name"].unique().tolist()) if not assigned.empty else []

        if assigned.empty:
            cell_ids = self.gdf_polygons["cell_id"].astype(str).tolist()
            df = pd.DataFrame(index=cell_ids)
        else:
            cell_gene_matrix = pd.crosstab(assigned["cell_id"], assigned["feature_name"])
            cell_gene_matrix = cell_gene_matrix.reindex(columns=all_genes, fill_value=0)
            cell_gene_counts = {
                str(cell): row.to_dict() for cell, row in tqdm(cell_gene_matrix.iterrows(), total=len(cell_gene_matrix))
            }
            df = pd.DataFrame.from_dict(cell_gene_counts, orient="index", columns=all_genes)
            df = df.reindex(columns=all_genes, fill_value=0)

        self.vdata = ad.AnnData(X=df.values if df.shape[1] > 0 else np.zeros((len(df.index), 0)))
        self.vdata.obs_names = pd.Index(df.index.astype(str), dtype=str)
        self.vdata.var_names = pd.Index(df.columns.astype(str), dtype=str)
        self.vdata.obs["cell_id"] = self.vdata.obs_names
        self.vdata.obs["region"] = pd.Categorical([self.shape_name] * self.vdata.n_obs)
        self.vdata.var["genes"] = self.vdata.var_names

        points_out = self.gdf_points[["x", "y", "z", "feature_name", "cell_id", "x_updated", "y_updated"]].copy()
        points_out["cell_id"] = points_out["cell_id"].where(points_out["cell_id"].notna(), "0").astype(str)
        n_partitions = max(1, min(16, len(points_out) // 50000 + 1))
        points_dd = dd.from_pandas(points_out, npartitions=n_partitions)
        self.points_model = PointsModel.parse(
            points_dd,
            coordinates={"x": "x", "y": "y", "z": "z"},
            feature_key="feature_name",
            instance_key="cell_id",
        )

        self.clear_gpu_memory()
        return self.points_model, self.vdata

    def update_spatialdata(self, proseg_refinement: bool = True):
        """Run full MERSCOPE pipeline and write updated Zarr."""
        self.process_masks_to_shapes()
        self.process_labels()
        self.process_points_and_tables()

        self.new_sdata = sd.SpatialData()
        dst_path = os.path.join(self.output_folder, f"updated_{self.zarr_name}")

        # Keep selected source image from the original object.
        self.new_sdata.images[self.image_name] = self.sdata.images[self.image_name]
        self.new_sdata.labels[self.label_name] = self.labels_model
        self.new_sdata.shapes[self.shape_name] = self.shapes_model
        self.new_sdata.points[self.point_name] = self.points_model

        table = TableModel.parse(
            self.vdata.copy(),
            region_key="region",
            region=self.shape_name,
            instance_key="cell_id",
        )
        self.new_sdata.tables["table"] = table
        self.new_sdata.write(dst_path, overwrite=True)

        if proseg_refinement:
            self.run_proseg()

        for var in (self.masks, self.flows, self.styles, self.gdf_points, self.gdf_polygons):
            try:
                del var
            except Exception:
                pass
        gc.collect()
        self.clear_gpu_memory()
        return self.new_sdata

    def run_proseg(
        self,
        samples: int = 1000,
        voxel_size: float = 0.5,
        voxel_layers: int = 2,
        nuclear_reassignment_prob: float = 0.2,
        diffusion_probability: float = 0.2,
        num_threads: int = 12,
    ):
        """Run ProSeg refinement for MERSCOPE output."""
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
        self.log = logging.getLogger(__name__)
        self.log.info("Starting ProSeg refinement pipeline for MERSCOPE data")

        if self.new_sdata is None:
            raise RuntimeError("update_spatialdata() must run before run_proseg().")

        points = self.new_sdata.points[self.point_name]
        transcripts_df = points.compute() if hasattr(points, "compute") else pd.DataFrame(points).copy()

        # Normalize required ProSeg input columns.
        if "feature_name" not in transcripts_df.columns and "gene" in transcripts_df.columns:
            transcripts_df["feature_name"] = transcripts_df["gene"].astype(str)
        if "x_updated" not in transcripts_df.columns and "x" in transcripts_df.columns:
            transcripts_df["x_updated"] = transcripts_df["x"]
        if "y_updated" not in transcripts_df.columns and "y" in transcripts_df.columns:
            transcripts_df["y_updated"] = transcripts_df["y"]
        if "cell_id" not in transcripts_df.columns:
            transcripts_df["cell_id"] = "0"
        transcripts_df["cell_id"] = transcripts_df["cell_id"].where(transcripts_df["cell_id"].notna(), "0").astype(str)

        output_path = Path("proseg_output.zarr")
        output_integrated_path = Path("integrated_proseg_output.zarr")
        reseg_params = _get_reseg_params()

        output_path = run_proseg_refinement(
            transcripts_df=transcripts_df,
            output_path=str(output_path),
            proseg_binary=reseg_params["proseg_binary"],
            x_col="x_updated",
            y_col="y_updated",
            z_col="z",
            gene_col="feature_name",
            cell_id_col="cell_id",
            samples=samples,
            voxel_size=voxel_size,
            voxel_layers=voxel_layers,
            nuclear_reassignment_prob=nuclear_reassignment_prob,
            diffusion_probability=diffusion_probability,
            num_threads=num_threads,
            overwrite=True,
            logger=self.log,
        )

        fix_zarr_metadata(output_path, self.log)
        fix_anndata_table(output_path, self.log)

        refined_shapes_gdf, adata_proseg, refined_transcripts_df = load_proseg_components(output_path, self.log)
        refined_sdata, adata_sanitized, refined_shapes_parsed, instance_key = validate_and_parse_components(
            output_path, refined_shapes_gdf, adata_proseg, self.log
        )
        integrated_sdata = create_integrated_spatialdata(
            self.new_sdata,
            output_path,
            refined_shapes_parsed,
            adata_sanitized,
            instance_key,
            transcripts_df,
            refined_transcripts_df,
            self.log,
        )

        create_visualizations(integrated_sdata, refined_shapes_gdf, self.log)
        generate_comparison_statistics(self.new_sdata, refined_sdata, integrated_sdata, self.log)
        export_data(refined_shapes_gdf, adata_proseg, refined_transcripts_df, self.log)

        self.integrated_sdata = integrated_sdata
        if output_integrated_path.exists():
            shutil.rmtree(output_integrated_path)
        self.integrated_sdata.write(output_integrated_path)
        self.log.info("ProSeg refinement complete")
        return self.integrated_sdata

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Batch MERSCOPE resegmentation workflow (Cellpose -> ProSeg)."""

import logging
import os
import shutil
import sys
import warnings
from pathlib import Path

import pandas as pd
import spatialdata as sd

warnings.filterwarnings("ignore")

# Load params from src/params.yaml.
src_dir = Path(__file__).resolve().parent.parent
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))
from params_loader import get_params

from reseg_merscope import Resegmentation_merscope


def _safe_move(src_name: str, output_dir: Path):
    """Move file or directory to output_dir, replacing destination if present."""
    if not os.path.exists(src_name):
        print(f"File not found: {src_name}")
        return
    destination = output_dir / Path(src_name).name
    if destination.exists():
        if destination.is_dir():
            shutil.rmtree(destination)
        else:
            destination.unlink()
    shutil.move(src_name, output_dir)


def main():
    logging.basicConfig(level=logging.WARNING)
    params = get_params("resegmentation_merscope")

    base_path = Path(params["base_path"])
    os.chdir(base_path)

    metadata = pd.read_csv(params["metadata_path"])
    zarr_dir = params["zarr_dir"]
    zarr_name_pattern = params["zarr_name_pattern"]
    sample_id_col = params["sample_id_column"]
    tissue_col = params["tissue_filter_column"]
    tissue_value = params["tissue_filter_value"]
    output_dir_pattern = params["output_dir_pattern"]
    output_files_to_move = params["output_files_to_move"]
    final_example_zarr = params.get("final_example_zarr")

    for _, row in metadata.iterrows():
        if tissue_value is not None and row[tissue_col] != tissue_value:
            continue

        sample_name = str(row[sample_id_col])
        output_dir = base_path / output_dir_pattern.format(sample_name=sample_name)
        output_dir.mkdir(exist_ok=True)

        zarr_name = zarr_name_pattern.format(sample_name=sample_name)
        zarr_path = Path(zarr_dir) / zarr_name
        if not zarr_path.exists():
            print(f"Missing input zarr: {zarr_path}")
            continue

        print(f"\n[Sample] {sample_name}")
        core = sd.read_zarr(str(zarr_path))
        print(f"Loaded: {zarr_path}")

        # Optional channel list from params (null -> all channels in projection).
        channels_to_use = params.get("channels_to_use")

        pipe = Resegmentation_merscope(
            zarr_dir=zarr_dir,
            zarr_name=zarr_name,
            output=params["preprocessed_image_output"],
            factor_rescale=params["factor_rescale"],
            image_name="morphology_focus",
            label_name=params["label_name"],
            shape_name=params["shape_name"],
            point_name=params["point_name"],
            z_start=params["z_start"],
            z_end=params["z_end"],
            image_key_prefix=params.get("image_key_prefix"),
            transcript_gene_col=params["transcript_gene_col"],
            transcript_x_col=params["transcript_x_col"],
            transcript_y_col=params["transcript_y_col"],
            transcript_z_col=params["transcript_z_col"],
            transcript_cell_id_col=params.get("transcript_cell_id_col"),
            min_transcript_qv=params.get("min_transcript_qv"),
        )

        pipe.preprocess_image(channels_to_use=channels_to_use)
        pipe.run_cellpose(
            flow_threshold=params["flow_threshold"],
            cellprob_threshold=params["cellprob_threshold"],
            tile_overlap=params["tile_overlap"],
        )

        run_proseg = bool(params.get("run_proseg", True))
        pipe.update_spatialdata(proseg_refinement=run_proseg)

        for filename in output_files_to_move:
            _safe_move(filename, output_dir)

        if run_proseg:
            _safe_move("proseg_output.zarr", output_dir)
            _safe_move("integrated_proseg_output.zarr", output_dir)

        del core
        print(f"Completed: {sample_name}")

    if final_example_zarr and os.path.exists(final_example_zarr):
        _ = sd.read_zarr(final_example_zarr)
        print(f"Loaded final example zarr: {final_example_zarr}")


if __name__ == "__main__":
    main()

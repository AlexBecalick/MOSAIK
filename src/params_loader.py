"""
Load MOSAIK parameters from src/params.yaml.
Use get_params() to obtain the full dict, or get_params(section) for a subsection.
"""
from pathlib import Path

import yaml

_PARAMS_PATH = Path(__file__).resolve().parent / "params.yaml"

_cached_params = None


def get_params(section: str | None = None) -> dict:
    """
    Load and return parameters from src/params.yaml.

    Parameters
    ----------
    section : str, optional
        If provided, return only this top-level section (e.g. "merscope_qc",
        "resegmentation_xenium"). Otherwise return the full params dict.

    Returns
    -------
    dict
        Parameter key-value map for the requested section or full file.
    """
    global _cached_params
    if _cached_params is None:
        if not _PARAMS_PATH.exists():
            raise FileNotFoundError(
                f"Params file not found: {_PARAMS_PATH}. "
                "Copy or create params.yaml in the src folder and edit paths for your setup."
            )
        with open(_PARAMS_PATH) as f:
            _cached_params = yaml.safe_load(f)
    if section is None:
        return _cached_params.copy()
    if section not in _cached_params:
        raise KeyError(f"Section {section!r} not found in params. Available: {list(_cached_params)}")
    return _cached_params[section].copy()


def get_params_path() -> Path:
    """Return the path to the params.yaml file (for reference or overrides)."""
    return _PARAMS_PATH

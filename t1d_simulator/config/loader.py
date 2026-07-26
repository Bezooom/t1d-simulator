"""
Configuration loader for t1d_simulator.
Wraps param_loader.py with convenient accessor helpers.
"""
import os
from t1d_simulator.param_loader import load_parameters

_PARAMS = None
_PARAMS_PATH = None


def _ensure_loaded(filepath=None):
    global _PARAMS, _PARAMS_PATH
    if filepath and filepath != _PARAMS_PATH:
        _PARAMS = None
        _PARAMS_PATH = filepath
    if _PARAMS is None:
        _PARAMS = load_parameters(_PARAMS_PATH)
    return _PARAMS


def load_params(filepath=None):
    """Load and return the full parameters dict."""
    return _ensure_loaded(filepath)


def get_param(key, default=None):
    """
    Dot-notation parameter access.
    get_param("mechanical.e_0_default_kpa") -> 50.0
    get_param("k_m") -> 0.5
    """
    params = _ensure_loaded()
    parts = key.split(".")
    obj = params
    for part in parts:
        if isinstance(obj, dict) and part in obj:
            obj = obj[part]
        else:
            return default
    return obj


def reload_params(filepath=None):
    """Force reload of parameters (clears cache)."""
    global _PARAMS, _PARAMS_PATH
    _PARAMS = None
    _PARAMS_PATH = filepath
    return load_params(filepath)

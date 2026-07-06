"""Tests for the simulation cache (raw/cache.py) — no pyechelle runs."""

from pathlib import Path

import numpy as np

from andes_simulator.raw.cache import SimCache

PROJECT_ROOT = Path(__file__).parent.parent


def make_cache(tmp_path, boost=10.0):
    return SimCache(tmp_path / "cache", PROJECT_ROOT, boost=boost)


def test_signature_stable_and_boost_sensitive(tmp_path):
    cache = make_cache(tmp_path)
    sig1 = cache.signature("Y", {"type": "fabry_perot"}, (1, 38, 75))
    sig2 = cache.signature("Y", {"type": "fabry_perot"}, (1, 38, 75))
    assert cache.path_for(sig1) == cache.path_for(sig2)

    other_boost = make_cache(tmp_path, boost=20.0)
    sig3 = other_boost.signature("Y", {"type": "fabry_perot"}, (1, 38, 75))
    assert other_boost.path_for(sig3) != cache.path_for(sig1)


def test_signature_includes_csv_file_stat(tmp_path):
    cache = make_cache(tmp_path)
    sig = cache.signature("R", {"type": "csv", "filepath": "SED/phoenix.csv"}, (1,))
    assert "source_file_stat" in sig


def test_store_load_roundtrip(tmp_path):
    cache = make_cache(tmp_path)
    sig = cache.signature("Y", {"type": "constant"}, (5, 6))
    rate = np.arange(12, dtype=np.float32).reshape(3, 4)
    cache.store(sig, rate)
    loaded = cache.load(sig)
    assert loaded is not None
    assert loaded.dtype == np.float32
    assert np.array_equal(loaded, rate)


def test_load_missing_returns_none(tmp_path):
    cache = make_cache(tmp_path)
    sig = cache.signature("Y", {"type": "lfc"}, (2,))
    assert cache.load(sig) is None


def test_resolve_fibers_drops_dead_ones(tmp_path):
    cache = make_cache(tmp_path)
    # YJH cal fibers 1,37,38,39,75; fibers 37 and 39 are dead in the model
    assert cache.resolve_fibers("Y", [1, 37, 38, 39, 75]) == (1, 38, 75)
    assert cache.resolve_fibers("R", [33, 34]) == (33, 34)

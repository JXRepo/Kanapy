"""
Exercise the file and plotting behavior of the complete EBSD example.
"""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import networkx as nx
import numpy as np
from orix.quaternion import symmetry


def load_example():
    """
    Load the example without launching the expensive RVE pipeline.

    Returns
    -------
    module
        Example entry point and its focused output helpers.
    """
    path = Path(__file__).parents[1] / "examples/EBSD_to_RVE/graph_to_rve.py"
    spec = importlib.util.spec_from_file_location("graph_to_rve", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_example_parameters_control_defaults_and_cli_overrides(tmp_path, monkeypatch):
    """
    Use the editable example settings unless the command line overrides them.

    This prevents the visible parameter block from becoming disconnected from
    the actual run configuration.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Alternate input and output paths.
    monkeypatch : pytest.MonkeyPatch
        Temporary edits to the example settings.
    """
    example = load_example()
    assert example.parse_args([]).connectivity == 4
    settings = {
        "fname": tmp_path / "measurement.ang", "output_dir": tmp_path / "results",
        "nvox": 12, "box_length": 30, "periodic": True, "seed": 42,
        "min_size": 20., "max_angle": 7.5, "vf_min": 0.05,
        "connectivity": 4, "show_figures": True, "extra_outputs": True,
    }
    for name, value in settings.items():
        monkeypatch.setattr(example, name, value)

    args = example.parse_args([])

    assert args.input == settings["fname"]
    assert args.output_dir == settings["output_dir"]
    assert (args.voxels, args.size, args.periodic, args.seed) == (12, 30, True, 42)
    assert (args.min_size, args.max_angle, args.vf_min) == (20., 7.5, 0.05)
    assert args.show and args.extra_outputs
    overridden = example.parse_args(["--connectivity", "8", "--min-size", "15", "--no-show"])
    assert (overridden.connectivity, overridden.min_size, overridden.show) == (8, 15., False)
    assert overridden.voxels == 12
    assert overridden.max_angle == 7.5


def test_manifest_tracks_only_current_artifacts(tmp_path):
    """
    Exclude stale files while retaining numeric metadata in JSON.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary output directory.
    """
    example = load_example()
    (tmp_path / "stale.png").write_bytes(b"old")
    current = tmp_path / "current.png"
    current.write_bytes(b"new")
    target = tmp_path / "manifest.json"
    example.write_manifest(target, {"dimensions": np.array([30, 30, 30])}, [current])
    manifest = json.loads(target.read_text())
    assert manifest["dimensions"] == [30, 30, 30]
    assert manifest["files"] == ["current.png"]
    assert (tmp_path / "stale.png").read_bytes() == b"old"


def test_statistics_preserve_raw_distributions_and_units(tmp_path):
    """
    Export experimental arrays as numbers with their statistical meaning.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary output directory.
    """
    example = load_example()
    data = {"index": 7, "name": "austenite", "vf": np.float64(0.9), "ngrains": 2}
    for prefix in ["gs", "ar", "om"]:
        data[f"{prefix}_data"] = np.array([1.0, 2.0])
        data[f"{prefix}_param"] = np.array([0.2, 0.0, 1.4])
        data[f"{prefix}_moments"] = [np.float64(1.4), np.float64(0.2)]
    result = example.knpy.EBSDGraphResult(
        graph=nx.path_graph(2), phase_data=data, phase_index=7,
        phase_name="austenite", rgb_image=np.zeros((1, 2, 3)),
        npx=2, sh_x=1, sh_y=2,
    )
    target = tmp_path / "statistics.json"
    example.write_statistics(target, result, {"Phase": {"Number": 0}})
    statistics = json.loads(target.read_text())
    assert statistics["experimental"]["gs_data"] == [1.0, 2.0]
    assert statistics["experimental"]["index"] == 7
    assert statistics["units"] == {"length": "um", "tilt_angle": "rad", "aspect_ratio": "1"}
    assert statistics["rve_descriptor"]["Phase"]["Number"] == 0


def test_mesh_orientation_colors_use_experimental_symmetry():
    """
    Keep cubic IPF colors when converting generated and imported Euler angles.
    """
    example = load_example()
    mesh = SimpleNamespace(grain_ori_dict={3: [0, 0, 0], 8: [0.2, 0.5, 0.8]})
    orientations = example.mesh_orientations(mesh, symmetry.Oh)
    assert orientations.symmetry == symmetry.Oh
    np.testing.assert_allclose(orientations.to_euler(), [[0, 0, 0], [0.2, 0.5, 0.8]], atol=1e-7)
    colors = example.knpy.get_ipf_colors(orientations)
    assert colors.shape == (2, 3)
    assert not np.allclose(colors[0], colors[1])

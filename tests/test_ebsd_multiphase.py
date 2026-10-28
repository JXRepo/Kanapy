"""
Exercise phase boundaries and shared EBSD figures while preserving single phase analysis.
"""

import copy
import json
import pickle
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pytest
from orix.crystal_map import CrystalMap, PhaseList
from orix.quaternion import Orientation, symmetry

from kanapy.graph_workflow import (
    EBSDGraphConfig, EBSDGraphOutputOptions, EBSDGraphResult, write_graph_result_outputs,
)
from kanapy.texture import (
    EBSDmap, _compute_all_edge_misorientations, _refresh_node_misorientation_cache,
    find_largest_neighbor, find_sim_neighbor, get_proper_symmetry_quaternions,
    merge_nodes,
)


def make_mixed_graph():
    """
    Create a misleading cross phase neighbor and a valid same phase candidate.

    The other phase is both larger and closer in orientation, so unconstrained
    neighbor searches select the wrong grain.

    Returns
    -------
    networkx.Graph
        Three grain graph with a shared boundary between different phases.
    """
    graph = nx.Graph(label_map=np.array([[1, 2, 3]]), dx=1., dy=1., symmetry=symmetry.Oh)
    for node, phase, size, angle in [(1, 2, 1, 0.), (2, 7, 100, 0.), (3, 2, 5, 0.2)]:
        graph.add_node(node, phase_id=phase, phase_name=str(phase),
                       pixels=np.array([node - 1]), npix=size,
                       ori_av=np.array([np.cos(angle / 2), 0., 0., np.sin(angle / 2)]))
    graph.add_edges_from([(1, 2), (1, 3)])
    return graph


def test_neighbor_searches_exclude_other_phases():
    """
    Prevent another phase from winning size or orientation based merge selection.
    """
    graph = make_mixed_graph()
    assert find_largest_neighbor(graph, 1) == 3
    assert find_sim_neighbor(graph, 1)[0] == 3
    graph.remove_edge(1, 3)
    with pytest.raises(ValueError):
        find_largest_neighbor(graph, 1)
    with pytest.raises(ValueError):
        find_sim_neighbor(graph, 1)


def test_cross_phase_merge_rejected_before_mutation():
    """
    Leave every graph field unchanged when a direct cross phase merge is requested.
    """
    graph = make_mixed_graph()
    before = copy.deepcopy(graph)
    with pytest.raises(ValueError, match="phase"):
        merge_nodes(graph, 1, 2)
    assert set(graph.nodes) == set(before.nodes)
    assert set(graph.edges) == set(before.edges)
    np.testing.assert_array_equal(graph.graph["label_map"], before.graph["label_map"])
    for node in before:
        for key, value in before.nodes[node].items():
            np.testing.assert_equal(graph.nodes[node][key], value)
    merge_nodes(graph, 1, 3)
    assert graph.nodes[3]["phase_id"] == 2
    assert graph.nodes[3]["phase_name"] == "2"


def test_tagged_grain_does_not_merge_with_unknown_phase():
    """
    Reject mixed known and unknown phase metadata while supporting legacy graphs.
    """
    graph = make_mixed_graph()
    del graph.nodes[3]["phase_id"]
    with pytest.raises(ValueError, match="phase"):
        merge_nodes(graph, 1, 3)
    del graph.nodes[1]["phase_id"]
    merge_nodes(graph, 1, 3)
    assert 1 not in graph


@pytest.mark.parametrize("refresh", [False, True])
def test_misorientation_cache_excludes_cross_phase_edges(refresh):
    """
    Clear stale cross phase angles instead of treating them as merge candidates.

    Parameters
    ----------
    refresh : bool
        Recompute all edges or only the selected node's incident edges.
    """
    graph = make_mixed_graph()
    graph[1][2]["misorientation"] = 0.
    if refresh:
        _refresh_node_misorientation_cache(graph, 1)
    else:
        _compute_all_edge_misorientations(graph, get_proper_symmetry_quaternions(symmetry.Oh))
    assert "misorientation" not in graph[1][2]
    assert graph[1][3]["misorientation"] == pytest.approx(0.2)


@pytest.fixture
def square_map():
    """
    Build two phases with two adjacent grains each at known image positions.

    Returns
    -------
    kanapy.texture.EBSDmap
        Map analyzed without grain statistics or interactive figures.
    """
    rows, columns = np.indices((8, 20))
    phases = np.where(columns.ravel() < 10, 2, 7)
    euler = np.zeros((160, 3))
    euler[:, 0] = np.deg2rad((columns.ravel() // 5) * 20)
    emap = CrystalMap(
        rotations=Orientation.from_euler(euler), phase_id=phases,
        x=columns.ravel().astype(float), y=rows.ravel().astype(float),
        phase_list=PhaseList(names=["A", "B"], point_groups=["m-3m", "6/mmm"], ids=[2, 7]),
    )
    data = SimpleNamespace(emap=emap, sh_x=8, sh_y=20, dx=1., dy=1., npx=160)
    return EBSDmap(ebsd_data=data, graph_only=True)


def make_hex_data():
    """
    Build two phases occupying opposite sides of a translated staggered scan.

    Returns
    -------
    types.SimpleNamespace
        Measurements and original grid metadata, before EBSD analysis.
    """
    xy = []
    for row in range(6):
        for column in range(8 - row % 2):
            xy.append([10. + column + 0.5 * (row % 2), 20. + row * np.sqrt(3.) / 2.])
    xy = np.asarray(xy)
    phases = np.where(xy[:, 0] < 13.75, 2, 7)
    euler = np.zeros((45, 3))
    euler[phases == 7] = [0.4, 0.3, 0.2]
    emap = CrystalMap(
        rotations=Orientation.from_euler(euler), phase_id=phases, x=xy[:, 0], y=xy[:, 1],
        phase_list=PhaseList(names=["A", "B"], point_groups=["m-3m", "6/mmm"], ids=[2, 7]),
        prop={"ci": np.where(phases == 2, 0.8, 0.4)},
    )
    return SimpleNamespace(emap=emap, sh_x=6, sh_y=8, dx=1., dy=np.sqrt(3.) / 2.,
                           npx=48, original_npx=45, hex_grid=True)


def test_nodes_keep_phase_identity_and_regions_do_not_cross_phases(square_map):
    """
    Retain phase labels through cleanup and exclude pixels of the other phase.

    Parameters
    ----------
    square_map : kanapy.texture.EBSDmap
        Two phase map with independently known grain counts.
    """
    for data in square_map.ms_data:
        graph = data["graph"]
        assert len(graph) == 2
        for node in graph.nodes.values():
            assert node["phase_id"] == data["index"]
            assert node["phase_name"] == data["name"]
            np.testing.assert_array_equal(square_map.emap.phase_id[node["pixels"]], data["index"])


def test_multiphase_hex_uses_disjoint_masks_on_one_full_grid():
    """
    Preserve the left and right phases instead of expanding both across the image.
    """
    ebsd = EBSDmap(ebsd_data=make_hex_data(), graph_only=True, gs_min=1)
    first, second = ebsd.ms_data
    assert np.intersect1d(first["idx"], second["idx"]).size == 0
    np.testing.assert_array_equal(np.union1d(first["idx"], second["idx"]), np.arange(48))
    labels = ebsd._hex_phase_grid
    np.testing.assert_array_equal(labels[:, 0], 2)
    np.testing.assert_array_equal(labels[:, -1], 7)
    np.testing.assert_allclose(ebsd._hex_x_grid[0], np.arange(10., 18.))
    for data in ebsd.ms_data:
        assert len(data["graph"]) == 1
        expected_ci = 0.8 if data["index"] == 2 else 0.4
        np.testing.assert_allclose(ebsd._hex_ci_grid.ravel()[data["idx"]], expected_ci)
        angles = [0., 0., 0.] if data["index"] == 2 else [0.4, 0.3, 0.2]
        quaternion = Orientation.from_euler(angles).data[0]
        np.testing.assert_allclose(np.abs(data["ori"].data @ quaternion), 1.)


def test_hex_filter_leaving_one_phase_preserves_existing_resampling():
    """
    Keep the old one phase grid and filling behavior after the other phase is excluded.
    """
    ebsd = EBSDmap(ebsd_data=make_hex_data(), graph_only=True, gs_min=1, vf_min=0.49)
    assert len(ebsd.ms_data) == 1
    assert ebsd.ms_data[0]["index"] == 2
    assert not hasattr(ebsd, "_hex_phase_grid")
    np.testing.assert_array_equal(ebsd.ms_data[0]["idx"], np.arange(48))
    assert (ebsd.dx, ebsd.dy) == (1., np.sqrt(3.) / 2.)


@pytest.mark.parametrize("minor_pixels, expected_phases", [(2, [2]), (3, [2, 7])])
def test_phase_fraction_cutoff_and_equal_orientation_boundaries(minor_pixels, expected_phases):
    """
    Retain exactly three percent and separate touching phases with identical orientations.

    Parameters
    ----------
    minor_pixels : int
        Number of minor phase pixels among one hundred measurements.
    expected_phases : list of int
        Independently expected retained phase IDs at the default cutoff.
    """
    rows, columns = np.indices((10, 10))
    phases = np.full(100, 2)
    phases[np.array([0, 1, 10])[:minor_pixels]] = 7
    emap = CrystalMap(
        rotations=Orientation.from_euler(np.zeros((100, 3))), phase_id=phases,
        x=columns.ravel().astype(float), y=rows.ravel().astype(float),
        phase_list=PhaseList(names=["A", "B"], point_groups=["m-3m"] * 2, ids=[2, 7]),
    )
    loaded = SimpleNamespace(emap=emap, sh_x=10, sh_y=10, dx=1., dy=1., npx=100)
    ebsd = EBSDmap(ebsd_data=loaded, graph_only=True, gs_min=1)
    assert [data["index"] for data in ebsd.ms_data] == expected_phases
    for data in ebsd.ms_data:
        assert len(data["graph"]) == 1
        for node in data["graph"].nodes.values():
            np.testing.assert_array_equal(phases[node["pixels"]], data["index"])


def test_combined_overlay_preserves_background_and_phase_colors(square_map):
    """
    Plot both phase graphs on one IPF image, including graph only nodes without centers.

    Parameters
    ----------
    square_map : kanapy.texture.EBSDmap
        Two phases with two graph nodes and one edge each.
    """
    fig, ax = square_map.plot_graph_overlay(iphase=None, show=False)
    try:
        background = np.asarray(ax.images[0].get_array())
        np.testing.assert_array_equal(background[:, :10], square_map.ms_data[0]["rgb_im"][:, :10])
        np.testing.assert_array_equal(background[:, 10:], square_map.ms_data[1]["rgb_im"][:, 10:])
        markers = [line for line in ax.lines if line.get_marker() == "o"]
        edges = [line for line in ax.lines if len(line.get_xdata()) == 2]
        assert len(markers) == 4
        assert len(edges) == 2
        assert markers[0].get_color() == markers[1].get_color() == edges[0].get_color()
        assert markers[2].get_color() == markers[3].get_color() == edges[1].get_color()
        assert markers[0].get_color() != markers[2].get_color()
        np.testing.assert_allclose([line.get_xdata()[0] for line in markers], [2., 7., 12., 17.])
        assert [text.get_text() for text in ax.get_legend().get_texts()] == ["A (phase 2)", "B (phase 7)"]
    finally:
        plt.close(fig)


def test_multiphase_handoff_includes_all_retained_graphs(square_map, tmp_path):
    """
    Keep both phase graphs in optional files while preserving the selected RVE phase.

    Parameters
    ----------
    square_map : kanapy.texture.EBSDmap
        Analyzed map with two retained phases.
    tmp_path : pathlib.Path
        Temporary export directory.
    """
    result = EBSDGraphResult.from_phase_data(
        square_map.ms_data[0], 160, 8, 20, None, EBSDGraphConfig(), ebsd=square_map,
    )
    written = write_graph_result_outputs(result, EBSDGraphOutputOptions(
        output_dir=tmp_path, write_pickle=True, write_summary_json=True,
        write_graph_overlay=True, write_ipf_map=True,
    ))
    with written["pickle"].open("rb") as handle:
        payload = pickle.load(handle)
    assert payload["phase"]["index"] == 2
    assert [row["index"] for row in payload["phases"]] == [2, 7]
    assert [len(row["graph"]) for row in payload["phases"]] == [2, 2]
    summary = json.loads(written["summary_json"].read_text())
    assert [row["index"] for row in summary["phases"]] == [2, 7]
    assert all(path.stat().st_size > 0 for path in written.values())

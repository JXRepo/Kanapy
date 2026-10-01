"""
Generate Kanapy geometry and figures for the grain boundary proposal.

Controlled examples illustrate topology and candidate edits. They are not
measurements of the professor's unavailable failing microstructure.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
import pickle
import random
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.ndimage import generate_binary_structure, label

import kanapy
from kanapy.core.apd_geometry import build_grain_geometry
from kanapy.core.plotting import plot_ellipsoids_3D, plot_polygons_3D, plot_voxels_3D
from kanapy.core.power_diagram import AnisotropicPowerDiagram
from kanapy.core.voxelization import nonmanifold_grain_boundary_voxels

GREEN = "#739743"
BLUE = "#194A70"
AMBER = "#C99958"
GREY = "#DCE2E5"
RED = "#A73738"
PALETTE = [GREEN, AMBER, GREY, "#6698A6", "#AC918D", "#9FACC0", "#B7B28F"]
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 13,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": .8,
    "legend.frameon": False,
})


def save_figure(fig, name):
    """
    Export a PNG figure with a fixed canvas.

    A fixed canvas makes separately placed Beamer panels comparable.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Finished figure.
    name : str
        Basename for the PNG image.
    """
    fig.canvas.draw()
    destination = HERE / name
    fig.savefig(destination.with_suffix(".png"), dpi=300, facecolor="white")
    plt.close(fig)


def style_3d(fig, limits=(0., 1.), camera=(19, -58)):
    """
    Set consistent geometry views without decorative axes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure containing one three dimensional axes.
    limits : tuple of float
        Common lower and upper coordinate limits.
    camera : tuple of float
        Elevation and azimuth in degrees.
    """
    fig.set_size_inches(4.5, 4.1)
    ax = fig.axes[0]
    low, high = limits
    ax.set(xlim=(low, high), ylim=(low, high), zlim=(low, high))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(*camera)
    ax.set_title("")
    ax.set_axis_off()
    ax.set_position([0, 0, 1, 1])


def draw_grain(geometry, grain_id, name, limits=(0., 1.), mark=None):
    """
    Render one complete grain hull from Kanapy triangles.

    Parameters
    ----------
    geometry : dict
        Output of Kanapy's grain geometry builder.
    grain_id : int
        Grain to show in the common green colour.
    name : str
        Export basename.
    limits : tuple of float
        Common spatial limits.
    mark : array_like or None
        Optional point annotation on the existing geometry.
    """
    fig = plt.figure(figsize=(4.5, 4.1))
    ax = fig.add_subplot(111, projection="3d")
    triangles = np.asarray(geometry["Grains"][grain_id]["Simplices"])
    vertices = geometry["Points"][triangles]
    if limits != (0., 1.):
        low, high = limits
        inside = np.all((vertices >= low) & (vertices <= high), axis=(1, 2))
        vertices = vertices[inside]
    ax.add_collection3d(Poly3DCollection(
        vertices, facecolors=GREEN, edgecolors=(.13, .25, .12, .26),
        linewidths=.25, alpha=1.
    ))
    if mark is not None:
        ax.scatter(*mark, color=RED, s=48, depthshade=False)
    style_3d(fig, limits)
    save_figure(fig, name)


def vertex_link(geometry, grain_id, point):
    """
    Extract the triangle link around one grain surface vertex.

    Opposite triangle edges form the local link. Connected components and
    degrees distinguish a simple cycle from the controlled pinch.

    Parameters
    ----------
    geometry : dict
        Kanapy grain geometry with complete exterior closure.
    grain_id : int
        Grain whose boundary is checked.
    point : array_like
        Exact vertex coordinate used in the controlled example.

    Returns
    -------
    dict
        Vertex ID, link edges, connected components and degree information.
    """
    points = geometry["Points"]
    distances = np.linalg.norm(points - point, axis=1)
    vertex = int(np.argmin(distances))
    if distances[vertex] > 1e-9:
        raise ValueError("The requested point is not a Kanapy surface vertex")
    graph = defaultdict(set)
    edges = set()
    for triangle in geometry["Grains"][grain_id]["Simplices"]:
        if vertex not in triangle:
            continue
        a, b = [int(node) for node in triangle if node != vertex]
        edges.add(tuple(sorted((a, b))))
        graph[a].add(b)
        graph[b].add(a)
    remaining = set(graph)
    components = []
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            node = pending.pop()
            if node in component:
                continue
            component.add(node)
            pending.extend(graph[node] - component)
        remaining -= component
        components.append(sorted(component))
    return {
        "vertex": vertex, "edges": sorted(edges), "components": components,
        "degrees": {int(node): len(neighbours) for node, neighbours in graph.items()},
        "one_cycle": len(components) == 1 and all(len(n) == 2 for n in graph.values()),
    }


def link_figure(geometry, grain_id, point, name):
    """
    Draw the normalized triangle link derived from a Kanapy grain hull.

    Parameters
    ----------
    geometry : dict
        Kanapy geometry.
    grain_id : int
        Boundary label to inspect.
    point : array_like
        Surface vertex at the centre of the drawing.
    name : str
        Export basename.

    Returns
    -------
    dict
        Combinatorial link report.
    """
    report = vertex_link(geometry, grain_id, point)
    fig = plt.figure(figsize=(4.5, 4.1))
    ax = fig.add_subplot(111, projection="3d")
    u = np.linspace(0, 2*np.pi, 32)
    v = np.linspace(0, np.pi, 20)
    ax.plot_wireframe(np.outer(np.cos(u), np.sin(v)),
                      np.outer(np.sin(u), np.sin(v)),
                      np.outer(np.ones_like(u), np.cos(v)),
                      color=(.6, .65, .68, .18), linewidth=.4, rstride=3, cstride=3)
    directions = geometry["Points"] - point
    norms = np.linalg.norm(directions, axis=1)
    valid = norms > 1e-12
    directions[valid] /= norms[valid, None]
    for a, b in report["edges"]:
        segment = directions[[a, b]]
        ax.plot(*segment.T, color=GREEN, linewidth=3.2)
    style_3d(fig, (-1.1, 1.1), (16, -58))
    save_figure(fig, name)
    return report


def controlled_apds():
    """
    Generate the repository's analytical pinch and two weight controls.

    Changing the weight moves the whole interface. The neck control is a
    topology illustration and is not a demonstration of local repair.

    Returns
    -------
    dict
        Geometry, volume and vertex link measurements for the controls.
    """
    records = {}
    for name, weight in [("pinch", 0.), ("neck", .025), ("split", -.025)]:
        diagram = AnisotropicPowerDiagram(
            [[.5, .5, .5]] * 2,
            [np.eye(3), np.diag([2., 2., .5])], [1., 1., 1.]
        )
        diagram.weights[1] = weight
        print(f"Building controlled APD: {name}", flush=True)
        geometry = build_grain_geometry(diagram, {1: 0, 2: 0}, resolution=8)
        draw_grain(geometry, 2, f"{name}_surface",
                   mark=[.5, .5, .5] if name == "pinch" else None)
        records[name] = {
            "weight_difference": weight,
            "grain_volumes": {str(gid): g["Volume"] for gid, g in geometry["Grains"].items()},
            "background_resolution": 8,
            "surface_triangles": len(geometry["Surface"].triangles),
        }
        if name == "pinch":
            records[name]["vertex_link"] = link_figure(
                geometry, 2, np.full(3, .5), "link_pinch"
            )
            diagram_report = diagram.check_topology(
                resolution=4, sphere_level=1, max_candidates=60
            )
            records[name]["sampled_diagnostics"] = diagram_report.to_dict()
            assert not records[name]["vertex_link"]["one_cycle"]
            assert len(records[name]["vertex_link"]["components"]) == 2
            draw_grain(geometry, 2, "pinch_closeup", limits=(.32, .68), mark=[.5]*3)
    regular = AnisotropicPowerDiagram(
        [[.25, .5, .5], [.75, .5, .5]], [np.eye(3)]*2, [1., 1., 1.]
    )
    geometry = build_grain_geometry(regular, {1: 0, 2: 0}, resolution=4)
    records["regular_link"] = link_figure(geometry, 2, np.full(3, .5), "link_regular")
    assert records["regular_link"]["one_cycle"]

    four = AnisotropicPowerDiagram(
        [[.2, .2, .2], [.8, .2, .2], [.5, .8, .2], [.5, .5, .8]],
        np.tile(np.eye(3), (4, 1, 1)), [1., 1., 1.]
    )
    geometry = build_grain_geometry(four, {i: i-1 for i in range(1, 5)}, resolution=4)
    internal = geometry["Boundary"].triangulate(include_exterior=False)
    visible_geometry = dict(geometry, Surface=internal, Points=internal.points)
    fig = plot_polygons_3D(visible_geometry, silent=True, phases=True,
                           cols=[GREEN, AMBER, "#6698A6", "#AC918D"],
                           ec=(.13, .25, .3, .15))
    for collection in fig.axes[0].collections:
        collection.set_alpha(.45)
    style_3d(fig)
    for curve in geometry["Boundary"].junction_curves:
        xyz = geometry["Boundary"].points[curve.vertices]
        fig.axes[0].plot(*xyz.T, color=BLUE, linewidth=2.)
    junctions = geometry["Boundary"].points[geometry["Boundary"].junction_vertices]
    if len(junctions):
        fig.axes[0].scatter(*junctions.T, color=BLUE, s=38, depthshade=False)
    save_figure(fig, "normal_junction")
    records["ordinary_junction"] = geometry["Boundary"].summary()
    records["ordinary_junction"]["view"] = "Internal interfaces only; box faces omitted."
    return records


def voxel_metrics(grains):
    """
    Measure face connectivity and Kanapy's voxel boundary defects.

    Parameters
    ----------
    grains : ndarray
        Space filling array of grain labels.

    Returns
    -------
    dict
        Counts and complete local topology measurements for each grain.
    """
    bad = nonmanifold_grain_boundary_voxels(grains, periodic=False)
    records = {}
    for grain_id in np.unique(grains):
        mask = grains == grain_id
        _, count = label(mask, generate_binary_structure(3, 1))
        records[str(grain_id)] = {
            "voxels": int(mask.sum()), "face_components": int(count),
            "defect_incident_voxels": int(np.count_nonzero(bad & mask)),
        }
    return records


def local_candidates():
    """
    Construct joint candidate edits and an intentionally rejected edit.

    Labels A and B are shown; surrounding grain C is omitted from the three
    dimensional view so the internal contact remains visible. All labels,
    including C, participate in the Kanapy topology checks.

    Returns
    -------
    dict
        Candidate connectivity, topology, changed cells and volume changes.
    """
    original = np.full((11, 11, 11), 3, dtype=int)
    original[1:5, 1:5, 1:5] = 1
    original[5:9, 5:9, 5:9] = 1
    original[1:10, 4, 5] = 2
    original[1:3, 1:4, 5:8] = 2
    original[8:10, 1:4, 5:8] = 2
    greedy = original.copy()
    greedy[4:6, 4:6, 4:6] = 1
    joint = greedy.copy()
    joint[3:7, 3, 5] = 2
    records = {}
    baseline = voxel_metrics(original)
    for name, grains in [("before", original), ("greedy", greedy), ("joint", joint)]:
        metrics = voxel_metrics(grains)
        changes = {}
        for grain_id in metrics:
            old = baseline[grain_id]["voxels"]
            changes[grain_id] = (metrics[grain_id]["voxels"] - old) / old
        records[name] = {
            "grains": metrics,
            "changed_voxels": int(np.count_nonzero(grains != original)),
            "relative_volume_change": changes,
            "all_local_links_valid": all(g["defect_incident_voxels"] == 0 for g in metrics.values()),
            "all_grains_face_connected": all(g["face_components"] == 1 for g in metrics.values()),
        }
        np.save(HERE / f"local_{name}.npy", grains)
        fig = plot_voxels_3D(grains, mask=grains != 3, silent=True,
                             clist=np.array([to_rgb(GREEN), to_rgb(AMBER), to_rgb(GREY)]))
        style_3d(fig, (0, 11), (24, -65))
        if name == "before":
            fig.axes[0].scatter(5., 5., 5., color=RED, s=36, depthshade=False)
        save_figure(fig, f"local_{name}")
    fig = plot_voxels_3D(original, mask=original != 3, silent=True,
                         clist=np.array([to_rgb(GREEN), to_rgb(AMBER), to_rgb(GREY)]))
    style_3d(fig, (0, 11), (24, -65))
    ax = fig.axes[0]
    corners = np.array(list(product([3., 7.], repeat=3)))
    for i, first in enumerate(corners):
        for second in corners[i+1:]:
            if np.count_nonzero(first != second) == 1:
                ax.plot(*np.stack([first, second]).T, color=RED, linewidth=1.4, linestyle=":")
    # Show the two interior locations through the surface as projected markers.
    # The grain geometry is unchanged; ordinary 3D scatter would hide these cores.
    for point in ((2.5, 2.5, 2.5), (7.5, 7.5, 7.5)):
        x_screen, y_screen, _ = proj3d.proj_transform(*point, ax.get_proj())
        ax.add_artist(Line2D([x_screen], [y_screen], linestyle="none", marker="o",
                             markersize=6, markerfacecolor=BLUE, markeredgecolor="white",
                             markeredgewidth=.7, transform=ax.transData, zorder=10000))
    save_figure(fig, "local_band")
    assert records["before"]["grains"]["1"]["face_components"] == 2
    assert records["greedy"]["grains"]["2"]["face_components"] == 2
    assert records["joint"]["all_local_links_valid"]
    assert records["joint"]["all_grains_face_connected"]
    assert records["joint"]["changed_voxels"] == 10
    assert [records["joint"]["grains"][str(g)]["voxels"] for g in (1, 2, 3)] == [134, 47, 1150]
    return records


def contextual_rve():
    """
    Generate one small packed RVE through the current Kanapy API.

    This adapts the statistical example with fixed seeds and a smaller domain.
    Its pictures provide context, not evidence that the proposed repair works.

    Returns
    -------
    dict
        Saved pipeline parameters and geometry counts.
    """
    cached = HERE / "context_rve_low_fill.pkl"
    metadata_path = HERE / "context_generation.json"
    if cached.is_file():
        with cached.open("rb") as handle:
            ms = pickle.load(handle)
    else:
        random.seed(20261001)
        np.random.seed(20261001)
        descriptor = {
            "Grain type": "Elongated",
            "Equivalent diameter": {"sig": .25, "scale": 10., "loc": 0.,
                                    "cutoff_min": 7., "cutoff_max": 14.},
            "Aspect ratio": {"sig": .25, "scale": 1.7, "loc": 0.,
                             "cutoff_min": 1., "cutoff_max": 2.7},
            "Tilt angle": {"kappa": 1., "loc": .5*np.pi,
                           "cutoff_min": 0., "cutoff_max": np.pi},
            "RVE": {"sideX": 24., "sideY": 24., "sideZ": 24.,
                    "Nx": 20, "Ny": 20, "Nz": 20},
            "Simulation": {"periodicity": False, "output_units": "um"},
        }
        print("Generating contextual packed RVE", flush=True)
        ms = kanapy.Microstructure(descriptor=descriptor, name="proposal_context")
        ms.init_RVE(nsteps=250)
        ms.pack(fill_factor=.25, relaxation_steps=500)
        ms.voxelize(fit_options={"n_samples": 16384, "seed": 20261001, "maxiter": 200})
        ms.generate_grains(resolution=6)
        with cached.open("wb") as handle:
            pickle.dump(ms, handle)
        (HERE / "context_descriptor.json").write_text(
            json.dumps(descriptor, indent=2), encoding="utf-8"
        )
        metadata_path.write_text(json.dumps({
            "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "kanapy_version": kanapy.__version__,
            "random_seed": 20261001,
        }, indent=2), encoding="utf-8")
    diagram = ms.mesh.apd
    np.savez(HERE / "context_apd.npz", centers=diagram.centers,
             matrices=diagram.matrices, weights=diagram.weights,
             box_size=diagram.box_size, grain_ids=diagram.grain_ids,
             target_volumes=diagram.target_volumes, voxel_labels=ms.mesh.grains)
    particle_colors = [GREY]
    particle_colors.extend(PALETTE[i % len(PALETTE)] for i in range(len(ms.particles)))
    cmap = ListedColormap(particle_colors)
    fig = plot_ellipsoids_3D(ms.particles, cmap=cmap, silent=True)
    style_3d(fig, (0, 24), (24, -58))
    save_figure(fig, "pipeline_ellipsoids")
    fig = plot_voxels_3D(ms.mesh.grains, silent=True, sliced=True,
                         clist=np.array([to_rgb(PALETTE[(int(g)-1) % len(PALETTE)])
                                         for g in np.unique(ms.mesh.grains)]))
    style_3d(fig, (0, 20), (24, -58))
    save_figure(fig, "pipeline_voxels")
    surface = ms.geometry["Surface"]
    indices = []
    for i, triangle in enumerate(surface.triangles):
        center = surface.points[triangle].mean(axis=0)
        if not np.all(center > 12.):
            indices.append(i)
    fig = plt.figure(figsize=(4.5, 4.1))
    ax = fig.add_subplot(111, projection="3d")
    colors = [PALETTE[(surface.face_grains[i][0]-1) % len(PALETTE)] for i in indices]
    ax.add_collection3d(Poly3DCollection(
        surface.points[surface.triangles[indices]], facecolors=colors,
        edgecolors=(.1, .2, .25, .22), linewidths=.23
    ))
    style_3d(fig, (0, 24), (24, -58))
    save_figure(fig, "pipeline_surface")
    original_count = sum(p.duplicate is None for p in ms.particles)
    return {
        "random_seed": 20261001, "box_um": [24, 24, 24],
        "generation_metadata": json.loads(metadata_path.read_text()) if metadata_path.is_file() else {
            "git_commit": None,
            "note": "Legacy cache has no generation metadata; regenerate it for a complete provenance record.",
        },
        "packing_fill_factor": .25, "packing_relaxation_limit": 500,
        "periodic": False, "original_particles": original_count,
        "voxel_shape": list(ms.mesh.grains.shape),
        "surface_background_resolution": 6,
        "surface_triangles": len(surface.triangles),
        "grain_volumes": {str(g): v["Volume"] for g, v in ms.geometry["Grains"].items()},
        "packing_relaxation": getattr(ms.simbox, "packing_relaxation", None),
        "visual_crop": "Upper coordinate octant omitted in voxel and surface views; source remains complete.",
    }


def periodic_figure():
    """
    Show a grain crossing the periodic box boundary using Kanapy labels.

    Returns
    -------
    dict
        Source diagram parameters for the slice.
    """
    diagram = AnisotropicPowerDiagram(
        [[.05, .5, .5], [.55, .5, .5]], [np.eye(3)]*2,
        [1., 1., 1.], periodic=True
    )
    edges = np.linspace(0, 1, 302)
    centers = (edges[:-1] + edges[1:]) / 2
    x, y = np.meshgrid(centers, centers, indexing="xy")
    points = np.column_stack((x.ravel(), y.ravel(), np.full(x.size, .5)))
    grains = diagram.labels(points).reshape(x.shape)
    fig, ax = plt.subplots(figsize=(6.5, 3.7))
    ax.pcolormesh(edges, edges, grains, cmap=ListedColormap([GREEN, AMBER]),
                  vmin=1, vmax=2, shading="flat", rasterized=True)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="x / L", ylabel="y / L")
    ax.set_xticks([0, .3, .8, 1])
    ax.set_yticks([0, 1])
    ax.text(.15, .5, "Grain A", ha="center", va="center", color="white", fontsize=14)
    ax.text(.55, .5, "Grain B", ha="center", va="center", color="white", fontsize=14)
    ax.text(.9, .5, "A", ha="center", va="center", color="white", fontsize=14)
    fig.subplots_adjust(left=.11, bottom=.19, right=.95, top=.94)
    save_figure(fig, "periodic_slice")
    return {"centers": diagram.centers.tolist(), "periodic": True, "slice_z": .5}


def json_default(value):
    """
    Convert NumPy values in the provenance report to JSON values.

    Parameters
    ----------
    value : object
        Object requested by the JSON encoder.

    Returns
    -------
    object
        Serializable representation.
    """
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if hasattr(value, "__dict__"):
        return vars(value)
    raise TypeError(f"Cannot serialize {type(value)}")


def main():
    """
    Generate selected figure sets and save their measurements.

    Independent figure sets can be regenerated without repeating packing.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=["all", "controls", "local", "context", "periodic"], default="all")
    args = parser.parse_args()
    report_path = HERE / "figure_provenance.json"
    records = json.loads(report_path.read_text()) if report_path.is_file() else {}
    records["kanapy_version"] = kanapy.__version__
    records["kanapy_source"] = str(Path(kanapy.__file__).resolve().relative_to(ROOT))
    records["git_commit"] = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    for option, function in [("controls", controlled_apds), ("local", local_candidates),
                              ("context", contextual_rve), ("periodic", periodic_figure)]:
        if args.only in ["all", option]:
            records[option] = function()
            report_path.write_text(json.dumps(records, indent=2, default=json_default), encoding="utf-8")
    print("Figure generation finished", flush=True)


if __name__ == "__main__":
    main()

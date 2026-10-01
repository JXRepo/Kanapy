"""
Draw a generic tetrahedral volume mesh for the presentation workflow.

This illustration does not use Kanapy geometry and is not a meshed version of
the neighbouring RVE. The cube has one material and no grain boundaries.
"""

from collections import defaultdict
from itertools import combinations, permutations, product
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


HERE = Path(__file__).resolve().parent


def cube_mesh(divisions=4):
    """
    Partition a unit cube into conforming tetrahedra.

    Every regular cell uses the same six tetrahedra around its ascending body
    diagonal. Neighbouring cells therefore share matching triangular faces.

    Parameters
    ----------
    divisions : int, optional
        Number of regular cells along each coordinate axis.

    Returns
    -------
    points : numpy.ndarray
        Shared vertex coordinates, with shape (n, 3).
    tetrahedra : numpy.ndarray
        Positively oriented vertex indices, with shape (m, 4).
    """
    grid = list(product(range(divisions + 1), repeat=3))
    indices = {coordinate: index for index, coordinate in enumerate(grid)}
    points = np.asarray(grid, dtype=float) / divisions
    tetrahedra = []

    for cell in product(range(divisions), repeat=3):
        for order in permutations(range(3)):
            corner = list(cell)
            tet = [indices[tuple(corner)]]
            for axis in order:
                corner[axis] += 1
                tet.append(indices[tuple(corner)])
            vertices = points[tet]
            if np.linalg.det(vertices[1:] - vertices[0]) < 0:
                tet[1], tet[2] = tet[2], tet[1]
            tetrahedra.append(tet)

    return points, np.asarray(tetrahedra, dtype=int)


def validate_mesh(points, tetrahedra):
    """
    Verify the original cube mesh before changing its display.

    Signed volumes and shared faces establish positive volume, complete cube
    coverage, and connectivity through faces for this regular construction.

    Parameters
    ----------
    points : numpy.ndarray
        Original shared vertex coordinates.
    tetrahedra : numpy.ndarray
        Tetrahedral vertex indices.

    Returns
    -------
    dict
        Counts and geometric checks recorded alongside the figure.

    Raises
    ------
    ValueError
        If volume, face incidence, boundary position, or connectivity fails.
    """
    vertices = points[tetrahedra]
    volumes = np.linalg.det(vertices[:, 1:] - vertices[:, :1]) / 6
    if not np.all(volumes > 0) or not np.isclose(volumes.sum(), 1.0):
        raise ValueError("Tetrahedra must have positive volume and fill the cube.")

    face_owners = defaultdict(list)
    for tet_id, tet in enumerate(tetrahedra):
        for face in combinations(tet, 3):
            face_owners[tuple(sorted(face))].append(tet_id)

    neighbours = [set() for _ in tetrahedra]
    boundary_count = 0
    for face, owners in face_owners.items():
        if len(owners) == 1:
            coordinates = points[list(face)]
            on_low = np.all(np.isclose(coordinates, 0), axis=0)
            on_high = np.all(np.isclose(coordinates, 1), axis=0)
            if not np.any(on_low | on_high):
                raise ValueError("An unmatched triangle lies inside the cube.")
            boundary_count += 1
        elif len(owners) == 2:
            left, right = owners
            neighbours[left].add(right)
            neighbours[right].add(left)
        else:
            raise ValueError("More than two tetrahedra share a face.")

    visited = {0}
    pending = [0]
    while pending:
        current = pending.pop()
        for neighbour in neighbours[current] - visited:
            visited.add(neighbour)
            pending.append(neighbour)
    if len(visited) != len(tetrahedra):
        raise ValueError("The mesh is not connected through shared faces.")

    return {
        "nodes": len(points),
        "tetrahedra": len(tetrahedra),
        "minimum_signed_volume": float(volumes.min()),
        "total_volume": float(volumes.sum()),
        "all_signed_volumes_positive": True,
        "boundary_triangles": boundary_count,
        "interior_shared_triangles": len(face_owners) - boundary_count,
        "all_unmatched_faces_on_cube_boundary": True,
        "connected_through_shared_faces": True,
    }


def mesh_demo():
    """
    Save a cutaway illustration and its explicit source record.

    A corner is hidden and the remaining tetrahedra are slightly shrunk around
    their centres to reveal volume elements. These changes affect display only;
    the validated mesh retains shared vertices and fills the whole unit cube.
    """
    divisions = 4
    shrink = 0.82
    points, tetrahedra = cube_mesh(divisions)
    checks = validate_mesh(points, tetrahedra)
    vertices = points[tetrahedra]
    centres = vertices.mean(axis=1)
    hidden = (centres[:, 0] > 0.5) & (centres[:, 1] < 0.5) & (centres[:, 2] > 0.5)
    visible_vertices = vertices[~hidden]
    visible_centres = centres[~hidden]
    display_vertices = visible_centres[:, None] + shrink * (
        visible_vertices - visible_centres[:, None]
    )

    faces = []
    colors = []
    base_color = np.asarray(to_rgb("#8CA8BA"))
    light = np.asarray([-0.5, -0.8, 1.6])
    light /= np.linalg.norm(light)
    for tet in display_vertices:
        centre = tet.mean(axis=0)
        for indices in combinations(range(4), 3):
            face = tet[list(indices)]
            normal = np.cross(face[1] - face[0], face[2] - face[0])
            normal /= np.linalg.norm(normal)
            if np.dot(normal, face.mean(axis=0) - centre) < 0:
                normal = -normal
            brightness = 0.72 + 0.27 * max(0, float(np.dot(normal, light)))
            faces.append(face)
            colors.append(base_color * brightness)

    fig = plt.figure(figsize=(4.5, 4.1))
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(
        faces,
        facecolors=colors,
        edgecolors="#344C5B",
        linewidths=0.48,
        antialiased=True,
    ))
    ax.set(xlim=(0, 1), ylim=(0, 1), zlim=(0, 1))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(24, -58)
    ax.set_axis_off()
    ax.set_position([0, 0, 1, 1])
    fig.savefig(HERE / "pipeline_volume_mesh.png", dpi=300, facecolor="white")
    plt.close(fig)

    provenance = {
        "figure": "pipeline_volume_mesh.png",
        "function": "mesh_demo()",
        "script": "mesh_demo.py",
        "classification": "Generic volume mesh illustration",
        "source": "Regular unit cube, generated by this report script; not Kanapy",
        "relationship_to_RVE": (
            "Not the RVE in the other pipeline panels. No grain boundaries are "
            "present, so this does not demonstrate a grain boundary conforming RVE mesh."
        ),
        "construction": {
            "cells_per_axis": divisions,
            "tetrahedra_per_cell": 6,
            "method": "Consistent ascending body diagonal subdivision",
            "material_count": 1,
        },
        "display_only": {
            "cutaway": "Hide tetrahedra whose centres satisfy x > 0.5, y < 0.5, z > 0.5",
            "hidden_tetrahedra": int(hidden.sum()),
            "visible_tetrahedra": int((~hidden).sum()),
            "shrink_about_each_tetrahedron_centre": shrink,
            "note": "Cutaway and gaps are display changes, not holes in the validated mesh.",
        },
        "validation_before_display_changes": checks,
    }
    destination = HERE / "pipeline_volume_mesh.json"
    destination.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    mesh_demo()

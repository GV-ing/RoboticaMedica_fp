#!/usr/bin/env python3
# tools/generate_meshes.py
"""Genera in meshes/ le mesh del progetto con i parametri nominali di modules/config.py.

Sono le stesse funzioni che usa la scena (modules/geometry_gen.py), che comunque rigenera
le mesh a ogni avvio: questo script serve a ottenerle senza SOFA, da ispezionare o consegnare.

  meshes/tube_surface.obj      superficie del tubo, modello visivo         (STL/OBJ, traccia)
  meshes/tube_collision.obj    superficie del tubo, modello di collisione  (STL/OBJ, traccia)
  meshes/tube_fem_surface.obj  superficie esterna della mesh FEM (controllo)
  meshes/tube_volume.vtk       mesh volumetrica tetraedrica del tubo      (meccanica, traccia)
  meshes/soft_volume.vtk       mesh volumetrica esaedrica dell'oggetto soffice
  meshes/object_surface.obj    superficie dell'oggetto di contatto (rigido e soffice)

I .vtk si aprono con ParaView. Non richiede SOFA (solo numpy):
    python3.12 tools/generate_meshes.py          # dopo source env.sh
    python3 tools/generate_meshes.py             # oppure Python di sistema con numpy
"""

import os
import sys

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

from modules.config import NOMINAL_PARAMS  # noqa: E402
from modules.geometry_gen import (  # noqa: E402
    box_hexa_mesh, box_surface, build_tube_meshes, check_tet_mesh, contact_object_axes,
    write_vtk_hexa, write_vtk_tetra,
)


def write_obj(path, positions, triangles):
    with open(path, "w", encoding="ascii") as f:
        f.write("# Oggetto di contatto - superficie\n")
        for p in positions:
            f.write(f"v {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        for t in triangles:
            f.write(f"f {t[0] + 1} {t[1] + 1} {t[2] + 1}\n")


def main():
    p = dict(NOMINAL_PARAMS)
    mesh_dir = os.path.join(PROJECT_DIR, "meshes")

    nodes, tetras, files = build_tube_meshes(p, mesh_dir)
    bad, vmin, vol = check_tet_mesh(nodes, tetras)
    files["volume"] = os.path.join(mesh_dir, "tube_volume.vtk")
    write_vtk_tetra(files["volume"], nodes, tetras)
    print(f"Tubo: {len(nodes)} nodi, {len(tetras)} tetraedri (h_flex = {p['mesh_h_flex']} mm), "
          f"invertiti/degeneri = {bad} (atteso 0), volume = {vol:.2f} mm^3")

    xs, ys, zs, _, _ = contact_object_axes(dict(p, contact_object="soft"))
    soft_nodes, soft_hexa = box_hexa_mesh(xs, ys, zs)
    files["soft_volume"] = os.path.join(mesh_dir, "soft_volume.vtk")
    write_vtk_hexa(files["soft_volume"], soft_nodes, soft_hexa)
    surf_pos, surf_tri = box_surface(xs, ys, zs)
    files["object_surface"] = os.path.join(mesh_dir, "object_surface.obj")
    write_obj(files["object_surface"], surf_pos, surf_tri)
    print(f"Oggetto soffice ({p['soft_mesh']}): {len(soft_nodes)} nodi, {len(soft_hexa)} esaedri")

    for name, path in files.items():
        print(f"  {name:15s} {os.path.relpath(path, PROJECT_DIR)}")
    if bad:
        sys.exit("Mesh del tubo non valida")


if __name__ == "__main__":
    main()

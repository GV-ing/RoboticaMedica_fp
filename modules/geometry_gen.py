# modules/geometry_gen.py
import os
import itertools
import numpy as np


def create_hollow_tube_meshes(
    obj_filepath,
    od=1.93,          # Diametro esterno [mm]
    thickness=0.22,   # Spessore parete [mm] -> ID = 1.49 mm
    length=100.0,     # Lunghezza totale [mm] (40 mm albero + 60 mm flessibile)
    n_x=30,           # Suddivisioni longitudinali
    n_theta=12        # Suddivisioni circonferenziali
):
    """
    Genera la mesh volumetrica tetraedrica (nodi + tetraedri) per TetrahedronFEMForceField
    e la mesh di superficie (.obj) per VisualModel e CollisionModel.

    Ogni cella (settore anulare tra ix e ix+1, ith e ith+1) e' un esaedro con coordinate
    locali (a, b, r): a = asse x, b = theta, r = 0 interno / 1 esterno.
    Viene scomposta con lo schema di Freudenthal: 6 tetraedri, uno per ciascuna permutazione
    degli assi, lungo i cammini da (0,0,0) a (1,1,1). Lo schema e' conforme tra celle adiacenti
    (diagonali di faccia coerenti) e non produce elementi degeneri.
    """
    if od <= 0 or thickness <= 0 or thickness >= od / 2:
        raise ValueError("Richiesti od > 0 e 0 < thickness < od/2")
    if length <= 0 or n_x < 1 or n_theta < 3:
        raise ValueError("Richiesti length > 0, n_x >= 1 e n_theta >= 3")

    r_out = od / 2.0
    r_in = r_out - thickness

    x_vals = np.linspace(0.0, length, n_x + 1)
    theta_vals = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)

    fem_nodes = []
    for ix, x in enumerate(x_vals):
        for ith, th in enumerate(theta_vals):
            # Nodo interno (r_in)
            fem_nodes.append([float(x), float(r_in * np.cos(th)), float(r_in * np.sin(th))])
            # Nodo esterno (r_out)
            fem_nodes.append([float(x), float(r_out * np.cos(th)), float(r_out * np.sin(th))])

    def node_idx(ix, ith, is_outer):
        return ix * (2 * n_theta) + (ith % n_theta) * 2 + (1 if is_outer else 0)

    nodes_np = np.array(fem_nodes)

    tetrahedra = []
    for ix in range(n_x):
        for ith in range(n_theta):
            def c(a, b, r):
                return node_idx(ix + a, ith + b, bool(r))

            for perm in itertools.permutations(range(3)):
                p = [0, 0, 0]
                path = [c(*p)]
                for ax in perm:
                    p[ax] = 1
                    path.append(c(*p))
                P = nodes_np[path]
                vol6 = np.dot(P[1] - P[0], np.cross(P[2] - P[0], P[3] - P[0]))
                if vol6 < 0:
                    path[1], path[2] = path[2], path[1]   # orientamento positivo
                tetrahedra.append([int(i) for i in path])

    # Generazione file .obj di superficie per VisualModel e CollisionModel
    output_dir = os.path.dirname(os.path.abspath(obj_filepath))
    os.makedirs(output_dir, exist_ok=True)
    faces = []
    for ix in range(n_x):
        for ith in range(n_theta):
            i_next = (ith + 1) % n_theta
            v_in_1  = node_idx(ix,     ith,    False) + 1
            v_in_2  = node_idx(ix,     i_next, False) + 1
            v_in_3  = node_idx(ix + 1, i_next, False) + 1
            v_in_4  = node_idx(ix + 1, ith,    False) + 1

            v_out_1 = node_idx(ix,     ith,    True) + 1
            v_out_2 = node_idx(ix,     i_next, True) + 1
            v_out_3 = node_idx(ix + 1, i_next, True) + 1
            v_out_4 = node_idx(ix + 1, ith,    True) + 1

            # Superficie esterna
            faces.append((v_out_1, v_out_2, v_out_3))
            faces.append((v_out_1, v_out_3, v_out_4))
            # Superficie interna
            faces.append((v_in_1, v_in_3, v_in_2))
            faces.append((v_in_1, v_in_4, v_in_3))

    # Chiude le due estremita' dell'anello, mantenendo vuoto il lume del tubo.
    for ith in range(n_theta):
        i_next = (ith + 1) % n_theta
        start_inner = node_idx(0, ith, False) + 1
        start_inner_next = node_idx(0, i_next, False) + 1
        start_outer = node_idx(0, ith, True) + 1
        start_outer_next = node_idx(0, i_next, True) + 1
        faces.extend([
            (start_inner, start_inner_next, start_outer_next),
            (start_inner, start_outer_next, start_outer),
        ])

        end_inner = node_idx(n_x, ith, False) + 1
        end_inner_next = node_idx(n_x, i_next, False) + 1
        end_outer = node_idx(n_x, ith, True) + 1
        end_outer_next = node_idx(n_x, i_next, True) + 1
        faces.extend([
            (end_inner, end_outer_next, end_inner_next),
            (end_inner, end_outer, end_outer_next),
        ])

    with open(obj_filepath, "w") as f:
        f.write(f"# Hollow Tube Surface Mesh (Tetra) - OD={od}mm, t={thickness}mm, L={length}mm\n")
        for pt in fem_nodes:
            f.write(f"v {pt[0]:.6f} {pt[1]:.6f} {pt[2]:.6f}\n")
        for fc in faces:
            f.write(f"f {fc[0]} {fc[1]} {fc[2]}\n")

    return fem_nodes, tetrahedra


def check_tet_mesh(nodes, tetrahedra):
    """Diagnostica: conta tetraedri invertiti/degeneri e restituisce il volume minimo [mm^3]."""
    n = np.array(nodes)
    t = np.array(tetrahedra)
    v = np.einsum("ij,ij->i", n[t[:, 1]] - n[t[:, 0]],
                  np.cross(n[t[:, 2]] - n[t[:, 0]], n[t[:, 3]] - n[t[:, 0]])) / 6.0
    return int((v <= 1e-12).sum()), float(v.min()), float(v.sum())


if __name__ == "__main__":
    out_obj = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "meshes", "tube_surface.obj")
    nodes, tetras = create_hollow_tube_meshes(out_obj)
    bad, vmin, vtot = check_tet_mesh(nodes, tetras)
    print(f"Mesh tetraedrica generata! Nodi FEM: {len(nodes)}, Tetraedri: {len(tetras)}")
    print(f"Tetraedri invertiti/degeneri: {bad} (atteso 0) | volume min: {vmin:.4f} mm^3 | volume tot: {vtot:.3f} mm^3")
    print(f"File OBJ salvato in: {os.path.abspath(out_obj)}")
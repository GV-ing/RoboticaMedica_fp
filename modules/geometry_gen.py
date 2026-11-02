# modules/geometry_gen.py
import os
import itertools
import numpy as np


def create_hollow_tube_meshes(
    obj_filepath,
    od=1.93,          # Diametro esterno [mm]
    thickness=0.22,   # Spessore parete [mm] -> ID = 1.49 mm
    length=100.0,     # Lunghezza totale [mm] (40 mm albero + 60 mm flessibile)
    n_x=30,           # Suddivisioni longitudinali (se x_vals non e' dato)
    n_theta=12,       # Suddivisioni circonferenziali
    x_vals=None,      # Coordinate assiali degli anelli [mm] (mesh graduata); prevale su n_x/length
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
    if x_vals is None:
        if length <= 0 or n_x < 1:
            raise ValueError("Richiesti length > 0 e n_x >= 1")
        x_vals = np.linspace(0.0, length, n_x + 1)
    else:
        x_vals = np.asarray(x_vals, dtype=float)
        if x_vals.ndim != 1 or x_vals.size < 2 or np.any(np.diff(x_vals) <= 0):
            raise ValueError("x_vals deve essere un vettore strettamente crescente (>= 2 valori)")
        n_x = x_vals.size - 1
        length = float(x_vals[-1] - x_vals[0])
    if n_theta < 3:
        raise ValueError("Richiesto n_theta >= 3")

    r_out = od / 2.0
    r_in = r_out - thickness
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


def graded_x_vals(length_total, length_shaft, h_flex, n_shaft=4):
    """Anelli assiali: n_shaft elementi sull'albero (vincolato) e passo ~h_flex sulla parte flessibile."""
    if not 0 < length_shaft < length_total:
        raise ValueError("Richiesto 0 < length_shaft < length_total")
    if h_flex <= 0 or n_shaft < 1:
        raise ValueError("Richiesti h_flex > 0 e n_shaft >= 1")
    n_flex = max(1, int(round((length_total - length_shaft) / h_flex)))
    shaft = np.linspace(0.0, length_shaft, n_shaft + 1)
    flex = np.linspace(length_shaft, length_total, n_flex + 1)[1:]
    return np.concatenate((shaft, flex))


def write_surface_obj(obj_filepath, od, thickness, x_start, x_end, h, n_theta=12):
    """Superficie chiusa del tubo tra x_start e x_end (passo ~h) per modelli visivo/collisione.

    I vertici giacciono sulla superficie della mesh FEM (stessi angoli theta), quindi la
    BarycentricMapping verso la mesh fine e' esatta. Restituisce il numero di vertici.
    """
    n = max(1, int(round((x_end - x_start) / h)))
    nodes, _ = create_hollow_tube_meshes(
        obj_filepath, od=od, thickness=thickness, n_theta=n_theta,
        x_vals=np.linspace(x_start, x_end, n + 1))
    return len(nodes)


def graded_axis(lo, hi, focus, h_min, h_max, growth=1.4):
    """Coordinate in [lo, hi] con passo h_min attorno a `focus`, crescente in progressione
    geometrica (fattore `growth`) fino a h_max. Contiene lo, hi e focus (se interno)."""
    if not lo < hi or not 0 < h_min <= h_max or growth < 1.0:
        raise ValueError("Richiesti lo < hi, 0 < h_min <= h_max, growth >= 1")
    focus = min(max(focus, lo), hi)

    def side(length):
        pts, x, h = [], 0.0, h_min
        while length - x > 1e-12:
            step = min(h, h_max)
            if length - x - step < 0.5 * step:   # ultimo intervallo assorbe il residuo
                step = length - x
            x += step
            pts.append(x)
            h *= growth
        return pts

    right = [focus + d for d in side(hi - focus)]
    left = [focus - d for d in side(focus - lo)]
    xs = sorted(set([lo, hi, focus] + right + left))
    out = [xs[0]]
    for a, b in zip(xs[:-1], xs[1:]):          # nessun intervallo oltre h_max
        n = int(np.ceil((b - a) / h_max - 1e-9))
        out.extend(np.linspace(a, b, n + 1)[1:].tolist())
    return np.clip(np.array(out, dtype=float), lo, hi)


def box_hexa_mesh(xs, ys, zs):
    """Griglia strutturata di esaedri su coordinate (anche non uniformi).

    Ordine dei nodi: indice = i + nx*(j + ny*k); ordine dei vertici dell'esaedro come in
    RegularGridTopology di SOFA: (i,j,k) (i+1,j,k) (i+1,j+1,k) (i,j+1,k) e poi a k+1.
    """
    nx, ny, nz = len(xs), len(ys), len(zs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pos = np.stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")], axis=1)

    def idx(i, j, k):
        return i + nx * (j + ny * k)

    hexa = []
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                hexa.append([idx(i, j, k), idx(i + 1, j, k), idx(i + 1, j + 1, k), idx(i, j + 1, k),
                             idx(i, j, k + 1), idx(i + 1, j, k + 1), idx(i + 1, j + 1, k + 1),
                             idx(i, j + 1, k + 1)])
    return pos.tolist(), hexa


def box_surface(xs, ys, zs):
    """Superficie triangolata (normali uscenti) delle 6 facce di una griglia strutturata.

    Restituisce (posizioni, triangoli); i vertici coincidono con i nodi di bordo della griglia.
    """
    nx, ny, nz = len(xs), len(ys), len(zs)
    index, positions = {}, []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if i in (0, nx - 1) or j in (0, ny - 1) or k in (0, nz - 1):
                    index[(i, j, k)] = len(positions)
                    positions.append([float(xs[i]), float(ys[j]), float(zs[k])])
    tris = []

    def quad(c):
        a, b, cc, d = (index[x] for x in c)
        tris.extend([[a, b, cc], [a, cc, d]])

    for j in range(ny - 1):
        for k in range(nz - 1):
            quad([(0, j, k), (0, j, k + 1), (0, j + 1, k + 1), (0, j + 1, k)])
            quad([(nx - 1, j, k), (nx - 1, j + 1, k), (nx - 1, j + 1, k + 1), (nx - 1, j, k + 1)])
    for i in range(nx - 1):
        for k in range(nz - 1):
            quad([(i, 0, k), (i + 1, 0, k), (i + 1, 0, k + 1), (i, 0, k + 1)])
            quad([(i, ny - 1, k), (i, ny - 1, k + 1), (i + 1, ny - 1, k + 1), (i + 1, ny - 1, k)])
    for i in range(nx - 1):
        for j in range(ny - 1):
            quad([(i, j, 0), (i, j + 1, 0), (i + 1, j + 1, 0), (i + 1, j, 0)])
            quad([(i, j, nz - 1), (i + 1, j, nz - 1), (i + 1, j + 1, nz - 1), (i, j + 1, nz - 1)])
    return positions, tris


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
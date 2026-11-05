"""Lettura della forza di contatto dalle forze dei vincoli del constraint solver.

Perche' non il campo ``lambda`` dei MechanicalObject: sui modelli di collisione mappati
(oggetto soffice, tubo) ``lambda`` puo' contenere la reazione contata due volte
(verificato durante lo sviluppo: rapporto esattamente 2.00
rispetto alla forza ricavata dalla cinematica). Le ``constraintForces`` del solver sono
invece la sorgente primaria e coincidono con il riferimento entro lo 0.5% sia sul rigido
(equilibrio) sia sul soffice (cinematica).

Convenzioni verificate in SOFA 26.06 con FreeMotionAnimationLoop:
  * i valori dei vincoli sono IMPULSI [N*s]: forza = valore / dt;
  * il CableConstraint occupa una riga, con valore uguale all'impulso imposto (T*dt);
  * ogni punto di contatto con FrictionContactConstraint occupa 3 righe (n, t1, t2):
    la forza normale e' la prima del blocco.
"""

import numpy as np

SOLVER_NAME = "ConstraintSolver"
CABLE_PATH = ("NitinolTube", "Tendon", "tendonCable")


def _find_cable(root):
    node = root
    for name in CABLE_PATH[:-1]:
        node = node.getChild(name) if node is not None else None
    return node.getObject(CABLE_PATH[-1]) if node is not None else None


def _cable_row(forces, cable_impulse):
    """Indice della riga del cavo: quella uguale all'impulso imposto, altrimenti l'ultima."""
    if cable_impulse is not None and abs(cable_impulse) > 0.0:
        diff = np.abs(forces - cable_impulse)
        idx = int(np.argmin(diff))
        if diff[idx] <= 1e-9 + 1e-6 * abs(cable_impulse):
            return idx
    return forces.size - 1


def contact_normal_impulses(root):
    """Impulsi normali di contatto [N*s] (uno per punto di contatto)."""
    solver = root.getObject(SOLVER_NAME)
    if solver is None:
        raise RuntimeError(f"Constraint solver '{SOLVER_NAME}' non trovato")
    forces = np.asarray(solver.constraintForces.value, dtype=float).ravel()
    if forces.size == 0:
        return forces
    if not np.all(np.isfinite(forces)):
        raise ValueError("constraintForces contiene valori non finiti")

    cable = _find_cable(root)
    has_cable = cable is not None
    if has_cable:
        impulse = float(np.ravel(cable.value.value)[0])
        forces = np.delete(forces, _cable_row(forces, impulse))
    if forces.size % 3 != 0:
        raise RuntimeError(
            f"Righe di contatto non multiple di 3 ({forces.size}, cavo "
            f"{'escluso' if has_cable else 'assente'}): verificare il tipo di risposta"
        )
    return forces[0::3]


def read_contact_force(root):
    """Forza normale di contatto totale [N] (somma delle normali / dt)."""
    dt = float(root.dt.value)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError(f"dt non valido: {dt}")
    normals = contact_normal_impulses(root)
    if normals.size == 0:
        return 0.0
    return float(abs(np.sum(normals)) / dt)

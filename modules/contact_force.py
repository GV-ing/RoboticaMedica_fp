"""Read contact reaction forces exposed by SOFA's Lagrangian solver."""

import numpy as np


def read_contact_force(root):
    """Return the contact's normal reaction in solver force units.

    The project places the contact object on the tube's +Y side, so the
    contact normal is parallel to Y. Read the mapped collision MechanicalObject
    lambda, following the installed SofaPython3 access_contact_forces example.
    """
    soft_node = root.getChild("SoftContactObject")
    rigid_node = root.getChild("RigidContactObject")
    contact_node = soft_node if soft_node is not None else rigid_node
    if contact_node is None:
        return 0.0

    collision_node = contact_node.getChild("CollisionModel")
    if collision_node is None:
        raise RuntimeError("CollisionModel del contatto non trovato")
    mechanical_object = collision_node.getObject("dofs")
    if mechanical_object is None:
        raise RuntimeError("MechanicalObject del contatto non trovato")

    values = np.asarray(
        mechanical_object.getData("lambda").value,
        dtype=float,
    )
    if values.size == 0:
        return 0.0
    if values.size % 3 != 0:
        raise RuntimeError(
            f"lambda del modello di contatto ha forma non valida: {values.shape}"
        )
    reaction_vectors = values.reshape((-1, 3))
    if not np.all(np.isfinite(reaction_vectors)):
        raise ValueError("lambda del contatto contiene valori non finiti")

    return float(abs(np.sum(reaction_vectors[:, 1])))

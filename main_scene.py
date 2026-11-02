# main_scene.py - SOFA v26.06 (SofaPython3, Python 3.12)
import os
import sys
import numpy as np
import Sofa
import Sofa.Core

# Aggiunge la directory corrente al path per importare il package 'modules'
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.config import NOMINAL_PARAMS
from modules.geometry_gen import (
    box_hexa_mesh,
    box_surface,
    create_hollow_tube_meshes,
    graded_axis,
    check_tet_mesh,
    graded_x_vals,
    write_surface_obj,
)
from modules.controller import TendonSpringController
from modules.contact_force import read_contact_force


_COLLISION_MODELS = {
    "triangle": "TriangleCollisionModel",
    "line": "LineCollisionModel",
    "point": "PointCollisionModel",
}


def _collision_model_names(names):
    """Converte ('Point', 'Triangle', ...) nei componenti SOFA corrispondenti."""
    if isinstance(names, str):
        names = (names,)
    try:
        models = [_COLLISION_MODELS[str(n).lower()] for n in names]
    except KeyError as exc:
        raise ValueError(f"Modello di collisione sconosciuto: {exc}; usare Triangle/Line/Point")
    if not models:
        raise ValueError("Serve almeno un modello di collisione")
    return models


def createScene(rootNode, custom_params=None):
    params = NOMINAL_PARAMS.copy()
    if custom_params is not None:
        params.update(custom_params)
    if params["contact_object"] not in {"soft", "rigid", "none"}:
        raise ValueError("contact_object deve essere 'soft', 'rigid' oppure 'none'")

    # 1. Configurazione Root Node (unita': mm, N, MPa, s)
    rootNode.dt = params["dt"]
    rootNode.gravity = [0.0, 0.0, 0.0]

    # 2. Plugin richiesti
    plugins = [
        "Sofa.Component.AnimationLoop",
        "Sofa.Component.Collision.Detection.Algorithm",
        "Sofa.Component.Collision.Detection.Intersection",
        "Sofa.Component.Collision.Geometry",
        "Sofa.Component.Collision.Response.Contact",
        "Sofa.Component.Constraint.Lagrangian.Correction",
        "Sofa.Component.Constraint.Lagrangian.Solver",
        "Sofa.Component.Constraint.Projective",
        "Sofa.Component.Engine.Select",
        "Sofa.Component.IO.Mesh",
        "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.Mapping.Linear",
        "Sofa.Component.Mass",
        "Sofa.Component.ODESolver.Backward",
        "Sofa.Component.Setting",
        "Sofa.Component.SolidMechanics.FEM.Elastic",
        "Sofa.Component.StateContainer",
        "Sofa.Component.Topology.Container.Constant",
        "Sofa.Component.Topology.Container.Dynamic",
        "Sofa.Component.Topology.Container.Grid",
        "Sofa.Component.Topology.Mapping",
        "Sofa.Component.Visual",
        "Sofa.GL.Component.Rendering3D",
        "SoftRobots",
    ]
    for p in plugins:
        rootNode.addObject("RequiredPlugin", name="req_" + p.replace(".", "_"), pluginName=p)

    # Camera e stile visivo
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject("VisualStyle", name="visualStyle",
                       displayFlags="showVisualModels showBehaviorModels showForceFields")

    # Loop di animazione e Constraint Solver
    rootNode.addObject("FreeMotionAnimationLoop", name="FreeMotionLoop")
    rootNode.addObject("BlockGaussSeidelConstraintSolver", name="ConstraintSolver",
                       maxIterations=500, tolerance=1e-7,
                       computeConstraintForces=True)

    # Pipeline di collisione e risposta con attrito mu = 0.2
    rootNode.addObject("CollisionPipeline", name="CollisionPipeline", depth=6)
    rootNode.addObject("BruteForceBroadPhase", name="BroadPhase")
    rootNode.addObject("BVHNarrowPhase", name="NarrowPhase")
    rootNode.addObject(
        "LocalMinDistance",
        name="LocalMinDistance",
        alarmDistance=params["alarm_distance"],
        contactDistance=params["contact_distance"],
        angleCone=0.01,
    )
    rootNode.addObject("CollisionResponse", name="Response",
                       response="FrictionContactConstraint",
                       responseParams="mu=" + str(params["mu_friction"]))

    # 3. Mesh del tubo: FEM graduata (fine sulla parte flessibile) e superfici separate,
    #    piu' grossolane, per i modelli visivo e di collisione (BarycentricMapping).
    mesh_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meshes")
    fem_obj_path = os.path.join(mesh_dir, "tube_fem_surface.obj")
    obj_path = os.path.join(mesh_dir, "tube_surface.obj")            # modello visivo
    collision_obj_path = os.path.join(mesh_dir, "tube_collision.obj")  # modello di collisione
    n_theta = int(params["mesh_n_theta"])
    fem_nodes, tetrahedra = create_hollow_tube_meshes(
        obj_filepath=fem_obj_path,
        od=params["od"],
        thickness=params["thickness"],
        n_theta=n_theta,
        x_vals=graded_x_vals(params["length_total"], params["length_shaft"],
                             params["mesh_h_flex"], int(params["mesh_n_shaft"])),
    )
    write_surface_obj(obj_path, params["od"], params["thickness"],
                      0.0, params["length_total"], params["visual_h"], n_theta)
    collision_x_start = max(
        params["length_shaft"],
        params["length_total"] - params["obj_size"] - params["collision_margin"],
    )
    write_surface_obj(collision_obj_path, params["od"], params["thickness"],
                      collision_x_start, params["length_total"], params["collision_h"], n_theta)
    bad, vmin, tube_volume = check_tet_mesh(fem_nodes, tetrahedra)
    if bad > 0:
        raise RuntimeError(f"Mesh tubo non valida: {bad} tetraedri invertiti/degeneri (vol min {vmin:.4g})")

    nodes_np = np.array(fem_nodes)
    tip_node_indices = np.where(np.isclose(nodes_np[:, 0], params["length_total"], atol=1e-3))[0].tolist()

    # 4. Nodo meccanico del tubo in Nitinol
    tubeNode = rootNode.addChild("NitinolTube")
    tubeNode.addObject("EulerImplicitSolver", name="odesolver", rayleighStiffness=0.1, rayleighMass=0.1)
    tubeNode.addObject("SparseLDLSolver", name="linearsolver", template="CompressedRowSparseMatrixMat3x3d")

    tubeNode.addObject("MeshTopology", name="topology", position=fem_nodes, tetrahedra=tetrahedra)
    tubeMeca = tubeNode.addObject("MechanicalObject", name="dofs", template="Vec3d")
    # Massa reale dalla densita' (mm-N-s -> tonnellate): ~7.6e-7 t = 0.76 g per il tubo baseline
    tubeNode.addObject("UniformMass", name="mass", totalMass=params["rho_tube"] * tube_volume)

    femFF = tubeNode.addObject(
        "TetrahedronFEMForceField",
        name="FEM",
        youngModulus=params["E_load"],
        poissonRatio=params["poisson"],
        method="large",
    )

    # Vincolo albero prossimale (primi 40 mm)
    tubeNode.addObject("BoxROI", name="shaftROI",
                       box=[-0.5, -2.0, -2.0, params["length_shaft"] + 0.01, 2.0, 2.0],
                       drawBoxes=True)
    tubeNode.addObject("FixedProjectiveConstraint", name="FixedShaft", indices="@shaftROI.indices")
    tubeNode.addObject("LinearSolverConstraintCorrection", name="LinearConstraintCorrection")

    # 5. Sotto-nodo Tendine (CableConstraint)
    tendonNode = tubeNode.addChild("Tendon")
    num_cable_pts = 15
    cable_x = np.linspace(params["length_shaft"], params["length_total"], num_cable_pts)
    cable_pts = [[float(x), float(params["yt_offset"]), 0.0] for x in cable_x]

    tendonNode.addObject("MechanicalObject", name="cablePts", position=cable_pts)
    cableConstraint = tendonNode.addObject(
        "CableConstraint",
        name="tendonCable",
        indices=list(range(num_cable_pts)),
        pullPoint=[params["length_shaft"], params["yt_offset"], 0.0],
        valueType="force",
        value=0.0,
        hasPullPoint=True,
    )
    tendonNode.addObject("BarycentricMapping", name="tendonMapping")

    # 6. Modello visivo
    visualNode = tubeNode.addChild("VisualModel")
    visualNode.addObject("MeshOBJLoader", name="objLoader", filename=obj_path)
    visualNode.addObject("OglModel", name="Visual", src="@objLoader", color=[0.75, 0.78, 0.82, 1.0])
    visualNode.addObject("BarycentricMapping", name="visualMapping")

    tubeCollisionNode = tubeNode.addChild("CollisionModel")
    tubeCollisionNode.addObject("MeshOBJLoader", name="collisionLoader", filename=collision_obj_path)
    tubeCollisionNode.addObject("MeshTopology", src="@collisionLoader")
    tubeCollisionNode.addObject("MechanicalObject", name="dofs", src="@collisionLoader", template="Vec3d")
    for model in _collision_model_names(params["tube_collision_models"]):
        tubeCollisionNode.addObject(model, group=0, simulated=True, moving=True)
    tubeCollisionNode.addObject("BarycentricMapping")

    # 7. Oggetto di contatto selezionabile.
    L = params["length_total"]
    contact_object = params["contact_object"]
    if contact_object != "none":
        cube_size = float(params["obj_size"])
        grid_n = int(params["obj_grid_n"])
        overlap_x = float(params["obj_overlap_x"])
        if cube_size <= 0 or overlap_x < 0 or grid_n < 2:
            raise ValueError("obj_size > 0, obj_overlap_x >= 0 e obj_grid_n >= 2 richiesti")

        x_min = L - cube_size
        x_max = L + overlap_x
        y_min = params["od"] / 2.0 + float(params["obj_gap"])
        y_max = y_min + cube_size
        z_min = -cube_size / 2.0
        z_max = cube_size / 2.0

        graded = contact_object == "soft" and params["soft_mesh"] == "graded"
        if graded:
            # Mesh graduata: passo soft_h_min sulla striscia di contatto (punta del tubo,
            # faccia verso il tubo, piano z = 0), crescente fino a soft_h_max.
            h_args = (params["soft_h_min"], params["soft_h_max"], params["soft_growth"])
            xs = graded_axis(x_min, x_max, L, *h_args)
            ys = graded_axis(y_min, y_max, y_min, *h_args)
            zs = graded_axis(z_min, z_max, 0.0, *h_args)
        else:
            xs = np.linspace(x_min, x_max, grid_n)
            ys = np.linspace(y_min, y_max, grid_n)
            zs = np.linspace(z_min, z_max, grid_n)
        collision_positions, collision_triangles = box_surface(xs, ys, zs)

        contactNode = rootNode.addChild(
            "SoftContactObject" if contact_object == "soft" else "RigidContactObject"
        )
        if contact_object == "soft":
            contactNode.addObject("EulerImplicitSolver", name="odesolver",
                                  rayleighStiffness=0.1, rayleighMass=0.1)
            contactNode.addObject("SparseLDLSolver", name="linearsolver",
                                  template="CompressedRowSparseMatrixMat3x3d")
            if graded:
                soft_pos, soft_hexa = box_hexa_mesh(xs, ys, zs)
                contactNode.addObject("HexahedronSetTopologyContainer", name="grid",
                                      position=soft_pos, hexahedra=soft_hexa)
                contactNode.addObject("MechanicalObject", name="state", template="Vec3d",
                                      position=soft_pos)
            else:
                contactNode.addObject(
                    "RegularGridTopology",
                    name="grid",
                    n=[grid_n, grid_n, grid_n],
                    min=[x_min, y_min, z_min],
                    max=[x_max, y_max, z_max],
                )
                contactNode.addObject("MechanicalObject", name="state", template="Vec3d")
            soft_volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
            contactNode.addObject("UniformMass", name="mass",
                                  totalMass=params["rho_soft"] * soft_volume)
            contactNode.addObject(
                "HexahedronFEMForceField",
                name="SoftFEM",
                topology="@grid",
                youngModulus=params["E_soft"],
                poissonRatio=params["poisson_soft"],
                method="large",
            )
            contactNode.addObject(
                "BoxROI",
                name="backROI",
                box=[
                    x_min - 0.1, y_max - 0.1, z_min - 0.1,
                    x_max + 0.1, y_max + 0.1, z_max + 0.1,
                ],
            )
            contactNode.addObject(
                "FixedProjectiveConstraint",
                name="FixedBack",
                indices="@backROI.indices",
            )
            contactNode.addObject(
                "LinearSolverConstraintCorrection",
                name="LinearConstraintCorrection",
            )

        collisionNode = contactNode.addChild("CollisionModel")
        collisionNode.addObject(
            "TriangleSetTopologyContainer",
            name="topology",
            position=collision_positions,
            triangles=collision_triangles,
        )
        collisionNode.addObject("TriangleSetTopologyModifier", name="modifier")
        collisionNode.addObject(
            "TriangleSetGeometryAlgorithms",
            name="geometry",
            template="Vec3d",
        )
        collisionNode.addObject(
            "MechanicalObject",
            name="dofs",
            template="Vec3d",
            position=collision_positions,
        )
        for model in _collision_model_names(params["object_collision_models"]):
            collisionNode.addObject(
                model,
                group=1,
                simulated=True,
                moving=(contact_object == "soft"),
            )
        if contact_object == "soft":
            collisionNode.addObject("BarycentricMapping")
        else:
            collisionNode.addObject(
                "FixedProjectiveConstraint",
                name="FixedRigidSurface",
                indices=list(range(len(collision_positions))),
            )

        visualNode = contactNode.addChild("VisualModel")
        visualNode.addObject(
            "OglModel",
            name="Visual",
            src="@../CollisionModel/topology",
            color=[0.85, 0.35, 0.25, 1.0],
        )
        if contact_object == "soft":
            visualNode.addObject("BarycentricMapping")

    # 8. Controller Python
    rootNode.addObject(
        TendonSpringController(
            name="ActuationController",
            root=rootNode,
            tube_meca=tubeMeca,
            fem_ff=femFF,
            cable=cableConstraint,
            params=params,
            tip_node_indices=tip_node_indices,
            contact_force_provider=params.get(
                "contact_force_provider",
                read_contact_force,
            ),
        )
    )

    return rootNode
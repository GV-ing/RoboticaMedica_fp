# main_scene.py - SOFA v26.06 (SofaPython3, Python 3.12)
import os
import sys
import numpy as np
import Sofa
import Sofa.Core

# Aggiunge la directory corrente al path per importare il package 'modules'
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.config import NOMINAL_PARAMS
from modules.geometry_gen import create_hollow_tube_meshes, check_tet_mesh
from modules.controller import TendonSpringController
from modules.contact_force import read_contact_force


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

    # 3. Mesh tetraedrica del tubo
    obj_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meshes", "tube_surface.obj")
    fem_nodes, tetrahedra = create_hollow_tube_meshes(
        obj_filepath=obj_path,
        od=params["od"],
        thickness=params["thickness"],
        length=params["length_total"],
    )
    bad, vmin, _ = check_tet_mesh(fem_nodes, tetrahedra)
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
    tubeNode.addObject("UniformMass", name="mass", totalMass=1e-3)

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
    tubeCollisionNode.addObject("MeshOBJLoader", name="collisionLoader", filename=obj_path)
    tubeCollisionNode.addObject("MeshTopology", src="@collisionLoader")
    tubeCollisionNode.addObject("MechanicalObject", name="dofs", src="@collisionLoader", template="Vec3d")
    tubeCollisionNode.addObject("TriangleCollisionModel", group=0, simulated=True, moving=True)
    tubeCollisionNode.addObject("LineCollisionModel", group=0, simulated=True, moving=True)
    tubeCollisionNode.addObject("PointCollisionModel", group=0, simulated=True, moving=True)
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
        y_min = params["od"] / 2.0 + params["obj_gap"]
        y_max = y_min + cube_size
        z_min = -cube_size / 2.0
        z_max = cube_size / 2.0

        x_coords = np.linspace(x_min, x_max, grid_n)
        y_coords = np.linspace(y_min, y_max, grid_n)
        z_coords = np.linspace(z_min, z_max, grid_n)
        boundary_indices = {}
        collision_positions = []
        for i in range(grid_n):
            for j in range(grid_n):
                for k in range(grid_n):
                    if i in (0, grid_n - 1) or j in (0, grid_n - 1) or k in (0, grid_n - 1):
                        boundary_indices[(i, j, k)] = len(collision_positions)
                        collision_positions.append([
                            float(x_coords[i]),
                            float(y_coords[j]),
                            float(z_coords[k]),
                        ])

        collision_triangles = []

        def add_collision_quad(corners):
            a, b, c, d = (boundary_indices[corner] for corner in corners)
            collision_triangles.extend([[a, b, c], [a, c, d]])

        for j in range(grid_n - 1):
            for k in range(grid_n - 1):
                add_collision_quad([(0, j, k), (0, j, k + 1),
                                    (0, j + 1, k + 1), (0, j + 1, k)])
                add_collision_quad([(grid_n - 1, j, k), (grid_n - 1, j + 1, k),
                                    (grid_n - 1, j + 1, k + 1), (grid_n - 1, j, k + 1)])
        for i in range(grid_n - 1):
            for k in range(grid_n - 1):
                add_collision_quad([(i, 0, k), (i + 1, 0, k),
                                    (i + 1, 0, k + 1), (i, 0, k + 1)])
                add_collision_quad([(i, grid_n - 1, k), (i, grid_n - 1, k + 1),
                                    (i + 1, grid_n - 1, k + 1), (i + 1, grid_n - 1, k)])
        for i in range(grid_n - 1):
            for j in range(grid_n - 1):
                add_collision_quad([(i, j, 0), (i, j + 1, 0),
                                    (i + 1, j + 1, 0), (i + 1, j, 0)])
                add_collision_quad([(i, j, grid_n - 1), (i + 1, j, grid_n - 1),
                                    (i + 1, j + 1, grid_n - 1), (i, j + 1, grid_n - 1)])

        contactNode = rootNode.addChild(
            "SoftContactObject" if contact_object == "soft" else "RigidContactObject"
        )
        if contact_object == "soft":
            contactNode.addObject("EulerImplicitSolver", name="odesolver",
                                  rayleighStiffness=0.1, rayleighMass=0.1)
            contactNode.addObject("SparseLDLSolver", name="linearsolver",
                                  template="CompressedRowSparseMatrixMat3x3d")
            contactNode.addObject(
                "RegularGridTopology",
                name="grid",
                n=[grid_n, grid_n, grid_n],
                min=[x_min, y_min, z_min],
                max=[x_max, y_max, z_max],
            )
            contactNode.addObject("MechanicalObject", name="state", template="Vec3d")
            contactNode.addObject("UniformMass", name="mass", totalMass=0.01)
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
        collisionNode.addObject(
            "TriangleCollisionModel",
            group=1,
            simulated=True,
            moving=(contact_object == "soft"),
        )
        collisionNode.addObject(
            "LineCollisionModel",
            group=1,
            simulated=True,
            moving=(contact_object == "soft"),
        )
        collisionNode.addObject(
            "PointCollisionModel",
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
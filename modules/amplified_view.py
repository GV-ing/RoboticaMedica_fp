# modules/amplified_view.py
"""Vista amplificata della deformazione, SOLO grafica (GUI).

Le deformazioni reali sono di decimi di millimetro (punta: ~0.3 mm a T = 10 N; indentazione
del soffice: ~0.1-0.3 mm) e sullo schermo non si vedono. Come nei post-processori FEM, un
secondo modello visivo mostra

    posizione_mostrata = riposo + A * (posizione_reale - riposo)

La meccanica non cambia: il modello visivo originale (con BarycentricMapping, come richiesto
dalla traccia) resta nella scena, reso semitrasparente come "ombra" della forma reale, e
fornisce le posizioni reali da amplificare.
"""

import numpy as np
import Sofa
import Sofa.Core


class AmplifiedView(Sofa.Core.Controller):
    """Copia ogni passo le posizioni del modello visivo reale, amplificate di `factor`."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pairs = kwargs.get("pairs")      # lista di (OglModel reale, OglModel amplificato)
        self.factor = float(kwargs.get("factor", 50.0))
        self.rest = [None] * len(self.pairs)
        self.enabled = True

    def onAnimateEndEvent(self, event):
        if not self.enabled:
            return
        for i, (real, shown) in enumerate(self.pairs):
            pos = np.asarray(real.position.value, dtype=float)
            if self.rest[i] is None:
                self.rest[i] = pos.copy()
                if len(shown.position.value) != len(pos):
                    print("[VISTA AMPLIFICATA] numero di vertici diverso: vista disattivata")
                    self.enabled = False
                    return
            with shown.position.writeableArray() as out:
                out[:] = self.rest[i] + self.factor * (pos - self.rest[i])


def add_amplified_view(root, factor, tube_obj_path, soft_surface=None):
    """Aggiunge i modelli amplificati per tubo e (se presente) oggetto soffice."""
    # Nota: in OglModel "color" vale solo alla creazione (non e' un Data modificabile),
    # quindi la trasparenza del modello reale si imposta in main_scene.py (ghost_alpha).
    pairs = []
    tube_real = root.getChild("NitinolTube").getChild("VisualModel").getObject("Visual")
    node = root.addChild("TubeAmplifiedView")
    node.addObject("MeshOBJLoader", name="loader", filename=tube_obj_path)
    tube_shown = node.addObject("OglModel", name="Visual", src="@loader",
                                color=[0.16, 0.47, 0.84, 1.0])
    pairs.append((tube_real, tube_shown))

    soft = root.getChild("SoftContactObject")
    if soft is not None and soft_surface is not None:
        soft_real = soft.getChild("VisualModel").getObject("Visual")
        positions, triangles = soft_surface
        node = root.addChild("SoftAmplifiedView")
        soft_shown = node.addObject("OglModel", name="Visual", position=positions,
                                    triangles=triangles, color=[0.92, 0.41, 0.20, 1.0])
        pairs.append((soft_real, soft_shown))

    root.addObject(AmplifiedView(name="AmplifiedView", pairs=pairs, factor=factor))
    print(f"[VISTA AMPLIFICATA] deformazioni mostrate x{factor:g} (solo grafica; "
          f"forma reale in trasparenza)")


# modules/controller.py (Aggiornato con logica a due livelli)
import Sofa
import Sofa.Core
import numpy as np
import logging

class TendonSpringController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.root = kwargs.get("root")
        self.tube_meca = kwargs.get("tube_meca")
        self.fem_ff = kwargs.get("fem_ff")
        self.cable = kwargs.get("cable")
        self.params = kwargs.get("params")
        self.tip_node_indices = kwargs.get("tip_node_indices")
        self.contact_force_provider = kwargs.get("contact_force_provider")
        self.force_scale = float(self.params["force_scale"])
        if not np.isfinite(self.force_scale) or self.force_scale <= 0.0:
            raise ValueError("force_scale deve essere un fattore positivo e finito")
        
        self.u = 0.0
        self.u_target = 0.0
        self.delta_Ls = 0.0
        self.prev_delta_Ls = 0.0
        self.T = 0.0
        self.direction = 1.0
        self.auto_direction = 1.0
        self.autonomous_mode = False
        
        self.contact_force_N = None
        self._warned_missing_cable_displacement = False

        self.history_u = []
        self.history_tip_y = []
        self.history_tendon_force_N = []
        self.history_contact_force_N = []
        self.initial_tip_y = None

    def get_tip_position(self):
        pos = self.tube_meca.position.value
        tip_pts = pos[self.tip_node_indices]
        return np.mean(tip_pts, axis=0)

    def get_tendon_displacement(self):
        tip_pos = self.get_tip_position()
        try:
            cable_disp = float(self.cable.displacement.value[0])
            return max(0.0, cable_disp)
        except (AttributeError, IndexError, TypeError, ValueError):
            if not self._warned_missing_cable_displacement:
                logging.warning(
                    "CableConstraint.displacement non disponibile; "
                    "uso lo spostamento geometrico approssimato."
                )
                self._warned_missing_cable_displacement = True
            delta_y = max(0.0, tip_pos[1])
            return (delta_y / self.params["length_flex"]) * self.params["yt_offset"]

    def read_contact_force(self):
        if self.contact_force_provider is None:
            return None
        force = float(self.contact_force_provider(self.root))
        if not np.isfinite(force) or force < 0.0:
            raise ValueError(f"Forza di contatto non valida: {force}")
        return force * self.force_scale

    def update_actuation_physics(self, dt):
        target_E = self.params["E_load"] if self.direction >= 0 else self.params["E_unload"]
        if abs(float(self.fem_ff.youngModulus.value[0]) - target_E) > 1.0:
            self.fem_ff.youngModulus.value = [target_E]

        delta_Lt = self.get_tendon_displacement()
        self.delta_Ls = max(0.0, self.u - delta_Lt)
        d_delta_Ls = (self.delta_Ls - self.prev_delta_Ls) / dt if dt > 0 else 0.0
        self.prev_delta_Ls = self.delta_Ls

        # Selezione rigidezza molla: Config (a) fissa vs Config (b) a due livelli
        # Se k_mode e' specificato nei parametri usiamo la logica switched, altrimenti k_spring fissa
        if self.params.get("switched_mode", False):
            if self.contact_force_provider is None:
                raise RuntimeError(
                    "La configurazione switched richiede contact_force_provider(root), "
                    "il cui valore viene convertito in N usando force_scale."
                )
            if (self.contact_force_N or 0.0) >= self.params["F_contact_th"]:
                k = self.params["k_high"] # 4.8 N/mm in contatto
            else:
                k = self.params["k_low"]  # 0.6 N/mm in moto libero[cite: 10]
        else:
            k = self.params["k_spring"]   # 1.2 N/mm fissa[cite: 10]

        c = self.params["c_spring"]
        
        if self.u <= 1e-5:
            self.T = 0.0
            self.delta_Ls = 0.0
        else:
            hysteresis_term = self.direction * (self.params["Fh_spring"] + self.params["Fh_mat"])
            raw_T = k * self.delta_Ls + c * d_delta_Ls + hysteresis_term
            self.T = float(np.clip(raw_T, 0.0, 15.0))

        self.cable.value.value = [self.T]
        
        tip = self.get_tip_position()
        if self.initial_tip_y is None:
            self.initial_tip_y = float(tip[1])
        self.history_u.append(self.u)
        self.history_tip_y.append(float(tip[1]) - self.initial_tip_y)
        self.history_tendon_force_N.append(self.T)

    def onAnimateEndEvent(self, event):
        self.contact_force_N = self.read_contact_force()
        self.history_contact_force_N.append(self.contact_force_N)

    def onAnimateBeginEvent(self, event):
        dt = self.root.dt.value
        if self.autonomous_mode:
            self.u_target += self.auto_direction * self.params["auto_rate"] * dt
            if self.u_target >= self.params["u_max"]:
                self.u_target = self.params["u_max"]
                self.auto_direction = -1.0
            elif self.u_target <= 0.0:
                self.u_target = 0.0
                self.auto_direction = 1.0

        rate = self.params["auto_rate"] if self.autonomous_mode else self.params["manual_rate"]
        max_step = rate * dt
        error = self.u_target - self.u
        previous_u = self.u
        self.u += float(np.clip(error, -max_step, max_step))
        if self.u > previous_u:
            self.direction = 1.0
        elif self.u < previous_u:
            self.direction = -1.0
        self.update_actuation_physics(dt)

    def onKeypressedEvent(self, event):
        key = event.get("key")
        step = self.params["u_step"]
        if key in ["\x13", "U", "u", "+", "\x15", "D", "d", "-", "R", "r"]:
            if self.autonomous_mode:
                self.autonomous_mode = False
                self.u_target = self.u

        if key in ["\x13", "U", "u", "+"]:
            self.u_target = min(self.params["u_max"], self.u_target + step)
        elif key in ["\x15", "D", "d", "-"]:
            self.u_target = max(0.0, self.u_target - step)
        elif key in ["R", "r"]:
            self.u_target = 0.0
        elif key in ["M", "m"]:
            self.autonomous_mode = not self.autonomous_mode
            self.u_target = self.u
            self.auto_direction = self.direction
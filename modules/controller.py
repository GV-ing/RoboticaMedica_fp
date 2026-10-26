# modules/controller.py
import logging
import numpy as np
import Sofa
import Sofa.Core


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

    # Variabili di stato cinematico e di attuazione
    self.u = 0.0  # Spostamento carrello motore [mm]
    self.u_target = 0.0  # Riferimento filtrato dal rate limiter [mm]
    self.delta_Ls = 0.0  # Allungamento istantaneo molla SEA [mm]
    self.prev_delta_Ls = 0.0
    self.T = 0.0  # Tensione al cavo [N]
    self.direction = 1.0  # Direzione reale di moto (+1 carico, -1 scarico)
    self.auto_direction = 1.0
    self.autonomous_mode = False

    # Variabili per la prova di stabilità di contatto (Target: 20 mN per 5 s)
    self.contact_force_N = 0.0
    self.contact_detected = False
    self.contact_hold_timer = 0.0
    self.hold_duration = 5.0  # Secondi di mantenimento richiesti

    self._warned_missing_cable_displacement = False

    # Cronologia e tracciamento
    self.history_u = []
    self.history_tip_y = []
    self.history_tendon_force_N = []
    self.history_contact_force_N = []
    self.initial_tip_y = None

  def get_tip_position(self):
    """Calcola la posizione cartesiana media dei nodi di punta [mm]."""
    pos = self.tube_meca.position.value
    tip_pts = pos[self.tip_node_indices]
    return np.mean(tip_pts, axis=0)

  def get_tendon_displacement(self):
    """Restituisce lo spostamento all'ancoraggio del cavo (u1) [mm]."""
    tip_pos = self.get_tip_position()
    try:
      cable_disp = float(self.cable.displacement.value[0])
      return max(0.0, cable_disp)
    except (AttributeError, IndexError, TypeError, ValueError):
      if not self._warned_missing_cable_displacement:
        logging.warning(
            "CableConstraint.displacement non accessibile; uso approssimazione"
            " geometrica."
        )
        self._warned_missing_cable_displacement = True
      delta_y = max(0.0, tip_pos[1])
      return (delta_y / self.params["length_flex"]) * self.params["yt_offset"]

  def read_contact_force(self):
    """Legge la forza di contatto normale convertita in Newton [N]."""
    if self.contact_force_provider is None:
      return 0.0
    try:
      force = float(self.contact_force_provider(self.root))
      if not np.isfinite(force) or force < 0.0:
        return 0.0
      return force * self.force_scale
    except Exception:
      return 0.0

  def print_state_feedback(self):
    """Stampa a console le variabili richieste dal capitolato a ogni comando."""
    tip = self.get_tip_position()
    delta_y = (
        float(tip[1]) - self.initial_tip_y
        if self.initial_tip_y is not None
        else 0.0
    )
    print(
        f"[CONTROLLER] u: {self.u:.3f} mm | T: {self.T:.3f} N | "
        f"ΔLs: {self.delta_Ls:.3f} mm | Punta: X={tip[0]:.2f}, Y={tip[1]:.2f}"
        f" (ΔY={delta_y:.2f}), Z={tip[2]:.2f} mm | F_contatto:"
        f" {self.contact_force_N * 1e3:.2f} mN"
    )

  def update_actuation_physics(self, dt):
    """Aggiorna le equazioni costitutive della molla in serie, isteresi e attuazione."""
    # 1. Aggiornamento Young Modulus asimmetrico (Nitinol: E_load vs E_unload)
    target_E = (
        self.params["E_load"]
        if self.direction >= 0
        else self.params["E_unload"]
    )
    current_E = float(self.fem_ff.youngModulus.value[0])
    if abs(current_E - target_E) > 10.0:
      self.fem_ff.youngModulus.value = [target_E]

    # 2. Cinematica SEA: allungamento molla Delta L_s
    delta_Lt = self.get_tendon_displacement()
    self.delta_Ls = max(0.0, self.u - delta_Lt)
    d_delta_Ls = (self.delta_Ls - self.prev_delta_Ls) / dt if dt > 0.0 else 0.0
    self.prev_delta_Ls = self.delta_Ls

    # 3. Logica di Rigidezza Molla: Configurazione (a) Fissa vs (b) Switched a due livelli
    if self.params.get("switched_mode", False):
      # Soglia di commutazione di specifica: 5 mN = 0.005 N
      if self.contact_force_N >= self.params.get("F_contact_th", 0.005):
        k = self.params["k_high"]  # 4.8 N/mm a contatto avvenuto
      else:
        k = self.params["k_low"]  # 0.6 N/mm in moto libero
    else:
      k = self.params["k_spring"]  # 1.2 N/mm costante

    c = self.params["c_spring"]

    # 4. Calcolo Tensione Tendine con regolarizzazione continua dell'Isteresi
    if self.u <= 1e-4:
      self.T = 0.0
      self.delta_Ls = 0.0
    else:
      # Regolarizzazione tangenziale per prevenire gradini istantanei di forza (anti-scatto)
      # sigma definisce la pendenza di transizione della frizione/isteresi interna
      hyst_magnitude = self.params["Fh_spring"] + self.params["Fh_mat"]
      hyst_smooth = hyst_magnitude * np.tanh(10.0 * self.direction)

      raw_T = k * self.delta_Ls + c * d_delta_Ls + hyst_smooth
      # Saturazione a T_max = 10.0 N come prescritto dalle specifiche generali
      self.T = float(np.clip(raw_T, 0.0, self.params.get("T_max", 10.0)))

    # Assegnazione del carico al vincolo di fune
    self.cable.value.value = [self.T]

    # Tracciamento telemetria
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

    # Gestione Ciclo Autonomo con Protocollo di Mantenimento Contatto
    if self.autonomous_mode:
      f_target = 0.020  # Setpoint di contatto stabile: 20 mN = 0.020 N
      f_thresh = self.params.get("F_contact_th", 0.005)  # 5 mN

      # Rilevamento contatto durante la fase di carico
      if self.auto_direction > 0 and self.contact_force_N >= f_thresh:
        self.contact_detected = True

      if self.contact_detected and self.contact_hold_timer < self.hold_duration:
        # Micro-regolazione per mantenere la forza a 20 mN evitando penetrazione e sovra-forza
        force_err = f_target - self.contact_force_N
        # Piccola correzione proporzionale sulla posizione
        self.u_target += 0.05 * force_err * dt
        self.contact_hold_timer += dt
      else:
        # Procedura ordinaria di ciclo carico/scarico
        self.u_target += self.auto_direction * self.params["auto_rate"] * dt

        if self.u_target >= self.params["u_max"]:
          self.u_target = self.params["u_max"]
          self.auto_direction = -1.0
          self.contact_detected = False
          self.contact_hold_timer = 0.0
        elif self.u_target <= 0.0:
          self.u_target = 0.0
          self.auto_direction = 1.0
          self.contact_detected = False
          self.contact_hold_timer = 0.0

    # Rate Limiter per evitare comandi a gradino rigido
    rate = (
        self.params["auto_rate"]
        if self.autonomous_mode
        else self.params["manual_rate"]
    )
    max_step = rate * dt
    error = self.u_target - self.u
    previous_u = self.u
    self.u += float(np.clip(error, -max_step, max_step))

    # Calcolo direzione per modello isteresi
    if self.u > previous_u + 1e-6:
      self.direction = 1.0
    elif self.u < previous_u - 1e-6:
      self.direction = -1.0

    self.update_actuation_physics(dt)

  def onKeypressedEvent(self, event):
    key = str(event.get("key", ""))
    step = self.params["u_step"]  # 0.35 mm da specifica baseline

    # Codici SOFA standard e tasti alfanumerici
    # Key 19: Freccia SU | Key 21: Freccia GIÙ
    is_up = key in ["19", "Up", "U", "u", "+"]
    is_down = key in ["21", "Down", "D", "d", "-"]
    is_reset = key in ["R", "r"]
    is_auto = key in ["M", "m", "A", "a"]

    if is_up or is_down or is_reset:
      if self.autonomous_mode:
        self.autonomous_mode = False
        self.u_target = self.u

    if is_up:
      self.u_target = min(self.params["u_max"], self.u_target + step)
      self.print_state_feedback()
    elif is_down:
      self.u_target = max(0.0, self.u_target - step)
      self.print_state_feedback()
    elif is_reset:
      self.u_target = 0.0
      self.u = 0.0
      self.contact_detected = False
      self.contact_hold_timer = 0.0
      self.cable.value.value = [0.0]
      print("[CONTROLLER] Reset eseguito: u = 0.0 mm")
      self.print_state_feedback()
    elif is_auto:
      self.autonomous_mode = not self.autonomous_mode
      self.u_target = self.u
      self.auto_direction = 1.0
      self.contact_detected = False
      self.contact_hold_timer = 0.0
      state = "ATTIVATA" if self.autonomous_mode else "DISATTIVATA"
      print(f"[CONTROLLER] Modalità autonoma: {state}")
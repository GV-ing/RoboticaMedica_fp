# modules/controller.py
"""Controller da tastiera e autonomo del tubo a tendine con molla in serie (SofaPython3).

Comandi (traccia: carico/scarico, reset, manuale/autonomo; a ogni comando stampa u, T, dLs, punta):
    Freccia SU / U / +   carico   (+u_step)          Freccia GIU' / D / -   scarico (-u_step)
    R                    reset                        M / A                  manuale <-> autonomo

Ciclo autonomo (macchina a stati): 'load' a velocita' costante finche' T raggiunge T_max
(o u = u_max) -> 'hold' per auto_hold_s secondi -> 'unload' fino a u = 0 -> 'done'
(oppure di nuovo 'load' se auto_repeat). La fase e' registrata per separare carico e scarico.

Convenzioni SOFA 26.06 verificate (tools/): al CableConstraint si impone l'impulso T*dt;
cambiare youngModulus richiede reinit() del TetrahedronFEMForceField.
"""

import csv
import logging
import os

import numpy as np
import Sofa
import Sofa.Core

from modules.spring_model import SeriesSpring

HISTORY_FIELDS = ("t", "phase", "u", "u_target", "u1", "dLs", "T", "k", "z",
                  "E", "tip_x", "tip_y", "tip_z", "dY", "F_contact")


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

    self.spring = SeriesSpring(self.params)
    self.u = 0.0
    self.u_target = 0.0
    self.u1 = 0.0
    self.T = 0.0
    self.delta_Ls = 0.0
    self.E_current = float(self.params["E_load"])
    self.contact_force_N = 0.0
    self.contact_read_errors = 0
    self.time = 0.0

    self.autonomous_mode = False
    self.phase = "manual"
    self.hold_timer = 0.0

    self.initial_tip = None
    self.history = {name: [] for name in HISTORY_FIELDS}
    self.listeners = []   # callable(sample: dict) chiamati a ogni passo (grafici live)
    self._warned_missing_cable_displacement = False

  # ------------------------------------------------------------------ misure
  def get_tip_position(self):
    """Posizione media dei nodi di punta [mm]."""
    pos = self.tube_meca.position.value
    return np.mean(pos[self.tip_node_indices], axis=0)

  def get_tendon_displacement(self):
    """Spostamento del tendine a valle della molla u1 [mm] (CableConstraint.displacement)."""
    try:
      # In SoftRobots 26.06 'displacement' e' uno scalare: np.ravel gestisce scalare e vettore
      return max(0.0, float(np.ravel(self.cable.displacement.value)[0]))
    except (AttributeError, IndexError, TypeError, ValueError):
      if not self._warned_missing_cable_displacement:
        logging.warning("CableConstraint.displacement non accessibile; uso u1 = 0.")
        self._warned_missing_cable_displacement = True
      return 0.0

  def read_contact_force(self):
    """Forza normale di contatto [N] dal provider (modules/contact_force.py)."""
    if self.contact_force_provider is None:
      return 0.0
    try:
      force = float(self.contact_force_provider(self.root))
    except Exception as exc:  # errore di lettura: segnalato, non confuso con "nessun contatto"
      self.contact_read_errors += 1
      if self.contact_read_errors == 1:
        logging.error("Lettura forza di contatto fallita (%s); uso 0 N. "
                      "Errori contati in contact_read_errors.", exc)
      return 0.0
    if not np.isfinite(force) or force < 0.0:
      raise ValueError(f"Forza di contatto non valida: {force}")
    return force * self.force_scale

  # ------------------------------------------------------------------ stampa
  def print_state_feedback(self):
    tip = self.get_tip_position()
    dY = float(tip[1]) - float(self.initial_tip[1]) if self.initial_tip is not None else 0.0
    print(
        f"[CONTROLLER] u: {self.u:.3f} mm (rif {self.u_target:.3f}) | T: {self.T:.3f} N | "
        f"dLs: {self.delta_Ls:.3f} mm | k: {self.spring.k:.2f} N/mm | "
        f"Punta: X={tip[0]:.3f}, Y={tip[1]:.4f} (dY={dY:.4f}), Z={tip[2]:.4f} mm | "
        f"F_contatto: {self.contact_force_N * 1e3:.2f} mN"
    )

  # ------------------------------------------------------------------ fisica
  def _set_young_modulus(self, E):
    if abs(E - self.E_current) > 1e-9:
      self.fem_ff.youngModulus.value = [E]
      # TetrahedronFEMForceField NON ricalcola la rigidezza al solo cambio del Data:
      # serve reinit() (verificato con tools/check_scene.py: rapporto 1.000 -> 1.185).
      self.fem_ff.reinit()
      self.E_current = E

  def update_actuation_physics(self, dt):
    """Molla in serie -> tensione T sul cavo; E di carico o scarico secondo il verso di u."""
    self.u1 = self.get_tendon_displacement()
    self.T = self.spring.step(self.u, self.u1, self.contact_force_N, dt)
    self.delta_Ls = self.spring.dL_s
    self._set_young_modulus(
        self.params["E_load"] if self.spring.direction >= 0 else self.params["E_unload"])
    # Con FreeMotionAnimationLoop il valore del vincolo e' un IMPULSO: si impone T*dt.
    self.cable.value.value = [self.T * dt]

  def _autonomous_step(self, dt):
    rate = self.params["auto_rate"]
    if self.phase == "load":
      self.u_target = min(self.params["u_max"], self.u_target + rate * dt)
      if self.T >= self.spring.T_max * (1.0 - 1e-6) or self.u >= self.params["u_max"] - 1e-9:
        self.u_target = self.u
        self.phase = "hold" if self.params["auto_hold_s"] > 0 else "unload"
        self.hold_timer = 0.0
    elif self.phase == "hold":
      self.hold_timer += dt
      if self.hold_timer >= self.params["auto_hold_s"] - 1e-9:
        self.phase = "unload"
    elif self.phase == "unload":
      self.u_target = max(0.0, self.u_target - rate * dt)
      if self.u <= 1e-9 and self.u_target <= 0.0:
        self.phase = "load" if self.params["auto_repeat"] else "done"

  def onAnimateBeginEvent(self, event):
    dt = float(self.root.dt.value)
    if self.initial_tip is None:
      self.initial_tip = np.array(self.get_tip_position(), dtype=float)
    if self.autonomous_mode:
      self._autonomous_step(dt)
      rate = self.params["auto_rate"]
    else:
      rate = self.params["manual_rate"]
    # Velocita' del carrello limitata anche da |dT/dt| <= max_force_rate (glossario)
    max_step = self.spring.max_motor_rate(rate) * dt
    self.u += float(np.clip(self.u_target - self.u, -max_step, max_step))
    self.update_actuation_physics(dt)

  def onAnimateEndEvent(self, event):
    dt = float(self.root.dt.value)
    self.time += dt
    self.contact_force_N = self.read_contact_force()
    tip = self.get_tip_position()
    sample = {
        "t": self.time, "phase": self.phase, "u": self.u, "u_target": self.u_target,
        "u1": self.u1, "dLs": self.delta_Ls, "T": self.T, "k": self.spring.k,
        "z": self.spring.z, "E": self.E_current,
        "tip_x": float(tip[0]), "tip_y": float(tip[1]), "tip_z": float(tip[2]),
        "dY": float(tip[1] - self.initial_tip[1]), "F_contact": self.contact_force_N,
    }
    for name in HISTORY_FIELDS:
      self.history[name].append(sample[name])
    for listener in self.listeners:
      try:
        listener(sample)
      except Exception as exc:  # un grafico che fallisce non deve fermare la fisica
        logging.warning("Listener disattivato dopo errore: %s", exc)
        self.listeners.remove(listener)

  # ------------------------------------------------------------------ comandi
  def start_autonomous_cycle(self, hold_s=None, repeat=None):
    if hold_s is not None:
      self.params["auto_hold_s"] = float(hold_s)
    if repeat is not None:
      self.params["auto_repeat"] = bool(repeat)
    self.autonomous_mode = True
    self.phase = "load"
    self.u_target = self.u

  def reset(self):
    self.autonomous_mode = False
    self.phase = "manual"
    self.u = self.u_target = 0.0
    self.spring.reset()
    self.T = 0.0
    self.cable.value.value = [0.0]
    self._set_young_modulus(self.params["E_load"])

  def export_csv(self, path):
    """Scrive la cronologia in CSV (una riga per passo)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
      w = csv.writer(f)
      w.writerow(HISTORY_FIELDS)
      w.writerows(zip(*(self.history[n] for n in HISTORY_FIELDS)))
    return path

  def onKeypressedEvent(self, event):
    key = str(event.get("key", ""))
    step = self.params["u_step"]
    is_up = key in ["19", "Up", "U", "u", "+"]
    is_down = key in ["21", "Down", "D", "d", "-"]
    is_reset = key in ["R", "r"]
    is_auto = key in ["M", "m", "A", "a"]

    if (is_up or is_down) and self.autonomous_mode:
      self.autonomous_mode = False
      self.phase = "manual"
      self.u_target = self.u
    if is_up:
      self.u_target = min(self.params["u_max"], self.u_target + step)
      self.print_state_feedback()
    elif is_down:
      self.u_target = max(0.0, self.u_target - step)
      self.print_state_feedback()
    elif is_reset:
      self.reset()
      print("[CONTROLLER] Reset eseguito: u = 0.0 mm")
      self.print_state_feedback()
    elif is_auto:
      if self.autonomous_mode:
        self.autonomous_mode = False
        self.phase = "manual"
        self.u_target = self.u
      else:
        self.start_autonomous_cycle()
      state = "ATTIVATA" if self.autonomous_mode else "DISATTIVATA"
      print(f"[CONTROLLER] Modalita' autonoma: {state}")
      self.print_state_feedback()

  # ------------------------------------------------------------------ compatibilita'
  @property
  def history_tip_y(self):
    return self.history["dY"]

  @property
  def history_tendon_force_N(self):
    return self.history["T"]

  @property
  def history_contact_force_N(self):
    return self.history["F_contact"]

# Progetto 5 — Incertezza parametrica: Monte Carlo e analisi di sensibilità

Robotica Medica 2026/27 (prof.ssa F. Ficuciello). Tubo in nitinol azionato da un tendine con
molla in serie, simulato in **SOFA v26.06 (SofaPython3, Python 3.12) su Ubuntu 22.04**.

Si confrontano due configurazioni della molla: **(a)** fissa con k = 1.2 N/mm e **(b)** a due
livelli, **0.6 N/mm in moto libero e 4.8 N/mm in contatto**, con commutazione istantanea a 5 mN
(`k_free_level = "low"` in `modules/config.py`). Per ognuna si simulano il caso nominale e 20
realizzazioni con i parametri perturbati del ±(10–20) % con segno casuale. Su ogni realizzazione
si calcolano l'isteresi H e la stabilità del contatto; una regressione lineare dà il diagramma
tornado.

## Pipeline completa

| Passo | Comando | Durata indicativa | Produce |
|---|---|---|---|
| 0. Ambiente | `source env.sh` | – | variabili di SOFA per `python3.12` |
| 1. Mesh | `python3.12 tools/generate_meshes.py` | secondi | `meshes/*.obj`, `meshes/*.vtk` |
| 2. Verifica della mesh | `python3.12 -u tools/mesh_convergence.py` | decine di minuti | tabella ΔY SOFA / teoria della trave |
| 3. Scena interattiva | `./run_gui.sh` | – | GUI, grafici live, `results/gui_<data>.csv` |
| 4. Prove nominali | `python3.12 -u tools/run_cycle.py --free --contact` | circa 1 ora | H nominali, prove di contatto, `results/*.csv` |
| 5. Monte Carlo | `python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_finale` | diverse ore | `results/mc_finale/` |
| 6. Figure | `env -u PYTHONPATH python3 tools/make_figures.py` | secondi | `results/figures/*.png` |

I passi 1, 2 e 4 sono indipendenti dal Monte Carlo: la scena rigenera le sue mesh a ogni
avvio e il Monte Carlo contiene già i casi nominali. Il minimo per rifare i risultati è
0 → 5 → 6.

### 0. Installazione e ambiente

Prerequisiti:
- SOFA v26.06 (binario per Python 3.12) in `~/SOFA_v26.06.00_Linux_Python3.12/`;
- questo progetto in `~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project`;
- il Python di sistema con matplotlib, usato per i grafici e le figure.

```bash
sudo apt install python3-matplotlib python3-tk
cd ~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project
chmod +x run_gui.sh
source env.sh                                    # va ripetuto in ogni nuovo terminale
python3.12 -c "import Sofa.Core, Sofa.Simulation; print('SOFA OK')"     # atteso: SOFA OK
```

`env.sh` e `run_gui.sh` danno per scontati i percorsi indicati sopra. Se SOFA si trova altrove,
basta correggere `SOFA_ROOT` in `env.sh`, oppure lanciare `SOFA_ROOT=/percorso ./run_gui.sh`.

### 1. Generazione delle mesh

```bash
python3.12 tools/generate_meshes.py
```

| File | Contenuto |
|---|---|
| `meshes/tube_surface.obj` | superficie del tubo, modello visivo |
| `meshes/tube_collision.obj` | superficie del tubo, modello di collisione (solo il tratto vicino alla punta) |
| `meshes/tube_fem_surface.obj` | superficie esterna della mesh FEM (controllo) |
| `meshes/tube_volume.vtk` | mesh volumetrica del tubo: 3720 nodi, 11088 tetraedri, più fine (0.4 mm) nella parte flessibile |
| `meshes/soft_volume.vtk` | mesh volumetrica dell'oggetto soffice: 432 esaedri, più fini nella zona di contatto |
| `meshes/object_surface.obj` | superficie dell'oggetto di contatto |

- I file `.vtk` si aprono con ParaView.
- Lo script usa le stesse funzioni della scena (`modules/geometry_gen.py`).
- La scena rigenera le sue mesh a ogni avvio in `meshes/scena/`, una cartella ignorata da git, con i parametri del momento, compreso il profilo GUI più leggero.

### 2. Verifica della mesh (facoltativa)

```bash
python3.12 -u tools/mesh_convergence.py
```

Lo script prova `mesh_h_flex` = 0.8, 0.4, 0.3 e 0.2 mm. Per ciascun valore confronta la
deflessione della punta con la teoria della trave: momento T·yt dal cavo e forza laterale in
punta. Con 0.4 mm, il valore scelto, ΔY a 10 N vale 0.304 mm contro 0.320 mm teorici (−4.8 %).

### 3. Scena interattiva (GUI)

```bash
./run_gui.sh                                          # soffice, molla fissa
./run_gui.sh --object rigid --spring switched         # rigido, molla a due livelli
./run_gui.sh --object none                            # moto libero
./run_gui.sh --spring direct                          # attuazione diretta k = 50 N/mm
```

Il launcher si può lanciare da qualunque cartella e imposta l'ambiente da solo. Opzioni:
- `--object soft|rigid|none`
- `--spring fixed|switched|direct`
- `--kfree low|high` (ordine dei livelli della molla b)

Si aprono la finestra di SOFA e quella dei grafici. Premere **Play**.

| Tasto | Azione |
|---|---|
| Freccia su / `U` / `+` (finestra SOFA) | carico: u + 0.35 mm |
| Freccia giù / `D` / `-` (finestra SOFA) | scarico: u − 0.35 mm |
| `R` (finestra SOFA) | reset di comando e cavo |
| `M` / `A` (finestra SOFA) | modo manuale ↔ ciclo autonomo |
| `R` (finestra dei grafici) | azzera i loop e l'energia dissipata |

- A ogni comando il terminale stampa u, T, ΔLs, X, Y e ΔY della punta e la forza di contatto.
- La finestra dei grafici segue la struttura del Modulo 06 del corso:
  - a sinistra, il loop T–ΔLs della molla (legge di Hooke, banda ±0.85 N, energia dissipata) e il loop ΔY–u della punta con H calcolato in tempo reale;
  - a destra, nel tempo: posizione della punta, u e ΔLs, T, forza di contatto (soglia 5 mN, range 10–50 mN) e k.
- Ogni passo viene salvato in `results/gui_<data>.csv`.

Da dichiarare nel video:
- **Profilo GUI** (`gui_fast`, `gui_overrides`): mesh più leggere (`mesh_h_flex` 0.8 mm, circa l'85 % della rigidezza convergente; soffice con griglia uniforme 5) e componenti `Parallel…`. `dt` = 0.01 s e il solver dei vincoli (500 iterazioni, tolleranza 1e-7) restano quelli di consegna. Con `"gui_fast": False` la GUI usa il modello completo.
- **Deformazione mostrata ×50** (`visual_amplification`): la forma reale resta visibile in trasparenza.
- La simulazione va più lenta del tempo reale, quindi va dichiarato l'eventuale fattore di accelerazione del video.

### 4. Prove nominali (facoltative)

```bash
python3.12 -u tools/run_cycle.py --free                       # H di fixed, switched e direct
python3.12 -u tools/run_cycle.py --contact                    # rigido e soffice, entrambe le molle
python3.12 -u tools/run_cycle.py --free --contact --configs switched   # solo la molla (b)
```

- **Ciclo libero:** carico fino a u_max = 15 mm a 1 mm/s (oppure fino a T = 10 N), poi scarico.
- **Prova di contatto:** carico fino a 30 mN o a fine corsa, poi u fermo; 2.5 s di assestamento e 5 s di valutazione della stabilità.
- Le serie vanno in `results/free_<config>.csv` e `results/contact_<oggetto>_<config>.csv`.

### 5. Campagna Monte Carlo (127 casi)

I 127 casi sono (1 nominale + 20 realizzazioni) × 2 molle × 3 prove (libera, rigido, soffice),
più 1 caso di attuazione diretta. Seed 2026.

```bash
python3.12 -u tools/run_monte_carlo.py --dry-run --outdir /tmp/prova_mc     # prova senza SOFA (secondi)
python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_finale 2>&1 | tee results/mc_finale.log
python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_finale --resume   # dopo un'interruzione
```

- **Perturbazioni:** sei gruppi di parametri, ciascuno con un solo fattore: E (E_load ed E_unload), Fh_mat, Fh_spring, c, k (k_spring, k_low, k_high) e yt. La soglia di 5 mN non viene ritarata.
- **Ripresa:** ogni caso viene aggiunto a `cases.jsonl` appena finisce, quindi `--resume` riparte dal primo caso mancante.
- **Worker:** in una VM conviene usarne 2–3; con più worker si sono visti blocchi della VM.
- **Output:**
  - `summary.json`: statistiche, tornado e target;
  - `cases.jsonl`: metriche di ogni caso;
  - `realizations.json`: parametri estratti;
  - `series/`: serie temporali. In git restano solo `series/nominal_*.csv`, che servono alle figure.

Riepilogo rapido dal terminale:

```bash
python3 -c "import json; s=json.load(open('results/mc_finale/summary.json'))
for c,v in s['configs'].items(): print(c, {k: v[k] for k in ('H_nominal_pct','H_mean_pct','H_std_pct','H_worst_pct','robustness_index','unstable_total','dominant_parameters','meets_target')})
print('H diretta:', s['H_direct_nominal_pct'], '| piu robusta:', s['most_robust_config'])"
```

### 6. Figure

```bash
env -u PYTHONPATH python3 tools/make_figures.py --campaign results/mc_finale --out results/figures
```

Figure prodotte: `fig_tornado`, `fig_H_realizzazioni`, `fig_forze_contatto`,
`fig_isteresi_nominale` e `fig_contatto_nominale`. Lo script usa il Python di sistema
(matplotlib), non serve SOFA.

## Risultati (campagna `results/mc_finale`, seed 2026)

| | (a) fissa 1.2 N/mm | (b) due livelli 0.6/4.8 N/mm |
|---|---:|---:|
| H nominale [% di ΔYmax] | 14.38 | 16.00 |
| H medio ± dev. std | 15.19 ± 1.59 | 18.01 ± 2.24 |
| H caso peggiore | 18.09 (r19) | 21.46 (r19) |
| Indice di robustezza (peggiore / nominale) | **1.258** | 1.341 |
| H nominale / H attuazione diretta (26.75 %) | 0.54 | 0.60 |
| Contatti instabili (rigido + soffice, su 40) | **0** | 0 |
| Due parametri dominanti (tornado) | Fh_mat, Fh_spring | k, Fh_spring |
| Target (RI ≤ 1.3 e 20/20 stabili) | **raggiunto** | non raggiunto |

La configurazione più robusta è la **(a) fissa**.

**Ordine dei livelli della molla (b).** È stato provato anche l'ordine inverso (4.8 N/mm in moto
libero, 0.6 N/mm in contatto) sui casi nominali. Il contatto sul rigido migliora, con 30.2 mN
invece di 53.5 mN, ma H sale a 27.3 %, cioè 1.02 volte l'attuazione diretta: in moto libero la
molla non porta più alcun vantaggio. Poiché H è la metrica del progetto, si è scelto 0.6 → 4.8.

## Requisiti della traccia

| Requisito | Implementazione |
|---|---|
| Superficie OBJ per i modelli visivo e di collisione, mesh volumetrica per la meccanica | `modules/geometry_gen.py`, `tools/generate_meshes.py` (passo 1) |
| FEM con E diverso in carico e scarico | `TetrahedronFEMForceField`; il controller passa da E_load a E_unload secondo il verso del moto (`reinit()`) |
| `BarycentricMapping`, `CableConstraint` | `main_scene.py`: tendine con yt = 0.6 mm, ancorato alla punta |
| Molla in serie T = k·ΔLs ± isteresi | `modules/spring_model.py`: isteresi di Dahl (σ0 = 10/mm, ±(Fh_mat + Fh_spring)), smorzamento c, \|dT/dt\| ≤ 5 N/s |
| Controller da tastiera (su, giù, reset, manuale/autonomo, stampa u, T, ΔLs, punta) | `modules/controller.py` (`Sofa.Core.Controller`) |
| Grafici in tempo reale + CSV | `tools/live_plot_viewer.py`, `modules/live_plot.py` (passo 3) |
| Oggetti rigido e soffice (5 kPa, ν 0.45), `FrictionContactConstraint` μ = 0.2 | `main_scene.py`; stessa distanza iniziale (0.03 mm) per i due oggetti |
| H = max ΔY carico−scarico a pari u, in % di ΔYmax, rispetto a k = 50 N/mm | `modules/metrics.py: hysteresis_H`; caso `direct` |
| Contatto stabile: picco-picco < 30 % della media per 5 s, senza perdita | `modules/metrics.py: contact_stability` |
| 2 configurazioni × 20 realizzazioni ±(10–20) %, commutazione non ritarata | `modules/monte_carlo.py`, `tools/run_monte_carlo.py` (passo 5) |
| Media, dev. std, caso peggiore, indice di robustezza, tornado, contatti instabili | `results/mc_finale/summary.json`, `tools/make_figures.py` (passo 6) |

## Struttura

```text
main_scene.py          scena SOFA (createScene); profilo GUI se lanciata da runSofa
run_gui.sh, env.sh     launcher della GUI / ambiente per gli script python3.12
meshes/                mesh nominali (passo 1); meshes/scena/ = mesh di lavoro della scena
modules/
  config.py            parametri nominali, gruppi incerti, profilo GUI
  geometry_gen.py      mesh del tubo e dell'oggetto, export OBJ e VTK
  spring_model.py      molla in serie: Dahl, smorzamento, commutazione di k
  controller.py        controller da tastiera e ciclo autonomo, registrazione dei dati
  contact_force.py     forza di contatto da constraintForces del solver
  metrics.py           H e stabilità del contatto
  protocols.py         ciclo libero e prova di contatto (usati da run_cycle e Monte Carlo)
  monte_carlo.py       campionamento, regressione e tornado, riepilogo
  mc_worker.py         esecuzione di un caso del Monte Carlo
  live_plot.py         invio dei dati ai grafici, registrazione CSV (solo GUI)
  amplified_view.py    vista della deformazione amplificata (solo GUI)
tools/                 generate_meshes, mesh_convergence, run_cycle, run_monte_carlo,
                       make_figures, live_plot_viewer
results/               mc_finale/, figures/ (create dai passi 5 e 6)
```

## Ipotesi e limiti

- **Unità:** mm, N, s, MPa, massa in tonnellate.
- **Convenzione SOFA:** con `FreeMotionAnimationLoop`, valore del cavo e forze di vincolo sono impulsi, quindi si impone T·dt e si legge λ/dt.
- **Forza di contatto:** si legge dalle `constraintForces` del `BlockGaussSeidelConstraintSolver`, perché il `lambda` dei modelli di collisione mappati conta la forza due volte.
- **Glossario contro traccia:** dove i valori differiscono (E_unload 66 contro 65 GPa, ν 0.30 contro 0.33, k 1.13 contro 1.2 N/mm, parete 0.215 contro 0.22 mm) si sono usati i valori della traccia.
- **Mesh del soffice non del tutto convergente:** dimezzando il passo la forza cambia di circa −15 %. I confronti tra configurazioni restano validi; le forze assolute sul soffice sono approssimate.
- **Passaggio E_load → E_unload istantaneo:** a inizio scarico la punta sale di colpo (circa +18 % di cedevolezza).
- **H alto con molle rigide:** con k ≥ 4.8 N/mm a governare la tensione è il limite \|dT/dt\| ≤ 5 N/s. Probabilmente contribuisce anche il ritardo dovuto allo smorzamento numerico (`rayleighStiffness` = 0.1, circa 0.1 s), che fa crescere H sia per k = 4.8 sia per l'attuazione diretta.
- **Smorzamento c quasi ininfluente su H:** il ciclo è lento, 1 mm/s.
- **Stabilità facile da soddisfare:** il modello è quasi statico, quindi la stabilità non discrimina molto tra le configurazioni.
- **Forza di contatto fuori dal range 10–50 mN:**
  - molla (b) sul rigido in 16 casi su 20 (nominale 53.5 mN), per il salto istantaneo a 4.8 N/mm;
  - soffice in 7 casi su 20 per configurazione (circa 10 mN).

  La traccia indica questo range come range di interesse: non è un target del progetto 5.

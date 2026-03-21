# Bohr-s-atomic-model
# ⚛️ Bohr's Atomic Model — Hydrogen Atom Simulation

A physics simulation of the hydrogen atom based on Bohr's model, built in C with a Python visualizer. The simulation numerically integrates the two-body Coulomb interaction between an electron and a proton using a Velocity Verlet integrator, and models probabilistic quantum state transitions.

---

## 📁 Project Structure

```
Bohr-s-atomic-model/
├── Simulation.c          # Core physics engine (C)
├── Animation.py          # Orbit & energy visualizer (Python/matplotlib)
├── run.sh                # One-shot build + run script
├── setuparch.sh          # Arch Linux dependency setup
└── sim_essen/
    ├── particle/
    │   └── particle.h    # Particle struct definition
    ├── acc_cal/
    │   └── acc.c / acc.h # Coulomb force & acceleration calculator
    └── verlet/
        └── integrator.c / integrator.h  # Velocity Verlet integrator
```

---

## ⚙️ How It Works

### Physics Engine (`Simulation.c`)
- Takes the **electron shell state n** as input (n = 1, 2, 3, ...)
- Places the electron at the correct Bohr radius: `r = 0.53 × n² Å`
- Sets the initial orbital velocity: `v = 2.18×10⁶ / n m/s`
- Runs a **two-body simulation** — the proton is not fixed; it receives a recoil velocity consistent with conservation of momentum
- Integrates motion using a **split Velocity Verlet scheme** for excellent energy conservation
- At each step, computes kinetic energy, potential energy, and total energy
- Logs everything to `orbit.csv`
- Includes a **probabilistic quantum jump mechanic**: the electron can spontaneously decay from state `n` to `n-1`, printing the emitted photon energy in eV

### Visualizer (`Animation.py`)
- Reads `orbit.csv` output from the C simulation
- Renders a **side-by-side matplotlib animation**:
  - **Left panel**: Electron and proton orbits with an electron trail
  - **Right panel**: Real-time total energy plot (in eV)

---

## 🚀 Getting Started

### Prerequisites

**C compiler and Python packages:**
```bash
# Debian/Ubuntu
sudo apt install gcc python3 python3-pip
pip install numpy pandas matplotlib

# Arch Linux (use the provided script)
bash setuparch.sh
```

### Build & Run

**Option 1 — Use the shell script (recommended):**
```bash
chmod +x run.sh
./run.sh
```
The script will ask for the electron's orbit (shell number), compile the C code, run the simulation, and launch the animation automatically.

**Option 2 — Manual:**
```bash
# Compile
gcc -I. sim_essen/acc_cal/acc.c sim_essen/verlet/integrator.c Simulation.c -o simulation -lm

# Run simulation
./simulation
# Enter state of electron: 3

# Visualize
python Animation.py
```

---

## 📊 Output

The simulation writes `orbit.csv` with columns:

| Column | Description |
|--------|-------------|
| `step` | Simulation step index |
| `ex`, `ey` | Electron x, y position (m) |
| `px`, `py` | Proton x, y position (m) |
| `ke` | Kinetic energy (J) |
| `pe` | Potential energy (J) |
| `te` | Total energy (J) |

Quantum transitions are printed to stdout:
```
Transition 3 -> 2 | Photon = 1.89 eV
```

---

## 🔬 Physical Constants Used

| Constant | Value |
|----------|-------|
| Proton mass | 1.67262193 × 10⁻²⁷ kg |
| Electron mass | 9.1093837 × 10⁻³¹ kg |
| Elementary charge | 1.60217663 × 10⁻¹⁹ C |
| Coulomb's constant k | 8.98755 × 10⁹ N·m²/C² |
| Bohr radius (n=1) | 0.53 × 10⁻¹⁰ m |

---

## 🧠 Design Decisions

- **Velocity Verlet integrator**: Chosen for its symplectic nature — it conserves energy over long integration times far better than simple Euler methods.
- **Two-body treatment**: The proton is given a recoil velocity rather than being fixed, making the simulation physically correct for a free hydrogen atom.
- **Modular `sim_essen` library**: Physics components (force calculation, integrator, particle struct) are separated from the main simulation logic for clarity and reusability.
- **Probabilistic quantum jumps**: Modeled using an exponential decay probability per timestep, loosely analogous to spontaneous emission.

---

## 📌 Known Limitations

- Bohr's model is a semiclassical approximation — it does not capture quantum mechanical effects like orbital shapes or spin.
- The quantum jump resets the electron position instantaneously (as the model requires), rather than as a continuous process.
- Only downward transitions (n → n−1) are currently modeled.

---

## 📄 License

This project is open source. Feel free to use, modify, and build upon it.

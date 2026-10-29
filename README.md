# Vapor-Injection-Heat-Pump-Optimization-for-Electric-Vehicle-Cold-Weather-Performance
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![CoolProp](https://img.shields.io/badge/CoolProp-6.4.3-green.svg)](http://www.coolprop.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

##  Project Overview
This repository contains a steady-state thermodynamic model evaluating a transcritical CO₂ (R744) heat pump designed for Electric Vehicle (EV) cabin heating in cold weather (down to -15 °C). 

The project evaluates and compares two system architectures:
1. **Case A (Baseline):** A conventional vapor-compression heat pump cycle.
2. **Case B (FT-VI):** A Flash-Tank Vapor Injection system that injects intermediate-pressure vapor into a two-stage compressor to improve heating capacity and COP while cooling the compressor discharge.

### Objectives & Optimization Framework
* **Decision Variables:** Intermediate flash-tank pressure ($P_{inj}$), high-side discharge pressure ($P_{dis}$), and compressor speed.
* **Objectives:** Maximize Heating COP, maximize heating capacity, and minimize compressor discharge temperature.
* **Constraints:** Maintain compressor discharge temperature $\le$ 130 °C, discharge pressure $\le$ 12 MPa, and ensure cabin heating capacity meets the 5.0 kW target.

---

##  Thermodynamic Methodology

### Phase 1: Baseline Model (Case A)
The conventional cycle models four main state points:
* **1 to 2:** Isentropic compression (adjusted by $\eta_{is} = 0.70$).
* **2 to 3:** Isobaric heat rejection in the gas cooler.
* **3 to 4:** Isenthalpic expansion via an Electronic Expansion Valve (EEV).
* **4 to 1:** Isobaric heat absorption in the evaporator.

*Heating Capacity:* $Q_h = \dot{m} (h_2 - h_3)$  
*Compressor Work:* $W = \dot{m} (h_2 - h_1)$

### Phase 2: FT-VI Model (Case B)
The FT-VI system utilizes two expansion stages and a flash tank to separate saturated liquid and vapor at an intermediate pressure ($P_{inj}$):
* **Mass Balance:** The vapor quality at the flash tank inlet dictates the injection ratio. Vapor fraction $x = (h_6 - h_f) / (h_g - h_f)$.
* **Injection Ratio:** $i = (Injected Vapor Mass Flow) / (Suction Mass Flow)$
* **Two-Stage Compression:** Suction gas is compressed to $P_{inj}$, mixed adiabatically with the injected saturated vapor, and then compressed to $P_{dis}$. This lowers the final discharge temperature and boosts total mass flow through the gas cooler.

### Phase 3:





###  Phase 4:




### Final Phase :







### To be Continued
---
## Project Structure

The repository is modularly designed with the following file structure:

```text
Vapor-Injection-Heat-Pump-Optimization-for-Electric-Vehicle-Cold-Weather-Performance/
│
├── run_all.py            # Main orchestration script: runs all sweeps, saves data, plots figures
├── baseline_cycle.py     # Conventional cycle & compressor displacement calibration
├── flash_tank_vi.py      # FT-VI cycle, mass/energy balances, P_inj sweeps
├── co2_props.py          # Fluid property backend: Wraps CoolProp & built-in Peng-Robinson fallback
├── config.py             # Global configuration: Project brief data, parameters, and assumptions
├── README.md             # Project documentation
│
├──RESULTS ├── /results/             # Output CSV files and JSON handoff data
    └── /figures/             #  P-h diagrams, T-s diagrams, maps, and schematics

```

##  Installation & Setup

To run this simulation locally, you will need **Python 3.8+** and a few scientific libraries.

**1. Clone the repository**
```bash
git clone https://github.com/durgeshnimse/Vapor-Injection-Heat-Pump-Optimization-for-Electric-Vehicle-Cold-Weather-Performance.git
cd Vapor-Injection-Heat-Pump-Optimization-for-Electric-Vehicle-Cold-Weather-Performance
```

**2. Create and activate a virtual environment (Recommended)**
```bash
* **Windows:**
```bash
python -m venv venv
venv\Scripts\activate

```


* **macOS/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate

```



**3. Install the required dependencies**
The primary dependencies are `numpy`, `scipy`, `matplotlib`, and `CoolProp`. Install them via pip:

```bash
pip install numpy scipy matplotlib CoolProp

```

*(Note: CoolProp is highly recommended for exact transcritical R744 property evaluations. If it fails to install, the code will seamlessly fall back to a built-in Peng-Robinson Equation of State model).*

---

##  Running the Code

The entire simulation suite is automated. To execute the models, generate the CSV data tables, and create the optimization plots, simply run:

```bash
python run_all.py

```

### Forcing a Specific Thermodynamic Backend

If you want to explicitly test the `CoolProp` (Span-Wagner EOS) backend vs. the built-in `Peng-Robinson` (PR) backend, use the `CO2_BACKEND` environment variable:


**Force CoolProp:**
* Windows (CMD): `set CO2_BACKEND=COOLPROP && python run_all.py`
* macOS/Linux: `CO2_BACKEND=COOLPROP python run_all.py`


**Force Built-in PR Fallback:**
* Windows (CMD): `set CO2_BACKEND=PR && python run_all.py`
* macOS/Linux: `CO2_BACKEND=PR python run_all.py`



---

##  Expected Outputs

Upon running `run_all.py`, the script prints a console summary of the optimal operating points and generates two directories:

### 1. `/results`

Contains CSV tables documenting every parametric sweep (e.g., intermediate pressure sweeps at -5 °C, -10 °C, and -15 °C) and a `handoff_students_3_4.json` file containing the calculated state points required for later project phases.

### 2. `/figures`

* `fig_ph.png` & `fig_ts.png`: Thermodynamic state diagrams mapping the Baseline and FT-VI cycles.
* `fig_schematic.png`: An automated system architecture diagram with state-point labels.
* `fig_B_map.png`: Heatmaps plotting FT-VI COP and Capacity against $P_{inj}$ and $P_{dis}$.
* `fig_A_Pdis.png` & `fig_B_Pinj.png`: Performance curves showing the effects of pressure variations.

---

## Configuration

To test different design limits or ambient temperatures, modify the dictionaries directly in `config.py`:

* `BRIEF`: Contains hard constraints (e.g., $T_{amb} = -10$, $Q_{load} = 5.0$ kW, $T_{dis\_max} = 130$ °C).
* `ASSUME`: Contains thermodynamic assumptions (e.g., superheat $= 5.0$ K, isentropic efficiency, volumetric efficiency factors).

---

##  References

1. Wang, D., et al. (2018). *Heating performance evaluation of a CO₂ heat pump system for an electrical vehicle at cold ambient temperatures.* Applied Thermal Engineering.
2. Mohaideen Abdul Khader, S. O., et al. (2025). *Comprehensive review of vapor injection heat pumps in electric vehicles.* Applied Thermal Engineering.
```

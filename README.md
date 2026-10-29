# Vapor-Injection-Heat-Pump-Optimization-for-Electric-Vehicle-Cold-Weather-Performance

# Transcritical R744 EV Heat Pump Optimization: Baseline vs. FT-VI

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains a steady-state thermodynamic model evaluating a transcritical CO2 (R744) Electric Vehicle heat pump. It compares a conventional cycle (Case A) against a Flash-Tank Vapour-Injection (FT-VI) system (Case B)[cite: 2, 5]. The toolset optimizes the intermediate injection pressure ($P_{inj}$) to maximize heating Coefficient of Performance (COP) and capacity while strictly maintaining the compressor discharge temperature below 130 °C.

![System Schematic](figures/fig_schematic.png)
*(Note: The system schematic above is generated automatically when the code is executed.)*

## Features
- **Parametric Sweeps:** Evaluates performance across compressor speeds (2500–8000 rpm), high-side pressures (8.0–12.0 MPa), and intermediate injection pressures (3.5–6.0 MPa)[cite: 1, 4].
- **Accurate Fluid Properties:** Defaults to `CoolProp` (Span-Wagner EOS) for precise transcritical CO2 evaluation, with a built-in Peng-Robinson cubic EOS fallback.
- **Automated Diagram Generation:** Automatically outputs P-h diagrams, T-s diagrams, system schematics, and Pareto contour maps.
- **Design Constraints:** Ensures system feasibility by checking bounds, such as $P_{evap} < P_{inj} < P_{dis}$ and $T_{dis} \le 130 ^\circ \text{C}$[cite: 4, 5].

---

## Repository Structure

* **`run_all.py`**: The orchestration script that executes all calculations for Students 1 and 2, writes CSV/JSON outputs, and generates all figures[cite: 1].
* **`baseline_cycle.py`**: Calculates the Student 1 baseline cycle and calibrates the compressor displacement for a target heating load[cite: 2].
* **`flash_tank_vi.py`**: Implements the Student 2 FT-VI model, including adiabatic flash tank mass/energy balances and optimal intermediate pressure sweeps[cite: 5].
* **`co2_props.py`**: The R744 property backend handling state-point queries. 
* **`config.py`**: Stores all project brief data (ambient temperatures, target loads) and modelling assumptions (isentropic efficiency, superheat, volumetric efficiency trends).

---

## Prerequisites

- **Python 3.8+**
- `numpy`[cite: 1, 3]
- `scipy`[cite: 1, 3]
- `matplotlib`[cite: 1]
- `CoolProp` (Highly recommended for the most accurate property calculations)[cite: 3]

## Installation

**1. Clone the repository:**
```bash
git clone [https://github.com/your-username/ev-heatpump-optimization.git](https://github.com/your-username/ev-heatpump-optimization.git)
cd ev-heatpump-optimization

python -m venv venv
venv\Scripts\activate

python3 -m venv venv
source venv/bin/activate

python3 -m venv venv
source venv/bin/activate

pip install numpy scipy matplotlib CoolProp
python run_all.py

# Real Gas Thermodynamic Property & Fugacity Engine (Cubic EOS)

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A lightweight computational thermodynamics solver in Python for predicting real-gas PVT behavior, energy departure functions ($H^R, S^R$), and pure-component fugacity ($\phi, f$) using classical cubic Equations of State (vdW, RK, SRK, PR).

Formulations and critical constants are implemented directly from **J. M. Smith, H. C. Van Ness, M. M. Abbott, and M. T. Swihart**, *Introduction to Chemical Engineering Thermodynamics* (8th Ed., McGraw-Hill, 2018).

---

## Features

- **4 Cubic EOS Models:** van der Waals (vdW, 1873), Redlich–Kwong (RK, 1949), Soave–Redlich–Kwong (SRK, 1972), and Peng–Robinson (PR, 1976).
- **10 Built-In Industrial Gases:** Methane, Ethane, Propane, *n*-Butane, $\text{CO}_2$, $\text{N}_2$, $\text{O}_2$, Ethylene, Ammonia, and Water (Appendix B, Table B.1).
- **Full Property Spectrum:** Gas compressibility factor $Z$, molar volume $V$, residual enthalpy $H^R$, residual entropy $S^R$, fugacity coefficient $\phi$, and fugacity $f$.
- **Automated Visualization:** Generates publication-ready 4-quadrant plots ($Z$, $f$, $H^R$, $S^R$ vs. $P$ from 1 to 100 bar).

---

## Governing Thermodynamic Relations

| Property | Formula | Reference |
| :--- | :--- | :--- |
| **Generic Cubic EOS** | $P = \frac{RT}{V - b} - \frac{a(T)}{(V + \epsilon b)(V + \sigma b)}$ | SVNA Eq. (3.41) |
| **Cubic in $Z$** | $Z^3 + [(\epsilon + \sigma - 1)\beta - 1]Z^2 + [\epsilon\sigma\beta^2 - (\epsilon + \sigma)\beta(1+\beta) + q\beta]Z - [\epsilon\sigma\beta^2(1+\beta) + q\beta^2] = 0$ | SVNA Eq. (3.48) |
| **Departure Integral ($I$)** | $I = \frac{1}{\sigma - \epsilon}\ln\left[\frac{Z + \sigma\beta}{Z + \epsilon\beta}\right] \quad (\sigma \neq \epsilon); \qquad I = \frac{\beta}{Z} \quad (\text{vdW})$ | SVNA Eq. (13.72) |
| **Residual Enthalpy** | $\frac{H^R}{RT} = Z - 1 + \left[\frac{d\ln\alpha}{d\ln T_r} - 1\right] qI$ | SVNA Eq. (13.75) |
| **Residual Entropy** | $\frac{S^R}{R} = \ln(Z - \beta) + \left[\frac{d\ln\alpha}{d\ln T_r}\right] qI$ | SVNA Eq. (13.76) |
| **Fugacity Coefficient** | $\ln\phi = \frac{G^R}{RT} = Z - 1 - \ln(Z - \beta) - qI, \qquad f = \phi P$ | SVNA Eqs. (10.33, 13.85) |

*Dimensionless parameters:* $\beta = \Omega \frac{P_r}{T_r}$ and $q = \frac{\Psi\alpha(T_r;\omega)}{\Omega T_r}$.

---

## Quickstart

### Installation
```bash
git clone https://github.com/your-username/cubic-eos-thermodynamics.git
cd cubic-eos-thermodynamics
pip install numpy matplotlib
```

### Run
```bash
python thermo_project.py
```

Follow the interactive prompts: choose a gas (e.g., `methane`), enter $T$ (K) and $P$ (bar), select an EOS model (`PR` or `ALL`), and optionally export the 4-panel diagnostic plot.

---

## Verification & Benchmark

Validated against published examples in SVNA (8th Ed.):

1. **Saturated *n*-Butane Vapor at $350\text{ K}$, $9.4573\text{ bar}$ (Ex. 3.9, p. 101):**
   - Experimental $V = 2482\text{ cm}^3\cdot\text{mol}^{-1}$
   - Calculated: $\text{vdW} = 2667$, $\text{RK} = 2555$, $\text{SRK} = 2520$, $\text{PR} = 2486\text{ cm}^3\cdot\text{mol}^{-1}$ (exact match).
2. **Supercritical *n*-Butane at $500\text{ K}$, $50\text{ bar}$ (Ex. 13.5 / Table 13.7, p. 490):**
   - Exact numerical reproduction of $Z$, $H^R$, and $S^R$ across all models (e.g., $\text{PR: } H^R = -4988\text{ J}\cdot\text{mol}^{-1}$, $S^R = -7.426\text{ J}\cdot\text{mol}^{-1}\cdot\text{K}^{-1}$).

---

## Authors & Course Context

- **Institution:** Indian Institute of Technology Delhi (IIT Delhi)
- **Course:** Chemical Engineering Thermodynamics $\cdot$ Term Computer Project
- **Instructor:** Prof. S. K. Gupta
- **Contributors:**
  - **Bansiwal Ankit** (`2025CH71777`) — `ch7251777@chemical.iitd.ac.in`
  - **Gupta Rohit Anil** (`2025CH71467`) — `ch7251467@chemical.iitd.ac.in`

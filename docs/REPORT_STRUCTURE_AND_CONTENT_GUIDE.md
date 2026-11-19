# CPR2 — Report Structure and Content Guide

**COMET: A Dual Swashplate Autonomous Coaxial Bi-Copter**

This companion guide maps the chapters of the report to their contents to facilitate academic review. Chapters 1–9 assess the published paper; Chapter 10 documents the independent MATLAB numerical reproduction. Paper measurements must be distinguished from simulated results.

## Part I — Critical Review of the Research Paper

| Chapter | Content |
| --- | --- |
| **1. Main Claims** | Research question; COMET architecture; maneuverability, propulsion efficiency, endurance, tracking and navigation claims. |
| **2. Significance** | Engineering relevance of dual-swashplate actuation and increased control authority. |
| **3. Evidence Provided** | Propulsion, torque, structural, endurance, tracking, saturation, power and navigation experiments. |
| **4. Assessment of Evidence** | Methodological validity, experimental comparison, uncertainty, repeatability and generalizability. |
| **5. Drawbacks & Improvements** | Limitations and potential improvements in aerodynamics, actuators, allocation, planning and robustness. |
| **6. Key Methods & Tools** | Swashplate mechanics, Newton–Euler dynamics, rotor wrench, differential-flatness-based planning, cascaded control, mixer, sensors and test tools. |
| **7. Critical Assessment** | Scientific contribution, strengths, limitations and supported conclusions. |
| **8. Relation to Course** | Underactuation, differential flatness, controllability, control allocation and saturation. |
| **9. Related Research** | FLOAT Drone, MOMAV and related control-allocation approaches. |

## Part II — Independent MATLAB Reproduction

| Section | Content |
| --- | --- |
| **10.1 Objectives and Methodology** | Scope and controlled Single/Dual evaluation. |
| **10.2 MATLAB Architecture** | Responsibilities of `comet_parameters.m`, `simulate_COMET.m` and `run_COMET_comparison.m`. |
| **10.3 Nonlinear Dynamic Model** | 19 numerical state components, Newton–Euler equations, quaternion attitude, rotor wrench, actuator dynamics, drag and power. |
| **10.4 Control Architecture** | Cascaded loops, feedforward, locally updated cyclic allocation, saturation and anti-windup. |
| **10.5 Numerical Implementation** | Trajectories, RK4, sampled feedback, delays and performance metrics. |
| **10.6 Calibration Campaign** | Parameter sweeps and final settings. |
| **10.7 Final Simulation Results** | Five trajectories in Single and Dual modes; tracking, actuation, saturation and power. |
| **10.8 Experimental Comparison** | Comparison with the paper's measurements and limitations. |
| **10.8A Supplementary Assessment** | Small baseline Single/Dual tracking gap and exploratory hybrid differential flatness. |
| **10.9 Discussion and Conclusions** | Capabilities, limitations and future work. |

## Supporting Material

The private repository includes executable MATLAB files under `MATLAB_MAIN_FILES/`, figures, a technical report and a MATLAB usage guide. Optional calibration and flatness studies are under `research/`.

**Interpretation:** The numerical baseline reproduces differing cyclic actuation demand more clearly than the experimental tracking gap. The paper's data are comparison targets, not simulated measurements.

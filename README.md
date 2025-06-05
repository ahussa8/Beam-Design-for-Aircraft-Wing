# 🛩️ ME-212 Solid Mechanics – Aircraft Wing Spar Design

This project focuses on designing a prismatic aircraft wing spar using analytical, numerical, and simulation methods. We analyzed distributed aerodynamic loads, verified compliance with FAA safety regulations, and optimized the spar geometry for strength and weight.

---

## 📌 Overview

Designed for a new aircraft with a 4m cantilevered wing, this project aimed to:
- Ensure structural compliance with FAA Part 25.303
- Maintain a safety factor of 1.5 against yielding
- Minimize spar weight using efficient cross-sectional design

---

## 🔧 Methods & Analysis

### Part 1: Prismatic Spar Design
- Modeled distributed lift using a parabolic load distribution equation
- Determined internal shear and moment using integration
- Designed I-beam cross-section with specified flange and web limits
- Verified stress using von Mises theory and ensured safety factor compliance

### Part 2: Simulation & Validation
- Modeled the spar in Autodesk Inventor and applied equivalent point loads
- Ran FEA stress simulations and probed stress at critical regions
- Compared theoretical and simulated von Mises stress values

### Part 3: Safety Factor Analysis
- Used MATLAB to plot internal shear force, bending moment, and safety factor over the span
- Explored variable flange width to maintain a uniform safety factor
- Proposed weight-optimized custom spar design (extra credit)

---

## 📊 Key Outputs

- Internal Shear and Moment Diagrams
- Stress Distribution in Autodesk Inventor
- Weight Calculation of Spar in Kilograms
- Safety Factor Profile Over Wing Span
- Optimization of flange width vs. location

---

## 🧰 Tools & Technologies

- MATLAB (analytical and numerical integration)
- Autodesk Inventor (FEA & CAD modeling)
- Excel (material database and force calculations)
- Wolfram Alpha & Symbolab (integral verification)

---

## 📁 Project Files

- `beam_design.m`: MATLAB code for shear/moment & safety factor plots
- `spar_model.ipt`: Autodesk Inventor CAD file
- `spar_simulation.idw`: FEA simulation results
- `Project_Report.pdf`: Full documentation and analysis

---

## 🛠️ Engineering Principles

- Solid Mechanics (Stress & Bending Theory)
- Structural Optimization
- Finite Element Analysis (FEA)
- Safety & Regulatory Compliance
- Design for Manufacturing

---

## 🚀 Extra Credit Features (Optional)

- Custom tapered spar design
- Manufacturing-aware redesign (sheet metal vs. machined aluminum)
- Multi-point loading simulations
- LaTeX report via Overleaf with figures and citations

---

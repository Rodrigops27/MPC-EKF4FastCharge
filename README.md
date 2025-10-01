# Lithium-Ion Cell-Level Optimal Fast Charge Algorithm using MPC and EKF with a Physics-Based Reduced Order Model
This repository is distributed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0).

This project is an implementation of a **fast-charging control algorithm** for lithium-ion cells using a **physics-based reduced-order model (ROM)**, an **Extended Kalman Filter (EKF)** for state estimation, and a **Model Predictive Controller (MPC)** with inequality constraints solved via **Hildreth’s algorithm**.

The code targets the NMC30 (3.7 V, 30 Ah) lithium-ion cell, using the physics-based reduced-order model developed by Prof. Gregory Plett and Prof. Scott Trimboli \[1], and demonstrates constrained MPC for safe and accelerated charging.

---

## Background

This repository builds on the computational framework described in \[2], extending it to a working MATLAB implementation for fast-charge protocols. The main goal is to:

* Exploit internal electrochemical states. Specifically the side-reaction overpotential (which signals lithium plating), rather than relying solely on voltage-based constraints..
* Enable **pseudo-minimum-time charging**, where the cost function is defined as the normed distance between the predicted output and the reference trajectory to approach the min-time solution. 
* Use an **Extended Kalman Filter (EKF)** to provide linearized state feedback of electrochemical variables, so that MPC can apply constraints on the unmeasurable electrochemical states.
* Provide a modular simulation structure combining model, estimator, and optimizer.

---

## Algorithm Description

The control loop consists of the following major steps each sampling instant:

1. **ROM Simulation (OB\_step)**

   * The reduced-order model (ROM) simulates cell voltage and internal states given the applied current.
   * Implemented in [`OB_step.m`](src/MPC-EKF4FastCharge/OB_step.m).

2. **State Estimation (iterEKF)**

   * The Extended Kalman Filter updates estimates of SOC, electrode stoichiometries, potentials, etc., based on ROM predictions and measured voltage.
   * Initialization via [`initKF.m`](src/UTILITY/initKF.m). Iterative update via [`iterEKF.m`](src/UTILITY/iterEKF.m).

3. **MPC Linearization (EKFmatsHandler)**

   * Converts EKF-estimated states into a linearized state-space form for MPC prediction.
   * Picks the **closest local ROM (no blending)** and builds the linearization around it.
   * Returns the state-space matrices required for MPC and constraint handling, **before Δu augmentation**.
   * Implemented in [`EKFmatsHandler.m`](src/MPC-EKF4FastCharge/EKFmatsHandler.m).

4. **Prediction Matrices (predMat)**

   * Build horizon-based prediction matrices Φ and G for SOC, voltage, and overpotential.
   * Implemented in [`predMat.m`](src/MPC-EKF4FastCharge/predMat.m).

5. **MPC Optimization (iterMPC)**

   * Solve the quadratic program to compute the optimal input sequence subject to constraints:

     * Current magnitude and slew limits
     * Voltage bounds
     * Side-reaction overpotential bound (plating constraint)
   * Implemented in [`iterMPC.m`](src/MPC-EKF4FastCharge/iterMPC.m), using [`constraintsMPC.m`](src/MPC-EKF4FastCharge/constraintsMPC.m) to assemble inequalities.

6. **QP Solver (Hildreth)**

   * If constraints are active, the quadratic program is solved using Hildreth’s dual coordinate ascent method.
   * Implemented in [`hildreth.m`](src/MPC-EKF4FastCharge/hildreth.m).

7. **Main Driver (runMPC)**

   * Executes the full loop (ROM → EKF → MPC → current command).
   * Implemented in [`runMPC.m`](src/MPC-EKF4FastCharge/runMPC.m).

---

## Repository Structure

* **Main simulation**

  * `runMPC.m` — entry point to run a full charging experiment
* **ROM simulation**

  * `outBlend.m` — full profile ROM simulation
  * `OB_step.m` — single-step ROM for MPC
* **Estimation**

  * `initKF.m` — initialize EKF
  * `iterEKF.m` — EKF iteration
* **MPC preparation**

  * `EKFmatsHandler.m` — builds augmented plant for MPC
  * `predMat.m` — constructs prediction matrices
* **MPC optimization**

  * `initMPC.m` — initialize MPC data structure
  * `iterMPC.m` — single MPC iteration
  * `constraintsMPC.m` — build inequality matrices
  * `hildreth.m` — dual QP solver

---

## Results

Simulation demonstrates:

* Charging from 10% to 95% SOC.
* Enforced safety constraints on **voltage, current, and overpotential**.
* Adaptive current profiles that balance charging speed with electrochemical limits.
* Maximum charge current constrained to **2C (≈ 59.72 A)**.
* Maximum cell voltage limited to **4.1 V**.
* Minimum side-reaction overpotential threshold set at **0.08 V**.
* MPC prediction horizon **Np = 5** and control horizon **Nc = 2**.

### Fast Charge Protocols Discussion
![Fast Charge Protocol with Side-Reaction Constraint](assets/MPCEKF1095Phise.png)
Fast Charge Protocol with applied current, voltage and side-reaction overpotential constraints.

![Fast Charge Protocol without Side-Reaction Constraint](assets/MPCEKF1095WOPhise.png)
Fast Charge Protocol with applied current and voltage but without side-reaction overpotential constraints. This plot is on development, where a PDE solver was used to complement the ROM solution, since the ROM becomes less accurate beyond 80% SOC.

When the **side-reaction overpotential constraint is enforced**, the MPC adaptively reduces the applied current as the electrode potentials approach the plating threshold, leading to a safer but slightly slower charge. At the beginning of the simulation, the charging profile resembles a **Constant-Current/Constant-Voltage (CC–CV)** trajectory. However, as the MPC foresees a potential violation of the side-reaction overpotential constraint, it proactively backs off the current to respect the limit — as can be seen in the figure. Subsequent current calculations then allow charging to continue as aggressively as possible while still honoring the imposed constraint.

By contrast, when **no overpotential constraint is imposed**, the controller follows a classic **constant-current/constant-voltage (CC–CV) charging profile**, which completes the 10–95% SOC charge in roughly **1,500 seconds**. While faster, this trajectory risks lithium plating.

In both cases, the **cell voltage limit is reached just prior to achieving the SOC target**, which is typical of lithium-ion fast-charge protocols. Importantly, these results illustrate that **operating within electrochemical limits extends lifetime but comes at the cost of additional charging time** compared to unconstrained CC–CV operation.

---

## Conclusions & Next Steps

* This framework enables **constraint-aware, physics-informed charging protocols**.
* Outlook: thermal coupling, experimental validation and more detailed degradation mechanisms.
* The potential viability of this method is inferred from both simulation and experimental results, which have also demonstrated four charging strategies similar in principle to the one presented here \[3].

---

## References

## References  

[1] G. L. Plett, & M. S. Trimboli, *Battery Management Systems, Volume III: Physics-Based Methods*. Artech House, 2024.  
[2] M. A. Xavier, A. K. de Souza, K. Karami, G. L. Plett, and M. S. Trimboli, *A Computational Framework for Lithium-Ion Cell-Level Model Predictive Control Using a Physics-Based Reduced-Order Model*, ACC 2021.  
[3] Y. Li, L. Aldrich, K. Stetzel, J. Lee, G. L. Plett, *Optimal fast charging of lithium-ion cells under electrochemical-thermal constraints*, Applied Energy, 269, 2020, 115127.  

---

## Acknowledgments & Copyright

The physics-based reduced-order model (ROM) for the NMC30 cell, as well as the Extended Kalman Filter (EKF) implementation (initKF.m and iterEKF.m), were originally developed by Prof. Gregory L. Plett and Prof. M. Scott Trimboli.

The ROM single-step simulator [OB_step.m] included here is derived from their [outBlend.m] framework.

This adaptation is provided with attribution for academic and educational purposes.
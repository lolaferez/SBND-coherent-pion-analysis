# Nuclear Effects in Coherent Pion Production at SBND

Bachelor's Thesis project developed as part of the Double Degree in Physics and Mathematics at the University of Granada.

This project studies the impact of **nuclear final-state interactions (FSI)** on the identification of coherently produced charged pions in neutrino–argon interactions at the **Short-Baseline Near Detector (SBND)** experiment at Fermilab.

## Project Overview

Coherent pion production is a neutrino–nucleus interaction in which the nucleus remains in its ground state. Experimentally identifying these events requires understanding how nuclear effects can modify the particles observed in the detector.

The project combines theoretical particle physics, mathematical methods and computational data analysis to investigate which observables are most sensitive to these nuclear effects.

The analysis compares two simulated datasets:

- **NO FSI:** 25,933 simulated events without final-state interactions
- **FSI:** 26,023 simulated events including nuclear final-state interactions

The simulations were generated with **GENIE**, using the **INTRANUKE hA2018** model for final-state interactions.

## Data Analysis

The computational analysis was performed using **C++ and CERN ROOT**.

The workflow included:

- Processing approximately **50,000 simulated neutrino events**
- Comparing samples with and without nuclear final-state interactions
- Exploring a large number of kinematic and particle-level variables
- Implementing event-by-event coordinate transformations relative to the incoming neutrino direction
- Analysing particle multiplicities and final-state composition
- Studying momentum and angular distributions
- Identifying observables sensitive to nuclear effects
- Producing statistical summaries and ROOT visualisations

The analysis focuses primarily on **charged-current coherent pion production**, where the expected final-state topology contains a charged lepton and a charged pion.

## Main Results

The analysis shows that final-state interactions primarily affect the **hadronic sector**, while the outgoing muon kinematics remain essentially unchanged.

Some of the main results are:

- **22.0%** of events in the FSI sample lose the charged pion through absorption inside the argon nucleus.
- After detector-motivated reconstruction requirements, **79.6%** of FSI events still retain a two-track topology, showing that track multiplicity alone is not sufficient to identify nuclear effects.
- Nuclear interactions broaden the **total transverse momentum** distribution.
- Using a threshold around `pT_tot = 0.2 GeV`, approximately **10.3%** of FSI events with a surviving charged pion migrate outside the region characteristic of the NO-FSI sample.
- FSI broaden the pion angular distribution relative to the incoming neutrino, particularly for pion momenta below approximately **1 GeV**.
- Using an angular threshold of `θ_νπ = 60°`, approximately **9.1%** of FSI events move outside the angular region characteristic of the NO-FSI sample.

Overall, the most sensitive observables were found to be **charged-pion survival, total transverse momentum and pion angular variables**.

## Physics Background

The thesis also develops the theoretical framework behind the analysis, including:

- Standard Model and neutrino physics
- Lie groups and Lie algebras
- SU(2) and SU(3) representations
- Quark model and meson multiplets
- Coherent pion production
- PCAC theorem
- Glauber model for pion–nucleus scattering
- Liquid Argon Time Projection Chambers (LArTPCs)
- The SBND experiment and Booster Neutrino Beam

## Technologies

- **C++**
- **CERN ROOT**
- **GENIE simulated data**
- Statistical data analysis
- Data visualisation
- Numerical and mathematical modelling

## Repository Structure

```text
.
├── Codes/        # C++ / ROOT analysis code
├── Data/            # Data information or sample files
├── plotis/         # Selected plots and analysis results
├── Ferez_Martinez_MDolores_TFG/          # Bachelor's Thesis PDF
├── prese.pdf/          # Bachelor's Thesis Presentation
└── README.md

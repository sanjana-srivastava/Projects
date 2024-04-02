# Simplified 6-DOF Trajectory Simulations for Stardust SRC and Space Shuttle
This project aims to develop a simplified six-degree-of-freedom (6-DOF) model for atmospheric reentry trajectory simulations of the Stardust Sample Return Capsule (SRC) and NASA's Space Shuttle. The equations of motion are derived, coded in Python, and solved using numerical integration. The results are then compared with existing literature to validate the simulations.

## Overview
The project focuses on understanding the various aspects involved in mission trajectory optimization, including deriving relevant equations, procuring accurate aerodynamic and structural databases, and solving the equations of motion to achieve accurate trajectory simulations. Trajectory optimization is a complex engineering problem, and this work contributes to the ongoing efforts in improving precision and accuracy in mission design.

## Features
6-DOF Equations of Motion: Implements a simplified, decoupled version of the 6-DOF equations of motion for atmospheric reentry.
Stardust SRC Simulation: Simulates the reentry trajectory of the Stardust Sample Return Capsule, which carried samples from comet Wild-2.
Space Shuttle Simulation: Simulates the reentry trajectory of NASA's Space Shuttle, a partially reusable low-Earth orbit spacecraft.
Numerical Integration: Solves the equations of motion using numerical integration (RK45 method) in Python.
Aerodynamic Data: Incorporates aerodynamic data (lift and drag coefficients) from existing literature for accurate simulations.
Result Validation: Compares the obtained results with published literature to validate the simulations.

## Authors
Sanjana Srivastava,
Nagachandra Nagaraj

## Acknowledgments
This project was carried out as part of the AE 598 Planetary Entry course at the University of Illinois at Urbana-Champaign, under the guidance of Professor Zachary Putnam.

## References
The project references the following publications:

Davis, P. (2021). Six Degree-of-Freedom Mission Planning for Reentry Trajectories. Air Force Institute of Technology.
Bollino, K., Oppenheimer, M., & Doman, D. (2006). Optimal Guidance Command Generation and Tracking for Reusable Launch Vehicle Reentry. AIAA Guidance, Navigation, and Control Conference and Exhibit.
Shoemaker, M., & Van der Ha, J. C. (2009). Trajectory estimation of the Hayabusa sample return capsule using optical sensors. 19th Workshop on Astrodynamics and Flight Mechanics.
Dilão, R., & Fonseca, J. (2016). Dynamic Guidance of Gliders in Planetary Atmospheres. Journal of Aerospace Engineering, 29(1), 04015012.
Desai, P., Mitcheltree, R., & Cheatwood, F. (1999). Entry Trajectory Issues for the Stardust Sample Return Capsule.
Ferrolho, H. (2020). Space Shuttle Reentry Trajectory. https://ferrolho.github.io/blog/2020-05-25/space-shuttle-reentry-trajectory

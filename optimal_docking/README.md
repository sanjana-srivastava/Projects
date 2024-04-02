# Spacecraft Optimal Docking Procedure Using LQR and Moving Targets
This project aims to design and implement an optimal control system for spacecraft docking maneuvers, utilizing a Linear Quadratic Regulator (LQR) controller. The controller handles the position and angle coordinates of the spacecraft, enabling docking with both stationary and moving targets.

## Overview
The project formulates the docking problem as an optimal control problem, where the objective is to minimize the time and control effort required for the spacecraft to dock with the target. The LQR controller adjusts the spacecraft's position and orientation by controlling the pitch, roll, and yaw angles, as well as the translational motion along the three axes.

Several scenarios are explored, including:

Docking with a stationary target
Docking with a non-rotating but predictably moving target
Docking with a rotating and predictably moving target
The report discusses the challenges involved, such as balancing optimality with safety constraints, ensuring adaptability to different target conditions, and managing computational complexity for moving targets.

## Features
LQR Controller: Implements a Linear Quadratic Regulator (LQR) controller for optimal spacecraft docking.
Multiple Scenarios: Explores docking scenarios with stationary, non-rotating moving, and rotating moving targets.
Adaptability: The controller adapts its approach based on the target's conditions and position relative to the spacecraft.
Constraints and Safety: Incorporates constraints and safety measures, such as velocity limits and collision avoidance, to ensure realistic and safe docking maneuvers.
Comparative Analysis: Compares the performance of the LQR controller with other methods, such as sliding mode control, for spacecraft docking.

## Authors
Derald Madson III,
Andrew Pelster,
Sanjana Srivastava

## Acknowledgments
This project was carried out as part of the course AE 504 Optimal Aerospace Systems under the guidance of Dr Negar Mehr at the University of Illinois at Urbana-Champaign.

## References
The project references the following publications:

Nandagopal, J. L., & Lethakumari, R. (2015). Optimal control of Spacecraft Docking System using integral LOR controller. In 2015 International Conference on Technological Advancements in Power and Energy (TAP Energy) (pp. 18-22). IEEE.
Guarnaccia, L., Bevilacqua, R., & Pastorelli, S. P. (2016). Suboptimal LQR-based spacecraft full motion control: Theory and experimentation. Acta Astronautica, 122, 114-136.
Yang, Y. (2012). Analytic LQR Design for Spacecraft Control System Based on Quaternion Model. Journal of Aerospace Engineering, 25(3), 448-453.

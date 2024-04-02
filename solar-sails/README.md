#Solar Sail Trajectory Optimization for Earth-Saturn Mission
This project aims to generate time-optimal trajectories for an interplanetary mission from Earth to Saturn, using solar sail propulsion. Solar sails harness solar radiation pressure for propulsion, offering a promising low-thrust alternative to traditional chemical propulsion for deep space exploration.

##Overview
The project formulates the trajectory optimization problem as a minimal time-of-flight problem, using the solar sail pitch angle as the control variable. The equations of motion, boundary conditions, cost function, Hamiltonian, and necessary conditions are derived systematically. Trajectories are optimized for scenarios with free and fixed final angular positions, and the effect of varying solar sail area is investigated.

##Features
Trajectory Optimization: Generates time-optimal trajectories for Earth-Saturn transfers using an ideal solar sail model.
Free and Fixed Final Position: Handles both free and fixed final angular position cases.
Solar Sail Area Study: Analyzes the impact of varying solar sail area on the minimum time of flight.
Validation: Results are validated against published literature for an Earth-Mars transfer.
Future Work: Discusses potential improvements, such as including perturbations (e.g., Jupiter's gravity) and non-ideal solar sail properties.

##Authors
Raman Singh
Robert Cazeau
Sanjana Srivastava

##Acknowledgments
This project was carried out as part of the AE508 Optimal Spacecraft Trajectories course at the University of Illinois at Urbana-Champaign, under the guidance of Dr. Robyn Woollands.

##References
The project references the following publications:

Zhang, M., Wang, G., & Sun, Y. (2010). Solar Sail Trajectory Optimization for Earth-Mars Mission. Chinese Control Conference.
Kim, M., & Hall, C. D. (2005). Symmetries in the optimal control of solar sail spacecraft. Celestial Mechanics and Dynamical Astronomy, 92(1-3), 25-45.
Jayaraman, T. S. (1980). Time-Optimal Orbit Transfer Trajectory for Solar Sail Spacecraft. Journal of Guidance, Control, and Dynamics, 3(6), 569-571.
Colasurdo, G., & Casalino, L. (2003). Optimal Control Law for Interplanetary Trajectories with Nonideal Solar Sail. Journal of Spacecraft and Rockets, 40(2), 260-265.

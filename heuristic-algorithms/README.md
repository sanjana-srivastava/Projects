# Heuristic Algorithms for Multidisciplinary Design Optimization
This project explores the application of heuristic optimization algorithms to solve a multidisciplinary design optimization (MDO) problem in aerospace engineering. The objective is to minimize the mass of an airship subject to various aerodynamic and mechanical constraints.

## Objective Function
The objective function for the airship mass minimization problem is formulated considering factors such as:

Geometry and shape of the airship
Lift and drag estimation
Power requirements and energy balance
Mass estimation of different components
The optimization is subject to constraints related to weight balance and energy balance to ensure stability of the airship.

## Heuristic Algorithms
Three heuristic optimization algorithms are implemented and compared:

Particle Swarm Optimization (PSO): Inspired by the social behavior of bird flocking or fish schooling, PSO optimizes a problem by iteratively improving candidate solutions based on a fitness function.
Genetic Algorithm (GA): GA is an evolutionary algorithm that mimics the process of natural selection. It evolves a population of candidate solutions using operators like selection, crossover, and mutation.
Simulated Annealing (SA): SA is a probabilistic technique inspired by the annealing process in metallurgy. It allows for a certain probability of accepting worse solutions during the search, enabling it to escape local optima.
Results and Analysis
Multiple iterations of each algorithm are performed with varying hyperparameters to analyze their performance and robustness. The results, including aerodynamic and power parameters, are tabulated and compared.

## Key findings:

PSO achieved consistent results across iterations with the lowest objective value of 235,903 kg.
GA showed more variability in results, with the lowest objective value of 260,086 kg.
SA faced issues with computation failures at higher initial values and performed the worst among the three algorithms.
The heuristic algorithms demonstrate potential for application in MDO problems, but careful selection of hyperparameters is crucial for their performance.

## Future Scope
Further fine-tuning of the algorithms, especially GA and SA, to improve their performance and robustness.
Testing the algorithms with different objective functions and engineering problems to gain more insights.
Exploring other heuristic and metaheuristic algorithms for MDO applications.
## Dependencies
MATLAB  
Global Optimization Toolbox  
Optimization Toolbox  
SRGTS Toolbox by Vienna  

## License
This project is licensed under the MIT License.

## Acknowledgements
Professor Manikandan M. for guidance and support throughout the project.
Manipal Institute of Technology for providing the necessary infrastructure and facilities.

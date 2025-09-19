A MATLAB-based tool for performance-based partitioning of large-scale district heating networks using distributed MPC.

This repository contains code for performing a branch-and-bound search of the partition of a large-scale district heating network (DHN). 
It also provides a 4-user case study to demonstrate the results.   

### Requirements

- MATLAB version: R2024b
- Toolboxes or dependencies:
  - Parallel Computing Toolbox
  - CasADi
  - ParforProgressbar

### To reproduce the main results:
The case study can be generated using 'generate_problem.m', and the results have been stored in sim_4user_params.mat for convenience.
The partitioning can be performed using 'partition_system.m', and the results have been stored in rslt_3_26.
The figures can be generated using 'process_results.m'
Additionally, the regression-based method for the case study can be explored using 'regression.m'.

Results have been compiled in 2 publications:
- Draft of paper on optimal partitioning: https://arxiv.org/abs/2507.02144
- Additional paper on regression accepted to MECC 2025.

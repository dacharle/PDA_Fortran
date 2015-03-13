# PDA_Fortran
The Population Dynamics Algorithm that was published by Charlebois et al. in 2011 in Communications in Computational Physics.

This parellel Fortran 90/95 program simulates the dynamics of a heterogeneous cell population using a by default the stochastic simulation algorithm (Gillespie, Journal of Physical Chemistry, 1977). Cell population size is kept constant using a constant-number Monte Carlo method (Mantzaris et al., Journal of Theoretical Biology, 2006). This version incorporates mutation at cell division.

Notes: 'CELL_PASSED.f90' is the main caller program which seeds the random number generator and calls 'SSA_Parallel.f90' in order to simulate the dynamics of the system and output the results at a user defined sampling interval. 'SSA_Parallel.f90' calls 'reactions_Parallel.f90' and 'results_Parallel.f90' as required. Global variables are stored in 'globals_Parallel.f90'. To implement parallelization, this program must be compiled using the '-fopenmp' flag.

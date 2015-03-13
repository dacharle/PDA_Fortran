!This module contains the global parameters.

MODULE globals_Parallel

real, parameter :: sample_intervl=90 !population sampling interval *choose sample_intervl and t_restore to be the same*
real, parameter :: t_restore=90 !cell population restore time 
real, parameter :: t_end=900 !simulation end time 
integer, parameter :: num_rand=100000 !number of random numbers
real, parameter :: log2=log10(2.), pi=3.14159265 !constants 
real, parameter :: init_cell_vol=1.0 !initial/following division cell volume 
real, parameter :: t0_div=90 !cell division time in absence of selection 
real, parameter :: D_init=1 !initial gene copy number on chromosome
integer, parameter :: n_threads=2 !number of threads
integer, parameter :: drug_env=0 !apply drug 0 or 1 (no or yes, respectively)
integer, parameter :: n_runs=10 !number of realizations
real, parameter :: K2_gene=150, n2_gene=4, K3_gene=150, n3_gene=4 !Hill parameters for gene expression
real, parameter :: K1=1000,n1=2,K2=1000,n2=2,K3=500,n3=2 !Hill parameters for fitness function (K is prop to drug conc)
real, parameter :: kP1=10,dP1=0.01,kP2_basal=0.1,kP2_act=10,dP2=0.01,kP3_basal=0.01,kP3_act=0.1,dP3=0.0001 !gene expression model parameters
real, parameter :: mutation_threshold=10**7 !threshold for a mutation to occur at cell division
integer, parameter :: mutation=0 !apply mutation 0 or 1 (no or yes, respectively)

END MODULE globals_Parallel

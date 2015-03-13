!Daniel Charlebois - Fall 2014 - Population Dynamics Algorithm (Charlebois et al. Communications in Computational Physics, 2011) 
!This Fortran 90/95 program simulates the dynamics of a heterogeneous cell population
!using a stochastic simulation algorithm (default SSA, Gillespie direct method, Journal of Physical Chemistry, 1977).
!Cell population size is kept constant using a constant-number Monte Carlo method (Mantzaris et al., Journal of Theoretical Biology, 2006).
!'CELL_PASSED.f90' is the main caller program which seeds the random number generator and calls the SSA (SSA_Parallel.f90) in order 
!to simulate the dynamics of the system and output the results at a user defined sampling interval. 
!'SSA_Parallel.f90' calls 'reactions_Parallel.f90' and 'results_Parallel.f90' as required. Global variables are stored in 'globals_Parallel.f90'.
!gfortran -fopenmp -o PDA_Parallel CELL_PASSED_Parallel.f95 SSA_Parallel.f95 globals_Parallel.f95 reactions_Parallel.f95 results_Parallel.f95

!***************************************************
!!!MAIN CALLER PROGRAM!!!                 
                                       
PROGRAM CELL_PASSED_Parallel  
use globals_Parallel                   
use SSA_Parallel                           
implicit none    
integer, allocatable :: iseed(:)  
integer :: idate(8), isize, k, stat_cnt
integer :: nc, ncmc !initial/real-time number of cells,fixed number of cells for CNMC          
real :: t_cumm !simulation start time                
real, dimension(n_runs,nint(t_end/t_restore)+1) :: mr1, mr2 !master results arrays
real :: tot_W,tot_frac,var_cnt_W,var_cnt_frac 
real, dimension(nint(t_end/t_restore)+1) :: avg_W, avg_frac
integer :: i,j
   
   open(3,file='data.dat')  
                                 
   call DATE_AND_TIME(VALUES=idate)        
   print *, "Simulation start date/time = ", idate

    mr1=0; mr2=0; stat_cnt=0
     
    do  k=1,n_runs      
    	t_cumm=0; nc=1000; ncmc=1000
		!set initial number of cells,fixed number of cells,simulation start time   
  		stat_cnt=stat_cnt+1                 
		!obtain random seed from date_and_time intrinsic                    
		call random_seed(SIZE=isize)	
		allocate(iseed(isize))
		call random_seed(GET=iseed)
		iseed = iseed * (idate(8)-500) !idate(8) contains millisecond
		call random_seed(PUT=iseed)
		deallocate(iseed)
	
		call main(nc,ncmc,t_cumm,mr1,mr2,stat_cnt)
		
		print*, 'run number =', k
	end do
	
	!obtain statistics over the realizations
	!calcuate mean
	tot_W=0; tot_frac=0; 
	do j=1,nint(t_end/t_restore)+1
		tot_W=0; tot_frac=0
		do i=1,n_runs
			tot_W=tot_W+mr1(i,j)
			tot_frac=tot_frac+mr2(i,j)
		end do
		avg_W(j)=tot_W/n_runs
		avg_frac(j)=tot_frac/n_runs
	end do	
	!calculate variance
	var_cnt_W=0; var_cnt_frac=0 
	do j=1,nint(t_end/t_restore)+1
		var_cnt_W=0; var_cnt_frac=0
		do i=1,n_runs
			var_cnt_W=var_cnt_W+((mr1(i,j)-avg_W(j))**2)
			var_cnt_frac=var_cnt_frac+((mr2(i,j)-avg_frac(j))**2)
		end do
		write(3,*) avg_W(j), avg_frac(j), sqrt(var_cnt_W/n_runs), sqrt(var_cnt_frac/n_runs)
	end do	
	
	close(3)

	call DATE_AND_TIME(VALUES=idate)
	print *, "Simulation end date/time = ", idate          
   
END PROGRAM CELL_PASSED_Parallel         
!***************************************************



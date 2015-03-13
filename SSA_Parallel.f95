!This module contains the SSA, which calls the module 'reactions_Parallel.f90' in order to simulate gene expression
!for a user defined genetic network and calculates the reaction propensities and determines the next reaction via subroutines 'propensities'
!and 'next_reaction', respectively. Module 'results_Parallel.f90' is called every 'sample_intervl' in order to calculate cell statistics 
!as well as output results.

MODULE SSA_Parallel
use reactions_Parallel
use globals_Parallel  
use results_Parallel
implicit none
contains 	

!**********************************
!!!MAIN ALGORITHM!!!
SUBROUTINE main(nc,ncmc,t_cumm,mr1,mr2,stat_cnt) 
use globals_Parallel  
real, dimension(n_threads,num_rand) :: rand_array1, rand_array2 !random number arrays for r1 and r2 
integer, dimension(n_threads) :: rand_cnt_array
real, dimension(ncmc,15) :: u !array to store data for mother cells
!KEY: [1:t_sample,2:M,3:P1,4:P2,5:P3,6:k,7:mu,8:cell_vol,9:time since last division,10:t_cumm,11:cell_age,12:mother_ID,13:fitness,14:t_div/fitness,15:cell_status]
real, dimension(15) :: u_temp_mothers, u_temp_daughters 
real, dimension(2,15) :: u_mothers_temp, u_daughters_temp
real, dimension(ncmc,15) :: u_daughters !array to store data for cells
real, dimension(1,4) :: s !array to store variables: A,R,M,P
real, dimension(1,9) :: a !array to store reaction propensities
real, dimension(nc)  :: t_last_div_array,cell_vol_array,init_N !array to store initial last division times,cell volumes,and rand array to determine initial molecular distribution
real :: t_sample=0,t_cumm,dt=0,t_marker !sample interval,simulation time,next reaction time
real :: M,P1,P2,P3,cell_age=0,mother_ID=0 !initial conditions
real :: a0=0,amu=0 !sum of reaction propensities and value holder for next reaction
integer :: k=0,mu=0,i=0,nc,nc_cnt,ncmc !cell I.D.,next reaction marker
integer :: TID, OMP_GET_THREAD_NUM !thread, get current thread (for OpenMP parallelization)
integer :: mother_sim,daughter_sim !indicates to SSA whether to simulate mother or daughter cells
integer :: daughter_cnt,daughter_cnt2,daughter_cnt_hld,d_count,d_count2,d_count3 !counters for number of daughter cells generated from mothers in the current population
!Constant Number Algorithm (CNA)
real, dimension(2*ncmc) :: CNA_array 
integer, dimension(ncmc) :: CNA_rand_index_array
real :: CNA_rand,w !fitness
integer :: ii,swap,t_restore_cnt,CNA_rand_index,CNA_rand_index_cnt,CNA_cnt1,CNA_cnt2,j !misc. for implementing CNA 
integer :: E1=1,E2=0,cell_status !environments,cell state marker
integer, dimension(1) :: oldest_idx !position of oldest daughter cell for CNMC replacement
!other
real, dimension(n_runs,nint(t_end/t_restore)+1) :: mr1,mr2
integer :: stat_cnt,stat_marker,stat_marker2

	!open(1,file='time_series_Parallel.dat') !file to store time series for selected cell(s)
	open(2,file='pop_stats_Parallel.dat')   !file to store population statistics

	!initialize arrays
	rand_array1=0; rand_array2=0; rand_cnt_array=1; a=0; s=0; u=0; u_temp_mothers=0; u_daughters=0; u_temp_daughters=0
	call random_number(rand_array1); call random_number(rand_array2); rand_cnt_array=1; call random_number(CNA_array)

	!initialize logicals
	mother_sim=1; daughter_sim=0

	!intialize misc.
	t_restore_cnt=1; CNA_rand_index_cnt=0; stat_marker=0; stat_marker2=0

	!initialize time since last division for population
	!uniform distribution 
	call random_number(t_last_div_array) 
	t_last_div_array=t_last_div_array*t0_div  
	!!synchronized population  
	!t_last_div_array=0
	!!initial cell volume
	do k=1,nc
		cell_vol_array(k)=init_cell_vol*EXP(log(2.)*(t_last_div_array(k)/t0_div)) !initial fitness here set to 1
	end do
	
!	!determine initial M and P concentrations and calculate initial fitness
!	call random_number(init_N)
!	if (E1==1) then
!		do k=1,nc
!			if (init_N(k)<=0.0) then
!				u(k,2)=10; u(k,3)=175  !Unfit in E1
!			else
!				u(k,2)=20; u(k,3)=350 !Fit in E1
!			end if
!			!intial cell volume and fitness
!			cell_vol_array(k)=init_cell_vol !initial cell volume
!			u(k,12)=((u(k,3)/cell_vol_array(k))**n)/(((u(k,3)/cell_vol_array(k))**n)+(K1**n)) !cell E1 fitness
!		end do
!	elseif (E2==1) then
!		do k=1,nc
!			if (init_N(k)<=0.0) then
!				u(k,2)=20; u(k,3)=350 !Unfit in E2 
!			else
!				u(k,2)=10; u(k,3)=175 !Fit in E2
!			end if
!			!intial cell volume and fitness
!			cell_vol_array(k)=init_cell_vol !initial cell volume
!			u(k,12)=(K2**n)/((K2**n)+((u(k,3)/cell_vol_array(k))**n)) !cell E2 fitness
!		end do
!	end if

	!initial values for cells
	do k=1,nc
		u(k,1)=0;u(k,2)=0;
		u(k,3)=kP1/dP1
		u(k,4)=kP2_act/dP2 
		u(k,5)=kP3_act/dP3
		u(k,6)=k;u(k,7)=0;u(k,8)=cell_vol_array(k);u(k,9)=t_last_div_array(k)
		u(k,10)=t_cumm;u(k,11)=0;u(k,12)=0
		u(k,13)=((K1**n1)/((K1**n1)+(u(1,3)**n1)))*((K2**n2)/((K2**n2)+(u(1,4)**n2)))*&
						((u(1,5)**n3)/((u(1,5)**n3)+(K3**n3)))
		u(k,14)=t0_div/u(k,13);u(k,15)=1 
	end do
	if (mutation==1) then
		u(1,2)=1
	end if

	!print*, "P1 = ", u(1,3), "P2 = ", u(1,4), "P3 = ", u(1,5) 

	d_count=0; daughter_cnt=0	

	!obtain initial values for time series and statistics for cell(s) 
	stat_marker2=stat_marker2+1
	!call time_series(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,stat_marker,stat_marker2)
	call pop_stats(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,E1,E2,mr1,mr2,stat_cnt,stat_marker,stat_marker2)		

	t_sample=sample_intervl
	do while (t_cumm < t_end) !*make sure t_sample is set such that t_cumm will always be < t_sample*
		t_marker=t_cumm !t_marker is time at begining of this t_sample (so that each cell starts at the beginning of interval each iteration)
		if (t_cumm < t_sample) then  
!			!changing environment
!			if (t_sample <= t_end/2) then
!				E1=1; E2=0
!			else 
!				E1=0; E2=1
!			end if

	!PARALLEL REGION: simulate dynamics of mother cells
	!$OMP PARALLEL PRIVATE(w,cell_status,M,P1,P2,P3,k,t_cumm,dt,s,a,a0,amu,mu,i,u_temp_mothers,u_temp_daughters) NUM_THREADS(n_threads) 
	!$OMP DO 
	do k=1,nc
		u_temp_mothers = u(k,1:15)
		call next_reaction(t_cumm,t_sample,t_marker,M,P1,P2,P3,rand_array1,rand_array2,k,nc,ncmc,rand_cnt_array,dt,s,a,a0,amu,mu,i, & 
		u_temp_mothers,u_temp_daughters,u_daughters,daughter_cnt,mother_sim,daughter_sim,d_count,E1,E2,w,cell_status)
		u(k,1:15) = u_temp_mothers
	end do
	!$OMP END DO
	!$OMP END PARALLEL

	!set logicals
	mother_sim=0;daughter_sim=1 

	!PARALLEL REGION: simulate dynamics of daughter cells
	if (daughter_cnt >= 1) then
	!$OMP PARALLEL PRIVATE(w,cell_status,M,P1,P2,P3,k,t_cumm,dt,s,a,a0,amu,mu,i,u_temp_mothers,u_temp_daughters) NUM_THREADS(n_threads) 
	!$OMP DO
		do d_count=1,daughter_cnt
			u_temp_daughters = u_daughters(d_count,1:15)
			call next_reaction(t_cumm,t_sample,t_marker,M,P1,P2,P3,rand_array1,rand_array2,k,nc,ncmc,rand_cnt_array,dt,s,a,a0,amu,mu,i, &
			u_temp_mothers,u_temp_daughters,u_daughters,daughter_cnt,mother_sim,daughter_sim,d_count,E1,E2,w,cell_status)		
			u_daughters(d_count,1:15) = u_temp_daughters
		end do
	!$OMP END DO
	!$OMP END PARALLEL
	end if
	
			!lethal drug application
			if (drug_env==1) then
				if (t_sample >= 1000) then !some t_sample
					stat_marker=1
					!call pop_stats(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,E1,E2,master_results,stat_cnt,stat_marker,stat_marker2)	
					nc_cnt=nc
					do k=1,nc_cnt
						!fitness for fixed volume model
						if (u(k,5) < K1) then
							u(k,1:15)=0; nc=nc-1
						end if						

						!fitness for cell growth model
						!u(k,13)=((u(k,5)/u(k,8))**n)/(((u(k,5)/u(k,8))**n)+(K1**n)) !new mother cells' fitness
						!if (u(k,13) < 0.5) then
						!	u(k,1:15)=0; nc=nc-1
						!end if
					end do
					daughter_cnt2=daughter_cnt
					do k=1,daughter_cnt2
						!fitness for fixed volume model
						if (u_daughters(k,5) < K2) then
							u_daughters(k,1:15)=0; daughter_cnt=daughter_cnt-1
						end if

						!fitness with cell growth model
						!u_daughters(k,13) = ((u_daughters(k,5)/u_daughters(k,8))**n)/(((u_daughters(k,5)/u_daughters(k,8))**n)+(K1**n)) !new daughter cells' fitness
						!if (u_daughters(k,13) < 0.5) then
						!	u_daughters(k,1:15)=0; daughter_cnt=daughter_cnt-1
						!end if
					end do
					!sort mother cell array (to eliminate any empty rows resulting from selective pressure)
					u_mothers_temp=0
					if (nc_cnt >= 2) then
						ii=1
						do while (ii <= nc_cnt-1)
							CNA_cnt1=u(ii,11)
							CNA_cnt2=u(ii+1,11)
							if (CNA_cnt2 > CNA_cnt1) then
								u_mothers_temp(1,1:15)=u(ii,1:15)	
								u_mothers_temp(2,1:15)=u(ii+1,1:15)
								u(ii,1:15)=u_mothers_temp(2,1:15)
								u(ii+1,1:15)=u_mothers_temp(1,1:15)		
								ii=1
							else
								ii=ii+1
							end if
						end do		
					end if			
					!sort daughter cell array
					u_daughters_temp=0
					if (daughter_cnt2 >= 2) then
						ii=1
						do while (ii <= daughter_cnt2-1)
							CNA_cnt1=u_daughters(ii,11)
							CNA_cnt2=u_daughters(ii+1,11)
							if (CNA_cnt2 > CNA_cnt1) then
								u_daughters_temp(1,1:15)=u_daughters(ii,1:15)	
								u_daughters_temp(2,1:15)=u_daughters(ii+1,1:15)
								u_daughters(ii,1:15)=u_daughters_temp(2,1:15)
								u_daughters(ii+1,1:15)=u_daughters_temp(1,1:15)		
								ii=1
							else
								ii=ii+1
							end if
						end do		
					end if		
				end if
			end if
			stat_marker=0

			!reset logicals
			mother_sim=1; daughter_sim=0

			t_cumm=t_sample	
			t_sample=t_sample+sample_intervl

		end if

		!CNMC 
		if (daughter_cnt >= 1) then
			if (t_restore*t_restore_cnt <= t_cumm) then

				d_count2=0

				!fill mother array with daughter cells until nc_fixed is reached (bug for d_count2 fixed)
				if (nc < ncmc) then
				    do while (d_count2 < daughter_cnt .and. nc+1 <= ncmc)
				    	d_count2=d_count2+1
					oldest_idx = MAXLOC(u_daughters(:,11))
					u(nc+1,1:15)=u_daughters(oldest_idx(1),1:15)
					u_daughters(oldest_idx(1),11)=0	
				    	nc=nc+1
				    end do
				end if

				CNA_array=0; call random_number(CNA_array); CNA_rand_index_cnt=0 
				do while (d_count2 < daughter_cnt)
					d_count2=d_count2+1
					!replace random cells in mother cell array with cells from daughter array (oldest daughter cells being inserted first)
					CNA_rand_index_cnt=CNA_rand_index_cnt+1	
		    			CNA_rand_index_array(d_count2)=nint(CNA_array(CNA_rand_index_cnt)*nc)	
					!ensure that random index array is non-zero (i.e. there must be a mother cell to replace!)  
					do while (CNA_rand_index_array(d_count2)==0) 
						CNA_rand_index_cnt=CNA_rand_index_cnt+1
						CNA_rand_index_array(d_count2)=nint(CNA_array(CNA_rand_index_cnt)*nc)
					end do	
					!do while (CNA_rand_index_array(d_count2)==1) !code to track first cell (so that it is never replaced) 
					!	CNA_rand_index_cnt=CNA_rand_index_cnt+1
					!	CNA_rand_index_array(d_count2)=nint(CNA_array(CNA_rand_index_cnt)*nc)
					!end do
					oldest_idx = MAXLOC(u_daughters(:,11))
					u(CNA_rand_index_array(d_count2),1:15)=u_daughters(oldest_idx(1),1:15)
					u_daughters(oldest_idx(1),11) = 0
				end do
				
				daughter_cnt_hld=daughter_cnt
				daughter_cnt=0; u_daughters=0
				t_restore_cnt=t_restore_cnt+1
			end if
		end if
		
		!if (t_cumm /= t_end/2) then
			!obtain time series for cell(s) 
			!call time_series(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,stat_marker,stat_marker2)
			!obtain population statistics for cells
			stat_marker2=stat_marker2+1		
			call pop_stats(u,u_daughters,nc,ncmc,daughter_cnt_hld,t_cumm,t_sample,E1,E2,mr1,mr2,stat_cnt,stat_marker,stat_marker2)		
		!end if
	end do 

	!close(1); 
	close(2); 

END SUBROUTINE main
!**********************************
   
!**********************************
SUBROUTINE next_reaction(t_cumm,t_sample,t_marker,M,P1,P2,P3,rand_array1,rand_array2,k,nc,ncmc,rand_cnt_array,dt,s,a,a0,amu,mu,i, &
u_temp_mothers,u_temp_daughters,u_daughters,daughter_cnt,mother_sim,daughter_sim,d_count,E1,E2,w,cell_status)
!use reactions_Parallel
use globals_Parallel  
!use globals_Parallel
implicit none
real, dimension(n_threads,num_rand) :: rand_array1, rand_array2
integer, dimension(n_threads) :: rand_cnt_array
real, dimension(ncmc,15) :: u_daughters
real, dimension(15) :: u_temp_mothers, u_temp_daughters 
real, dimension(1,4) :: s 
real, dimension(1,6) :: a 
real :: dt,t_cumm,t_sample,t_marker
real :: M,P1,P2,P3,a0,amu,w 
integer :: k,i,nc,ncmc,mu,idate(8), TID, OMP_GET_THREAD_NUM 
integer :: mother_sim,daughter_sim,daughter_cnt,d_count 
integer :: E1,E2,cell_status
integer, parameter :: cell_growth=1, cell_div=1 !cell growth/division ON (1) or OFF (0)

	!Mothers
	if (mother_sim==1) then
	 	if (u_temp_mothers(15) == 1) then
	 		cell_status=1 !cell is alive
			M=u_temp_mothers(2)
			P1=u_temp_mothers(3)
			P2=u_temp_mothers(4)	
			P3=u_temp_mothers(5) 
			w=u_temp_mothers(13) 
			t_cumm=t_marker
		else
			cell_status=0 !cell is dead
		end if
	end if

	!Daughters
	if (daughter_sim==1) then
		if (u_temp_daughters(15) == 1) then
			cell_status=1
			M=u_temp_daughters(2)
			P1=u_temp_daughters(3)
			P2=u_temp_daughters(4)	
			P3=u_temp_daughters(5) 
			w=u_temp_daughters(13)
			t_cumm=u_temp_daughters(10)
		else
			cell_status=0	
		end if	
	end if

if (cell_status==1) then
 
	!TID = 0
	TID = OMP_GET_THREAD_NUM()  
	!print*, "TID = ", TID

	do while (t_cumm < t_sample) 
		s(1,1)=M
		s(1,2)=P1
		s(1,3)=P2
		s(1,4)=P3

		call propensities(s,a,E1,E2,w)

		a0 = sum(a)

		!check and replenish random number arrays as required
		if (rand_cnt_array(TID+1)+6 >= num_rand) then	
			!print *, 'Replenishing rand_arrays', k
			call random_number(rand_array1(TID+1,num_rand))	
			!call random_number(rand_array2(TID+1,num_rand))			
			rand_cnt_array(TID+1)=1	
		end if

		!calculate time of next reaction
		dt=(log(1./rand_array1(TID+1,rand_cnt_array(TID+1)))/a0)
		rand_cnt_array(TID+1)=rand_cnt_array(TID+1)+1

		!update time and time since last division
		t_cumm=t_cumm+dt
		
		if (mother_sim==1) then
			u_temp_mothers(9)=u_temp_mothers(9)+dt
			u_temp_mothers(10)=u_temp_mothers(10)+dt
			u_temp_mothers(11)=u_temp_mothers(11)+dt
		elseif (daughter_sim==1) then
			u_temp_daughters(9)=u_temp_daughters(9)+dt
			u_temp_daughters(10)=u_temp_daughters(10)+dt
			u_temp_daughters(11)=u_temp_daughters(11)+dt
		end if

	if (cell_growth==1) then
		if (mother_sim==1) then		
			!!fix cell volume (for time based cell division)
			u_temp_mothers(8)=init_cell_vol
			!!calculate cell volume
			!u_temp_mothers(8)=init_cell_vol*EXP(log(2.)*(u_temp_mothers(9)/u_temp_mothers(14))) !exponential growth 
			!u_temp_mothers(8)=init_cell_vol*(1+(u_temp_mothers(9)/u_temp_mothers(14))) !linear growth 
			if (cell_div==1) then
				!if (u_temp_mothers(8) >= 2*init_cell_vol) then !division based on cell volume
				!if (u_temp_mothers(9) >= t0_div) then !division based on division time (t_div)
				if (u_temp_mothers(9) >= u_temp_mothers(14)) then
					daughter_cnt=daughter_cnt+1
	
					!!P1 distribution at cell division 	
					!if (mod(u_temp_mothers(3),2.)==0) then
					!	u_temp_mothers(3) = real(nint(u_temp_mothers(3)/2)); P1=u_temp_mothers(3)
					!	u_daughters(daughter_cnt,3)=u_temp_mothers(3)
					!else
					!	if (nint(rand_array1(TID+1,rand_cnt_array(TID+1)))==0) then
					!		u_temp_mothers(3)=nint(0.5*u_temp_mothers(3)); P1=u_temp_mothers(3)
					!		u_daughters(daughter_cnt,3)=P1-1
					!	else
					!		u_temp_mothers(3)=nint(0.5*u_temp_mothers(3))-1; P1=u_temp_mothers(3)
					!		u_daughters(daughter_cnt,3)=P1+1
					!	end if	
					!	rand_cnt_array(TID+1)=rand_cnt_array(TID+1)+1
					!end if
					!!P2 distribution at cell division 	
					!if (mod(u_temp_mothers(4),2.)==0) then
					!	u_temp_mothers(4) = real(nint(u_temp_mothers(4)/2)); P2=u_temp_mothers(4)
					!	u_daughters(daughter_cnt,4)=u_temp_mothers(4)
					!else
					!	if (nint(rand_array1(TID+1,rand_cnt_array(TID+1)))==0) then
					!		u_temp_mothers(4)=nint(0.5*u_temp_mothers(4)); P2=u_temp_mothers(4)
					!		u_daughters(daughter_cnt,4)=P2-1
					!	else
					!		u_temp_mothers(4)=nint(0.5*u_temp_mothers(4))-1; P2=u_temp_mothers(4)
					!		u_daughters(daughter_cnt,4)=P2+1
					!	end if	
					!	rand_cnt_array(TID+1)=rand_cnt_array(TID+1)+1
					!end if	
					!!P3 distribution at cell division 	
					!if (mod(u_temp_mothers(5),2.)==0) then
					!	u_temp_mothers(5) = real(nint(u_temp_mothers(5)/2)); P3=u_temp_mothers(5)
					!	u_daughters(daughter_cnt,5)=u_temp_mothers(5)
					!else
					!	if (nint(rand_array1(TID+1,rand_cnt_array(TID+1)))==0) then
					!		u_temp_mothers(5)=nint(0.5*u_temp_mothers(5)); P3=u_temp_mothers(5)
					!		u_daughters(daughter_cnt,5)=P3-1
					!	else
					!		u_temp_mothers(5)=nint(0.5*u_temp_mothers(5))-1; P3=u_temp_mothers(5)
					!		u_daughters(daughter_cnt,5)=P3+1
					!	end if	
					!	rand_cnt_array(TID+1)=rand_cnt_array(TID+1)+1
					!end if	
					
					if (mutation==1) then
						if (u_temp_mothers(2)==0) then
							if (rand_array1(TID+1,rand_cnt_array(TID+1))*mutation_threshold<=1) then
								rand_cnt_array(TID+1)=rand_cnt_array(TID+1)
								if (rand_array1(TID+1,rand_cnt_array(TID+1))*mutation_threshold<=1) then
								u_temp_mothers(2)=1; M=1
								end if
							end if
							rand_cnt_array(TID+1)=rand_cnt_array(TID+1)
						end if
					end if
				  
					u_temp_mothers(3)=P1;
					u_temp_mothers(4)=P2;
					u_temp_mothers(5)=P3;  

					u_temp_mothers(8)=init_cell_vol
					u_temp_mothers(9)=0

					!FOR ASYMMETRIC DIVISION RECALCULATE FITNESS
					!if (E1==1) then !update E1 fitness    
	 				!     u_temp_mothers(13) = 1 ! (K1**n)/((K1**n)+((u_temp_mothers(5)/u_temp_mothers(8))**n))
					!elseif (E2==1) then !update E2 fitness
				             !u_temp_mothers(13) = ((u_temp_mothers(5)/u_temp_mothers(8))**n)/(((u_temp_mothers(5)/u_temp_mothers(8))**n)+(K1**n)) 	
					!end if
					!u_temp_mothers(14)= t0_div/u_temp_mothers(13) !calc. new t_div  
	
					u_daughters(daughter_cnt,1)=u_temp_mothers(1)
					u_daughters(daughter_cnt,2)=u_temp_mothers(2)
					u_daughters(daughter_cnt,3)=u_temp_mothers(3)
					u_daughters(daughter_cnt,4)=u_temp_mothers(4)
					u_daughters(daughter_cnt,5)=u_temp_mothers(5)	
					u_daughters(daughter_cnt,6)=daughter_cnt !correct daughter cell number (may not be same as 'k')
					u_daughters(daughter_cnt,7)=u_temp_mothers(7)
					u_daughters(daughter_cnt,8)=u_temp_mothers(8)
					u_daughters(daughter_cnt,9)=0
					u_daughters(daughter_cnt,10)=u_temp_mothers(10) !t_cumm
					u_daughters(daughter_cnt,11)=0
					u_daughters(daughter_cnt,12)=k
					u_daughters(daughter_cnt,13)=u_temp_mothers(13)
					u_daughters(daughter_cnt,14)=u_temp_mothers(14)
					u_daughters(daughter_cnt,15)=1
					!snap shot at moment of division
					!print*, 'mother', u_temp_mothers; print*, 'daughter', u_daughters(daughter_cnt,1:15) 
		end if  
	end if
		elseif (daughter_sim==1) then
				!calculate cell volume
				u_temp_daughters(8)=init_cell_vol
				!u_temp_daughters(8)=init_cell_vol*(1+(u_temp_daughters(9)/u_temp_daughters(14))) linear growth (Swain)
				!u_temp_daughters(8)=init_cell_vol*EXP(log(2.)*u_temp_daughters(9)/u_temp_daughters(14)) !exponential growth (Kaern)
                            
				!if (u_daughters(d_count,8) >= 2*init_cell_vol) then 
				if (u_temp_daughters(9) >= u_temp_daughters(14)) then !division based on fixed division time (t_div)
					write(1,*) 'DAUGHTER CELL HAS DIVIDED! SET t_restore LESS THAN SHORTEST CELL DIVISION TIME!!!' 
				end if       
		end if	
		end if

		!calculate which reaction occurs next (if outside of dimension of rand_array2 then replenish with random numbers)
		i=1; mu=0; amu=0 !set i=1 for each cycle			
		do while (amu < rand_array1(TID+1,rand_cnt_array(TID+1))*a0)		       		
		        mu=mu+1
		    	do while (i <= mu)
			    	amu=amu+a(1,i)                     
			    	i=i+1
		    	end do
		end do
		rand_cnt_array(TID+1)=rand_cnt_array(TID+1)
		
		!carry out reaction
		call rxns(mu,M,P1,P2,P3) !call rxns subroutine in 'reactions_Parallel' module		
		     
		!save cell specific data to u_temp
		if (mother_sim==1) then
			u_temp_mothers(1) = t_sample
			u_temp_mothers(2) = M
			u_temp_mothers(3) = P1
			u_temp_mothers(4) = P2
			u_temp_mothers(5) = P3
			u_temp_mothers(6) = k
			u_temp_mothers(7) = mu 		
			!if (E1==1) then !update E1 fitness 			
			!	u_temp_mothers(13) = 1 ! (K1**n)/((K1**n)+((u_temp_mothers(5)/u_temp_mothers(8))**n))
			!elseif (E2==1) then !update E2 fitness
			!	u_temp_mothers(13) = ((u_temp_mothers(5)/u_temp_mothers(8))**n)/(((u_temp_mothers(5)/u_temp_mothers(8))**n)+(K2**n)) 
			!end if
			u_temp_mothers(13) = ((K1**n1)/((K1**n1)+(u_temp_mothers(3)**n1)))*((K2**n2)/((K2**n2)+(u_temp_mothers(4)**n2)))*&
						((u_temp_mothers(5)**n3)/((u_temp_mothers(5)**n3)+(K3**n3)))
			u_temp_mothers(14) = t0_div/u_temp_mothers(13)
		end if	

		!save data to u_daughters
		if (daughter_sim==1) then
			u_temp_daughters(1) = t_sample
			u_temp_daughters(2) = M
			u_temp_daughters(3) = P1
			u_temp_daughters(4) = P2
			u_temp_daughters(5) = P3
			u_temp_daughters(6) = d_count
			u_temp_daughters(7) = mu  
			!if (E1==1) then !update E1 fitness				 
			!	u_temp_daughters(13) = 1 ! (K1**n)/((K1**n)+((u_temp_daughters(5)/u_temp_daughters(8))**n))
			!elseif (E2==1) then !update E2 fitness
			!	u_temp_daughters(13) = ((u_temp_daughters(5)/u_temp_daughters(8))**n)/(((u_temp_daughters(5)/u_temp_daughters(8))**n)+(K2**n)) 
			!end if
			u_temp_daughters(13) = ((K1**n1)/((K1**n1)+(u_temp_daughters(3)**n1)))*((K2**n2)/((K2**n2)+(u_temp_daughters(4)**n2)))*&
				((u_temp_daughters(5)**n3)/((u_temp_daughters(5)**n3)+(K3**n3)))
			u_temp_daughters(14) = t0_div/u_temp_daughters(13)
		end if

	end do 

end if

END SUBROUTINE next_reaction
!**********************************

!**********************************
SUBROUTINE propensities(s,a,E1,E2,w)
use globals_Parallel
implicit none
real, dimension(1,4) :: s !array to store variables
real, dimension(1,6) :: a !propensities vector
real ::  M,P1,P2,P3,w
integer :: E1, E2 

	!conversion into Non and Noff from 's' array
	M=s(1,1)
	P1=s(1,2)
	P2=s(1,3)
	P3=s(1,4) 
	
	if (M==0) then
		a(1,1)=kP1
		a(1,2)=dP1*P1
		a(1,3)=kP2_basal+(kP2_act-kP2_basal)*(P1**n2_gene/((K2_gene**n2_gene)+(P1**n2_gene)))
		a(1,4)=dP2*P2
		a(1,5)=kP3_basal+(kP3_act-kP3_basal)*(P2**n3_gene/((K3_gene**n3_gene)+(P2**n3_gene)))
		a(1,6)=dP3*P3
	elseif (M==1) then
		a(1,1)=0
		a(1,2)=dP1*P1
		a(1,3)=kP2_basal+(kP2_act-kP2_basal)*(P1**n2_gene/((K2_gene**n2_gene)+(P1**n2_gene)))
		a(1,4)=dP2*P2
		a(1,5)=kP3_basal+(kP3_act-kP3_basal)*(P2**n3_gene/((K3_gene**n3_gene)+(P2**n3_gene)))
		a(1,6)=dP3*P3
	end if

END SUBROUTINE propensities
!**********************************

END MODULE SSA_Parallel





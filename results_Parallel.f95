!Obtain time-series and calculate statistics.

MODULE results_Parallel
implicit none
contains 

!**********************************
SUBROUTINE time_series(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,stat_marker,stat_marker2)
use globals_Parallel
implicit none
real, dimension(ncmc,15) :: u
real, dimension(ncmc,15) :: u_daughters 
real :: t_cumm,t_sample
integer :: j,k,nc,ncmc,daughter_cnt,stat_marker,stat_marker2
	
	!obtain time series for mother cells
	!write(1,*) 'Mother Cells'
	do k=1,nc
		if (k==1) then !obtain values for a single cell (comment to obtain values for all cells)
		!if (k==etc.or.3.or.2.or.1) then !obtain values for select cells 
			write(1,*) (u(k,1)), (u(k,2)), (u(k,3)), (u(k,4)), (u(k,5)), (u(k,6)), (u(k,7)), &
			(u(k,8)), (u(k,9)), (u(k,10)), (u(k,11)), (u(k,12)), (u(k,13)), (u(k,14)), (u(k,15))    
		end if
	end do

	!obtain time series for daughter cells
	if (daughter_cnt >= 1) then
		!write(1,*) 'Daughter Cells'
		do j=1,daughter_cnt
			if (j==1) then !obtain values for a single cell (comment to obtain values for all cells)
			!if (j==etc.or.3.or.2.or.1) then !obtain values for select cells 
				write(1,*) (u_daughters(j,1)), (u_daughters(j,2)), (u_daughters(j,3)), (u_daughters(j,4)), &
				(u_daughters(j,5)), (u_daughters(j,6)), (u_daughters(j,7)), (u_daughters(j,8)), (u_daughters(j,9)), &
				(u_daughters(j,10)), (u_daughters(j,11)), (u_daughters(j,12)), (u_daughters(j,13)), &
				(u_daughters(j,14)), (u_daughters(j,15))
			end if
		end do
	end if

END SUBROUTINE time_series
!**********************************

!**********************************
SUBROUTINE pop_stats(u,u_daughters,nc,ncmc,daughter_cnt,t_cumm,t_sample,E1,E2,mr1,mr2,stat_cnt,stat_marker,stat_marker2)
use globals_Parallel
implicit none
real, dimension(ncmc,15) :: u
real, dimension(ncmc,15) :: u_daughters 
real :: t_cumm,t_sample
real :: avg_protein1,tot_protein1
real :: var_cnt_protein1,var_protein1,noise_protein1
real :: avg_protein2,tot_protein2
real :: var_cnt_protein2,var_protein2,noise_protein2
real :: avg_protein3,tot_protein3
real :: var_cnt_protein3,var_protein3,noise_protein3
real :: avg_fitness, tot_fitness, avg_fitness2
real :: var_cnt_fitness, var_fitness  
real :: mutant_cnt, mutant_fraction
integer :: E1,E2,nc,ncmc,i,j,k,daughter_cnt,protein_num,stat_cnt,stat_marker,stat_marker2
real, dimension(n_runs,nint(t_end/t_restore)+1) :: mr1,mr2

	!calculate mean population protein level  
	tot_protein1=0; avg_protein1=0
	tot_protein2=0; avg_protein2=0
	tot_protein3=0; avg_protein3=0
	tot_fitness=0; avg_fitness=0
	do i=1,nc
		tot_protein1=tot_protein1+u(i,3)
		tot_protein2=tot_protein2+u(i,4)
		tot_protein3=tot_protein3+u(i,5)
		tot_fitness=tot_fitness+u(i,13)
	end do	
	avg_protein1=tot_protein1/real(nc)
	avg_protein2=tot_protein2/real(nc)
	avg_protein3=tot_protein3/real(nc)
	avg_fitness=tot_fitness/real(nc)

	!calculate variance population protein level
	var_cnt_protein1=0; var_protein1=0
	var_cnt_protein2=0; var_protein2=0
	var_cnt_protein3=0; var_protein3=0
	var_cnt_fitness=0; var_fitness=0
	do i=1,nc
		var_cnt_protein1=var_cnt_protein1+((u(i,3)-avg_protein1)**2)
		var_cnt_protein2=var_cnt_protein2+((u(i,4)-avg_protein2)**2)
		var_cnt_protein3=var_cnt_protein3+((u(i,5)-avg_protein3)**2)
		var_cnt_fitness=var_cnt_fitness+((u(i,13)-avg_fitness)**2)
	end do	
	var_protein1=var_cnt_protein1/real(nc)
	var_protein2=var_cnt_protein2/real(nc)
	var_protein3=var_cnt_protein3/real(nc)
	var_fitness=var_cnt_fitness/real(nc)
 
	mutant_cnt=0; mutant_fraction=0
	do i=1,nc
		mutant_cnt=mutant_cnt+u(i,2)
	end do
	mutant_fraction=mutant_cnt/real(nc)	

	!calculate noise population protein level 
	noise_protein1=0; noise_protein1=sqrt(var_protein1)/avg_protein1
	noise_protein2=0; noise_protein2=sqrt(var_protein2)/avg_protein2
	noise_protein3=0; noise_protein3=sqrt(var_protein3)/avg_protein3

	!print*, avg_fitness, real(daughter_cnt)/real(nc), mutant_fraction	

	!if (stat_marker==0) then
		!write statistics to file
		write(2,*) t_cumm, avg_protein1, avg_protein2, avg_protein3, avg_fitness, mutant_fraction   		
		mr1(stat_cnt,stat_marker2)=avg_fitness
		mr2(stat_cnt,stat_marker2)=mutant_fraction
	!else !get avg. [P] just before drug application
	!	master_results(1,stat_cnt)=avg_protein3
	!end if

END SUBROUTINE pop_stats 
!**********************************

END MODULE results_Parallel



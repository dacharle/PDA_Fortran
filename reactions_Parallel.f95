!This module contains the reactions for the user defined genetic network.

MODULE reactions_Parallel
implicit none
contains 

!**********************************
SUBROUTINE rxns(mu,M,P1,P2,P3)
use globals_Parallel  
implicit none
integer :: mu
real :: M,P1,P2,P3

	if (mu==1) then
		P1=P1+1 
		end if
		if (mu==2) then 
		P1=P1-1 
		end if
		if (mu==3) then
		P2=P2+1 
		end if
		if (mu==4) then
		P2=P2-1
		end if
		if (mu==5) then
		P3=P3+1
		end if 
		if (mu==6) then
		P3=P3-1
	end if 

END SUBROUTINE rxns
!**********************************

END MODULE reactions_Parallel

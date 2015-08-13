!--------------------------------------------------------
!---------------Wave function correlation----------------
!--------------------------------------------------------

subroutine f_f(q,ans)

implicit none

real :: q				
real :: ans				
		
! select case(N_bands)
! 	case(1)
! 
! 	case(2)
! 
! 	case(3)
! 
! 	case(4)
! 
! 	case(5)
! 
! 	case(6)
! 
! 	case(7)
! 
! 	case(8)
! 
! 	case(9)
! 
! end select

ans=(-5.0387e-28)*q**3+(1.3601e-18)*q**2+(-1.5339e-09)*q+(0.9954)

return
end subroutine f_f


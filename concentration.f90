!--------------------------------------------------------
!----------Concentration computation subroutine----------
!--------------------------------------------------------

subroutine concentration(k,l_k,ne,conc)

implicit none

real :: pi
integer :: j				! counter

integer :: l_k				! length of k array
real :: stk				! length of k array
real, dimension(l_k) :: k		! k-space (1st particle)
real, dimension(l_k) :: Ee		! k-space (1st particle)
real, dimension(l_k) :: ne		! k-space (1st particle)

real :: conc				! concentration


pi=3.14
stk=k(3)-k(2)
conc=sum(ne*k*stk/(pi))

open(7,file='conc.dat')
write (7,'(e15.7)'), conc

return
end subroutine concentration


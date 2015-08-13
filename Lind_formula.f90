!-----------------------------------------------------
!-----dielectrical constant for plasma screening------
!------------------Lindhard formula-------------------
!-----------------------------------------------------

subroutine Lind_formula(q,k,l_k,Eg,me,mh,Tempr,Ef_h,Ef_e,eps)

implicit none

real :: q
real :: Eg
real :: me
real :: mh
real :: Tempr
real :: eps
real :: eps1
real :: eps2
real :: Ef_e
real :: Ef_h
real :: Ee
real :: Ee1
real :: Eh
real :: Eh1
real :: ne
real :: ne1
real :: nh
real :: nh1
integer :: l_k
real, dimension(l_k) :: k
real :: k1
real, dimension(10002) :: phi
real :: h				! Planck constant
real :: e				! elementary charge
real :: m0				! electron mass
real :: pi=3.14				! pi	
real :: eps0				! dielectric constant
real :: kb				! Boltzmann constant
real :: c				! light speed in vacuum
complex :: ij				! imaginary unit
integer :: j
integer :: j1

h=1.05e-34				!Planck constant
 c=3e8					!Light speed
e=1.6e-19				!Electron charge
m0=9.1e-31				!Electon mass
eps0=8.85e-12				!dielectric constant
ij=sqrt(CMPLX(-1.0))			!imaginary init
kb=1.38e-23				!Boltzmann constant

phi=(/(j1*2*pi/10001,j1=0,(10001))/)

eps1=0.0
eps2=0.0

do j1=1,size(k)
    do j=1,10001

	k1=sqrt(k(j1)**2+q**2+2.0*k(j1)*q*cos(phi(j)))

	Ee=Eg+h*(k(j1)**2.0)/(2.0*me)*h
	Ee1=Eg+h*(k1**2.0)/(2.0*me)*h
	Eh=h*(k(j1)**2.0)/(2.0*mh)*h
	Eh1=h*(k1**2.0)/(2.0*mh)*h

	ne=1/(1+exp((Ee-Ef_e)/kb/Tempr))
	ne1=1/(1+exp((Ee1-Ef_e)/kb/Tempr))
	nh=1/(1+exp((Eh-Ef_h)/kb/Tempr))
	nh1=1/(1+exp((Eh1-Ef_h)/kb/Tempr))

	eps1=eps1+k(j1)*(ne1-ne)/(Ee1-Ee)*(phi(3)-phi(2))*(k(3)-k(2))
	eps2=eps2+k(j1)*(nh1-nh)/(Eh1-Eh)*(phi(3)-phi(2))*(k(3)-k(2))
    enddo
enddo

eps=2/q*(eps1+eps2)

end subroutine Lind_formula
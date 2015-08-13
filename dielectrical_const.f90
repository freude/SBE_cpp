!-----------------------------------------------------
!-----dielectrical constant for plasma screening------
!--------------Plasmon-Pole approximation-------------
!-----------------------------------------------------

subroutine dielectrical_const(q,mh,me,nh,ne,conc,Tempr,eps)

implicit none

integer :: l_k


real :: h			!Planck constant
real :: el			!Electron charge
real :: eps0			!dielectric constant
real :: eps_b			!dielectric constant
real :: kb			!Boltzmann constant
real :: pi

real :: q
real :: mh
real :: me
real :: nh
real :: ne
real :: omega
real :: Tempr

real :: cf
real :: mr
double precision :: kappa
double precision :: omega_pl
real :: omega_q
real :: conc
real :: eps

!----------------Definition of constants-----------------

h=1.05e-34				!Planck constant
el=1.6e-19				!Electron charge
eps0=8.85e-12				!dielectric constant
kb=1.38e-23				!Boltzmann constant
pi=3.14
eps_b=12.3

!---------------------------------------------------------

 cf=0.017

!---------------------------------------------------------

mr=me/(me+mh)*mh

kappa=2*mh/eps0/h*el/h*el*ne+2*me/eps0/h*el/h*el*(1.0-nh)
omega_pl=2*pi*conc*el*el*q/eps0/mr

omega_q=omega_pl*(1+q/kappa)+cf*q**4
eps=1.0-omega_pl/omega_q

!open(17,file='eps.dat')
!write (17,'(e15.3)'), eps
!close(17)

!open(17,file='eps1.dat')
!write (17,'(e10.3)'), omega_q
!close(17)

! subroutine appr(q,l_q,eps,4,coef)

end subroutine dielectrical_const
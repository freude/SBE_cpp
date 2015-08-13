subroutine screening(l_k,k1,conc,gway,ccc)

use fitting

implicit none

integer, parameter :: dp = selected_real_kind(15, 307)
integer :: l_k						!length of k array

real, parameter :: h=1.054e-34		! Planck constant
real, parameter :: e=1.602e-19		! elementary charge
real, parameter :: m0=9.109e-31		! electron mass
real, parameter :: kb=1.38e-23		! electron mass
real, parameter :: pi=3.14		! electron mass
real, parameter :: eps0=8.85e-12				!dielectric constant

complex :: ij						! imaginary unit

integer :: j
integer :: j1
integer :: j2
integer :: j3

!---------------------------k-space-----------------------------

real :: k_min			
real :: k_max			
real :: stk				
real, dimension(l_k) :: k1			! k-space (1st particle)
real, dimension(l_k) :: k			! k-space (1st particle)
real, dimension(l_k) :: q			! wave vector difference for interacted particles
real, dimension(10002) :: phi
real :: dphi


real :: kq

real :: Ee1			! Conduction band structure
real, dimension(l_k) :: Ee			! Conduction band structure
real, dimension(l_k) :: Ee0			! Conduction band structure
real :: Eh1			! Valence band structure
real, dimension(l_k) :: Eh			! Valence band structure
real, dimension(l_k) :: Eh0			! Valence band structure

real :: ne1			! Conduction band structure
real :: ne			! Conduction band structure

real :: nh1			! Conduction band structure
real :: nh			! Conduction band structure

!----------------Material parameters---------------------

real :: me							! Electrons' effective mass
real :: mh							! Holes' effective mass
real :: Tempr						! temperature
real :: conc						! concentration of carriers
real :: Eg							! Band gap
real :: delta							! Band gap

complex, dimension(l_k) :: I
complex, dimension(l_k) :: I1
complex, dimension(l_k) :: eps
real :: Ef_e
real :: Ef_h
real :: epsb
real(dp), dimension(l_k) :: aaa
real(dp), dimension(l_k) :: bbb
real(dp), dimension(9) :: ccc
real(dp), dimension(3) :: ccc0
character (len=100) :: gway

!------------------------------------------------------------

ij=sqrt(CMPLX(-1.0))

Eg=1.52*e								!nominal band gap
delta=0.0001*e							!nominal band gap
me=0.02*m0								!electrons effective mass
mh=0.5*m0								!holes effective mass

Tempr=300.0
!conc=8.25e16
epsb=11.5

!----------------Formation of arrays---------------------

k_min=0.0					! min kx/ky
k_max=4.0e9					! max kx/ky
stk=(k_max-k_min)/l_k				! step in k-grid

k=(/(j*stk+k_min,j=0,int(l_k-1))/)		! k array
q=k

phi=(/(j1*2*pi/10001,j1=0,(10001))/)

stk=k(4)-k(3)
dphi=phi(5)-phi(4)

!------------------------------------------------------

open(7,file=trim(gway)//'E1.dat')
read (7,'(e15.5)'), Eh0
close(7)

Eh0=Eh0*e/h

open(7,file=trim(gway)//'Ee1.dat')
read (7,'(e15.5)'), Ee0
 close(7)

Ee0=Ee0*e/h

aaa=k1
bbb=Eh0
ccc0=polyfit(aaa,bbb,2)
Eh=ccc0(1)+ccc0(2)*k+ccc0(3)*(k**2.0)
bbb=Ee0
ccc0=polyfit(aaa,bbb,2)
Ee=ccc0(1)+ccc0(2)*k+ccc0(3)*(k**2.0)

!------------------------------------------------------

call Fermi_level(k,l_k,1,'c',Tempr,conc,Ef_e,h*Ee)
call Fermi_level(k,l_k,1,'v',Tempr,conc,Ef_h,h*Eh)

!------------------------------------------------------
  
do j1=1,size(q)	
	do j2=1,size(k)	
		do j3=1,size(phi)
			
			kq=sqrt(k(j2)**2.0+q(j1)**2.0-2.0*k(j2)*q(j1)*cos(phi(j3)))
		
			Ee1=Eg/h+h*(kq**2.0)/(2.0*me)
			Eh1=-h*(kq**2.0)/(2.0*mh)
		
			ne1=1.0/(1+exp((Ee1*h-Ef_e)/kb/Tempr))
			ne=1.0/(1+exp((Ee(j2)*h-Ef_e)/kb/Tempr))
		
			nh1=1.0/(1+exp((Eh1*h-Ef_h)/kb/Tempr))
			nh=1.0/(1+exp((Eh(j2)*h-Ef_h)/kb/Tempr))
		
			I(j1)=I(j1)+k(j2)/q(j1)*(ne1-ne)/(Ee1*h-Ee(j2)*h+ij*delta)*stk*dphi			
			I1(j1)=I1(j1)+k(j2)/q(j1)*(nh1-nh)/(Eh1*h-Eh(j2)*h+ij*delta)*stk*dphi
			
		enddo
	enddo

		!write(*,'(A)',ADVANCE='NO')'*'
		
		print *,'j1=',j1
enddo

eps=1.0-2.0*(I+I1)*(e**2.0)/(8.0*pi*pi*eps0*epsb)
eps(1)=eps(2)

aaa=q
bbb=1.0/real(eps)


ccc=polyfit(aaa,bbb,8)
print *,'ccc=',ccc

open(7,file='I.dat')
write (7,'(e15.7)'), bbb

open(7,file='eps_8d25e16.dat')
write (7,'(e15.7)'), real(eps)

open(7,file='eps_8d25e16i.dat')
write (7,'(e15.7)'), imag(eps)

end subroutine screening

!--------------------------------------------------------
!-----------Fermi level computation subroutine-----------
!--------------------------------------------------------

subroutine Fermi_level(k,l_k,N_bands,flag,Tempr,conc,Ef,Ee1,Ee2,Ee3,Ee4,Ee5,Ee6,Ee7,Ee8,Ee9)

implicit none

real :: pi
real :: h
real :: m0
real :: me
real :: mh
integer :: j				! counter

integer :: l_k				! length of k array
real :: stk				! length of k array
real, dimension(l_k) :: k		! k-space (1st particle)
real, dimension(l_k) :: Ee		! k-space (1st particle)
real, dimension(l_k) :: ne		! k-space (1st particle)

real Tempr
real :: kb
real :: el
real :: conc				! concentration
real :: conc1				! concentration
real,dimension(2) :: diff				! concentration

real :: Ef				! concentration

real :: Ef_min				! concentration
real :: Ef_max				! concentration
real :: st_Ef				! concentration
real,dimension(1000) :: Ef_pr		! concentration
real,dimension(1000) :: a		! concentration
integer :: N_bands
character :: flag

real, dimension(l_k) :: Ee1		! Conduction band structure
real, dimension(l_k) :: Ee2		! Conduction band structure
real, dimension(l_k) :: Ee3		! Conduction band structure
real, dimension(l_k) :: Ee4		! Conduction band structure
real, dimension(l_k) :: Ee5		! Conduction band structure
real, dimension(l_k) :: Ee6		! Conduction band structure
real, dimension(l_k) :: Ee7		! Conduction band structure
real, dimension(l_k) :: Ee8		! Conduction band structure
real, dimension(l_k) :: Ee9		! Conduction band structure

integer, parameter :: l_E=300
real :: stE
real :: E_max
real :: E_min
real, dimension(l_E) :: dos_el
real, dimension(l_E) :: dos_h
real, dimension(l_E) :: E
real :: damp
real :: lz

kb=1.38e-23
el=1.6e-19
h=1.05e-34				!Planck constant
m0=9.1e-31				!Electon mass
pi=3.14

me=0.0665*m0				!electrons effective mass
mh=0.234*m0				!holes effective mass

E_min=-0.5					! min kx/ky
E_max=0.5					! max kx/ky
stE=(E_max-E_min)/l_E				! step in k-grid
damp=0.005

E=(/(j*stE+E_min,j=0,int(l_E-1))/)		! k array

!--------------------------------------------------------

stk=k(3)-k(2)
lz=1.0

!--------------------------------------------------------

diff(1)=conc**2

if (flag .EQ. 'c') then
	Ef_min=-0.5						! min kx/ky
	Ef_max=0.5						! max kx/ky
	st_Ef=(Ef_max-Ef_min)/1000				! step in k-grid
	Ef_pr=(/(j*st_Ef+Ef_min+minval(Ee1)/el,j=0,(1000-1))/)	! k array
	do j=1,1000		

		select case(N_bands)
		case(1)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))
		case(2)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))
		case(3)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))
		case(4)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))

		case(5)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee5-Ef_pr(j)*el)/kb/Tempr))
		case(6)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee6-Ef_pr(j)*el)/kb/Tempr))
		case(7)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee7-Ef_pr(j)*el)/kb/Tempr))
		case(8)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee7-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee8-Ef_pr(j)*el)/kb/Tempr))

		case(9)
			ne=1/(1+exp((Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee7-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee8-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((Ee9-Ef_pr(j)*el)/kb/Tempr))
		end select

		conc1=2/lz*sum(ne*k*stk/2/pi)
		diff(2)=abs(conc-conc1)
		a(j)=diff(2)

		if (diff(1)>diff(2)) then
			diff(1)=diff(2)
			Ef=Ef_pr(j)
		endif
	enddo
else
	Ef_min=-0.5					! min kx/ky
	Ef_max=0.5					! max kx/ky
	st_Ef=(Ef_max-Ef_min)/1000			! step in k-grid
	Ef_pr=(/(j*st_Ef+Ef_min+minval(-Ee1)/el,j=0,(1000-1))/)		! k array
	do j=1,1000
		select case(N_bands)
		case(1)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))
		case(2)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))
		case(3)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))
		case(4)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))
		case(5)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee5-Ef_pr(j)*el)/kb/Tempr))
		case(6)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee6-Ef_pr(j)*el)/kb/Tempr))
		case(7)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee7-Ef_pr(j)*el)/kb/Tempr))
		case(8)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee7-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee8-Ef_pr(j)*el)/kb/Tempr))
		case(9)
			ne=1/(1+exp((-Ee1-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee2-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee3-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee4-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee5-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee6-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee7-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee8-Ef_pr(j)*el)/kb/Tempr))+&
			&1/(1+exp((-Ee9-Ef_pr(j)*el)/kb/Tempr))
		end select

		conc1=2/lz*sum(ne*k*stk/2/pi)
		diff(2)=abs(conc-conc1)
		a(j)=diff(2)

		if (diff(1)>diff(2)) then
			diff(1)=diff(2)
			Ef=-Ef_pr(j)
		endif
	enddo
endif

open(7,file='ef.dat')
write (7,'(e15.7)'), Ef
close(7)

open(8,file='diff.dat')
write (8,'(e15.7)'), diff
close(8)

open(9,file='a.dat')
write (9,'(e15.7)'), a
close(9)

open(9,file='dos_el.dat')
write (9,'(e15.7)'), dos_el
close(9)

open(9,file='dos_h.dat')
write (9,'(e15.7)'), dos_h
close(9)

print *,'Fermi level=',Ef

Ef=Ef*el

return
end subroutine Fermi_level 
subroutine polarization(k,l_k,Eh,Ee,mu,Ef_h,Ef_e,Tempr,l_f,fff,PSr)

implicit none

!----------------Matrix sizes----------------------------

integer, parameter :: l_t=15000		!length of time array
integer :: l_k				!length of k array
integer :: l_f

!----------------Constantes------------------------------

real :: h				! Planck constant
real :: e				! elementary charge
real :: m0				! electron mass
real :: pi=3.14				! pi	
real :: eps0				! dielectric constant
real :: kb				! Boltzmann constant
real :: c				! light speed in vacuum
complex :: ij				! imaginary unit

!----------------Material parameters---------------------

real :: eps				! Permittivity
real :: Eg				! Band gap
real :: Ef_e				! Fermi level
real :: Ef_h				! Fermi level
real :: me				! Electrons' effective mass
real :: mh				! Holes' effective mass
real :: Tempr				! temperature
real :: conc				! concentration of carriers
real :: Vol				! active volume
real :: n_reff				! refractive index

!---------------------------Counters----------------------------

integer :: j1
integer :: j2
integer :: j  

!----------------------------time-------------------------------

real :: t_min
real :: t_max
real :: stt

real, dimension(l_t) :: t

!---------------------------k-space-----------------------------

real :: k_min				!
real :: k_max				!
real :: stk				!

real, dimension(l_k) :: k		! k-space (1st particle)
real, dimension(l_k) :: k1		! k-space (2nd particle)


!---------------------Fourier frequencies------------------------

real, dimension(l_f) :: fff
!----------------Single-particles characteristics---------------

real, dimension(l_k) :: mu		! dipole matrix element
real, dimension(l_k) :: ne		! Electrons distribution function
real, dimension(l_k) :: nh		! Holes distribution function
real, dimension(l_k) :: omega		! Single-paticle transition frequency
real, dimension(l_k) :: Ee		! Conduction band structure
real, dimension(l_k) :: Eh		! Valence band structure

!---------------------Additional variables----------------------

complex :: RS

complex :: kk1				! Runge-Kutta RS
complex :: kk2				! Runge-Kutta RS
complex :: kk3				! Runge-Kutta RS
complex :: kk4				! Runge-Kutta RS

!-----------------Many-particle characteristics-----------------

real, dimension(l_k) :: exce		! Exchange energy
real :: damp				! dephasing
real :: q				! diff. between momenta of interacting carriers
real, dimension(l_k,l_k) :: V		! interaction matrix

!-------------------------Polarization--------------------------

complex, dimension(l_t,l_k) :: pp
complex, dimension(l_t) :: P
complex, dimension(l_f) :: PS
real, dimension(l_f) :: PSr

!------------------------Electric field-------------------------

complex, dimension(l_t) :: E_ft
complex, dimension(l_f) :: ES
complex :: E_field

!----------------Definition of constants-----------------

h=1.05e-34				!Planck constant
 c=3e8					!Light speed
e=1.6e-19				!Electron charge
m0=9.1e-31				!Electon mass
eps0=8.85e-12				!dielectric constant
ij=sqrt(CMPLX(-1.0))			!imaginary init
kb=1.38e-23				!Boltzmann constant

!------------------Material parameters-------------------

!eps=12.3				!permitivity
!n_reff=3.61				!refraction index
!Eg=2.5679*e				!nominal band gap
!me=0.148*m0				!electrons effective mass
!mh=0.234*m0				!holes effective mass

eps=12.3				!permitivity
n_reff=3.61				!refraction index
!Eg=1.519*e				!nominal band gap
me=0.0665*m0				!electrons effective mass
mh=0.234*m0				!holes effective mass

!----------------Formation of arrays---------------------


t_min=0.0					! min time
t_max=0.15e-11					! max time					
stt=(t_max-t_min)/l_t				! step in t-grid

t=(/ ((j*stt+t_min),j=0,(int(l_t)-1)) /)	! time array
stk=k(4)-k(3)
k1=k						! k' array

!--------------------------------------------------------

damp=0.01*e					!damping

!----------------------Plasma parameters-----------------

!Ef_e=Eg+0.02*e
!Ef_h=-0.003*e
!Tempr=300
!conc=1.0e14
!conc=0.0

!---------------------Transition frequency----------------------

omega=Ee-Eh
Eg=h*omega(1)

!-----------------Distribution functions--------------------

ne=1.0/(1+exp((Ee*h-Ef_e)/kb/Tempr))
nh=1.0-1.0/(1+exp(-(Eh*h-Ef_h)/kb/Tempr))

!ne=0.0
!nh=1.0

call concentration(k,l_k,1.0-nh,conc)
print *,'concentration=',conc

!-----------------------------------------------------------

!call dielectrical_const(k,l_k,mh,me,0.0*minval(Ee),nh(1),ne(1),Tempr,conc,eps)

!----------------------Interaction matrix--------------------

 call int_matrix(k,k1,l_k,eps,V,mh,me,Tempr,conc,ne(1),nh(1))
 call exchange(k,l_k,ne,1.0-nh,V,exce)

!-----------Solving semiconductor Bloch equations------------

print *,'Eg=',Eg/e
mu=mu*1.0e-27
vol=1.0
P=0.0

do j2=2,size(t)	
     do j1=1,size(k)
        RS=-ij*(omega(j1)-Eg/h-exce(j1)/h)*pp(j2-1,j1)-ij*(ne(j1)-nh(j1))*(mu(j1)*E_field(t(j2-1))+&
	&sum(V(j1,:)*pp(j2-1,:)*k*stk))/h-damp*pp(j2-1,j1)/h-ne(j1)*(1.0-nh(j1))*E_field(t(j2-1))/h

	kk1=RS

	kk2=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk1/2)-ij*(ne(j1)-nh(j1))*&
	&(mu(j1)*E_field(t(j2-1)+stt/2)+&
	&sum(V(j1,:)*pp(j2-1,:)*k*stk))/h-damp*(pp(j2-1,j1)+stt*kk1/2)/h-&
	&ne(j1)*(1.0-nh(j1))*E_field(t(j2-1)+stt/2)/h

	kk3=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk2/2)-ij*(ne(j1)-nh(j1))*&
	&(mu(j1)*E_field(t(j2-1)+stt/2)+&
	&sum(V(j1,:)*pp(j2-1,:)*k*stk))/h-damp*(pp(j2-1,j1)+stt*kk2/2)/h-&
	&ne(j1)*(1.0-nh(j1))*E_field(t(j2-1)+stt/2)/h

	kk4=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk3)-ij*(ne(j1)-nh(j1))*&
	&(mu(j1)*E_field(t(j2-1)+stt)+&
	&sum(V(j1,:)*pp(j2-1,:)*k*stk))/h-damp*(pp(j2-1,j1)+stt*kk3)/h-&
	&ne(j1)*(1.0-nh(j1))*E_field(t(j2-1)+stt)/h

        pp(j2,j1)=pp(j2-1,j1)+(stt/6)*(kk1+2*kk2+2*kk3+kk4)
!        pp(j2,j1)=pp(j2-1,j1)+(stt)*kk1
  	P(j2)=P(j2)+2.0*pi/vol*mu(j1)*k(j1)*pp(j2,j1)*stk
    enddo
    E_ft(j2)=E_field(t(j2-1))



enddo

!-----------------Fourie transformation------------------

do j=1,l_f  
	ES(j)=sum(E_ft*exp(ij*fff(j)*t)*stt)
	PS(j)=sum(P*exp(ij*fff(j)*t)*stt)/(4.0*pi*eps0*eps)
enddo

!PS=(fff)*(PS/(ES*Vol))/(c*n_reff)
!PS=(PS/(4.0*eps0*eps*ES))

PSr=(fff+Eg/h)*abs(PS)/(c*n_reff)

!------------------------Data output---------------------

!ne=ne-n

open(7,file='p.dat')
write (7,'(e10.3)'), p
close(7)
open(9,file='ps11r.dat')
write (9,'(e10.3)'), PSr
close(9)
open(10,file='ne.dat')
write (10,'(e10.5)'), E_ft
close(10)
open(11,file='nh.dat')
write (11,'(e12.5)'), nh
close(11)
open(11,file='ne.dat')
write (11,'(e12.5)'), ne
close(11)
open(12,file='Ee.dat')
write (12,'(e12.3)'), Ee
close(12)
open(15,file='Eh.dat')
write (15,'(e12.3)'), Eh
close(15)

end subroutine polarization

subroutine polarization1(n_charact,k,l_k,Eh,Ee,mu,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,PSr)

implicit none

!----------------Matrix sizes----------------------------

integer, parameter :: l_t=5000		!length of time array
integer :: l_k						!length of k array
integer :: l_f

!----------------Constantes------------------------------

real :: h					! Planck constant
real :: e					! elementary charge
real :: m0					! electron mass
real :: pi=3.14				! pi	
real :: eps0				! dielectric constant
real :: kb					! Boltzmann constant
real :: c					! light speed in vacuum
complex :: ij				! imaginary unit

!----------------Material parameters---------------------

real :: eps					! Permittivity
real :: Eg					! Band gap
real :: Ef_e				! Fermi level
real :: Ef_h				! Fermi level
real :: me					! Electrons' effective mass
real :: mh					! Holes' effective mass
real :: Tempr				! temperature
real :: conc				! concentration of carriers
real :: conc1				! concentration of carriers
real :: Vol					! active volume
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
real, dimension(l_k*l_k) :: Co		! interaction matrix
real, dimension(l_k*l_k) :: Co1		! interaction matrix
real, dimension(l_k,l_k) :: Corr		! interaction matrix
real, dimension(l_k,l_k) :: Corr1		! interaction matrix
integer :: n_charact
integer :: Ncut

!----------------Definition of constants-----------------

h=1.05e-34				!Planck constant
 c=3.0e8					!Light speed
e=1.6e-19				!Electron charge
m0=9.1e-31				!Electon mass
eps0=8.85e-12				!dielectric constant
ij=sqrt(CMPLX(-1.0))			!imaginary init
kb=1.38e-23				!Boltzmann constant
pi=3.14;

!------------------Material parameters-------------------

!eps=12.3				!permitivity
!n_reff=3.61				!refraction index
!Eg=2.5679*e				!nominal band gap
!me=0.148*m0				!electrons effective mass
!mh=0.234*m0				!holes effective mass

eps=12.3				!permitivity
n_reff=3.61				!refraction index
!Eg=1.519*e				!nominal band gap
me=0.067*m0				!electrons effective mass
mh=0.377*m0				!holes effective mass

!!----------------Formation of arrays---------------------

t_min=0.0					! min time
t_max=0.25e-11					! max time
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
nh=1.0/(1+exp((Eh*h-Ef_h)/kb/Tempr))

ne=0.0
nh=1.0

call concentration(k,l_k,1.0-nh,conc1)
print *,'concentration1=',conc1

!-----------------------------------------------------------

! call dielectrical_const(k,l_k,mh,me,0.0*minval(Ee),nh(1),ne(1),Tempr,conc,eps)

!----------------------Interaction matrix--------------------

!	call int_matrix(k,k1,l_k,eps,V,mh,me,Tempr,conc,ne(1),nh(1))
call exchange(k,l_k,ne,1.0-nh,V,exce)

!-----------Solving semiconductor Bloch equations------------

print *,'Eg=',Eg/e
mu=sqrt(mu)
vol=1.0
P=0.0

! open(7,file='/home/mk/results/InGaN_2nm_0d1_2003/dependencies_s/0d7_0d7_0d4_7d6/Corr.dat')
! read (7,'(e15.5)'), Co
!  close(7)
! 
! 	Co=Co/0.25e20
! 
! open(7,file='/home/mk/results/InGaN_2nm_0d1_2003/dependencies_s/0d7_0d7_0d4_7d6/Corr1.dat')
! read (7,'(e15.5)'), Co1
!  close(7)
! 
! 	Co1=Co1/0.25e20
! 
! do j1=1,size(k)
! 	Corr(:,j1)=Co((1+(j1-1)*size(k)):(j1*size(k)))
! 	Corr1(:,j1)=Co1((1+(j1-1)*size(k)):(j1*size(k)))
! enddo
! 
! open(7,file='Corrr.dat')
! write (7,'(400e15.5)'), Corr
! close(7)

Ncut=200;

do j2=2,size(t)	
!print *,'t=',j2
     do j1=1,size(k)
     
		if ((j1-Ncut)<1) then
		
			RS=-ij*(omega(j1)-Eg/h-exce(j1)/h)*pp(j2-1,j1)-ij*(ne(j1)-nh(j1))*(mu(j1)*E_field(t(j2-1))+&
			&sum(V(j1,1:(j1+Ncut))*pp(j2-1,1:(j1+Ncut))*k(1:(j1+Ncut))*stk))/h-damp*pp(j2-1,j1)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk1=RS

			kk2=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk1/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,1:(j1+Ncut))*pp(j2-1,1:(j1+Ncut))*k(1:(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk1/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk3=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk2/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,1:(j1+Ncut))*pp(j2-1,1:(j1+Ncut))*k(1:(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk2/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk4=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk3)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt)+&
			&sum(V(j1,1:(j1+Ncut))*pp(j2-1,1:(j1+Ncut))*k(1:(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk3)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			pp(j2,j1)=pp(j2-1,j1)+(stt/6)*(kk1+2*kk2+2*kk3+kk4)
	!        pp(j2,j1)=pp(j2-1,j1)+(stt)*kk1
			P(j2)=P(j2)+mu(j1)*k(j1)*pp(j2,j1)*stk		
		endif

		if ((j1+Ncut)>size(k)) then
		
			RS=-ij*(omega(j1)-Eg/h-exce(j1)/h)*pp(j2-1,j1)-ij*(ne(j1)-nh(j1))*(mu(j1)*E_field(t(j2-1))+&
			&sum(V(j1,(j1-Ncut):size(k))*pp(j2-1,(j1-Ncut):size(k))*k((j1-Ncut):size(k))*stk))/h-damp*pp(j2-1,j1)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk1=RS

			kk2=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk1/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,(j1-Ncut):size(k))*pp(j2-1,(j1-Ncut):size(k))*k((j1-Ncut):size(k))*stk))/h-damp*(pp(j2-1,j1)+stt*kk1/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk3=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk2/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,(j1-Ncut):size(k))*pp(j2-1,(j1-Ncut):size(k))*k((j1-Ncut):size(k))*stk))/h-damp*(pp(j2-1,j1)+stt*kk2/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk4=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk3)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt)+&
			&sum(V(j1,(j1-Ncut):size(k))*pp(j2-1,(j1-Ncut):size(k))*k((j1-Ncut):size(k))*stk))/h-damp*(pp(j2-1,j1)+stt*kk3)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			pp(j2,j1)=pp(j2-1,j1)+(stt/6)*(kk1+2*kk2+2*kk3+kk4)
	!        pp(j2,j1)=pp(j2-1,j1)+(stt)*kk1
			P(j2)=P(j2)+mu(j1)*k(j1)*pp(j2,j1)*stk		
		endif
		
		if (((j1-Ncut)>=1).AND.((j1+Ncut)<=size(k))) then
		
			RS=-ij*(omega(j1)-Eg/h-exce(j1)/h)*pp(j2-1,j1)-ij*(ne(j1)-nh(j1))*(mu(j1)*E_field(t(j2-1))+&
			&sum(V(j1,(j1-Ncut):(j1+Ncut))*pp(j2-1,(j1-Ncut):(j1+Ncut))*k((j1-Ncut):(j1+Ncut))*stk))/h-damp*pp(j2-1,j1)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk1=RS

			kk2=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk1/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,(j1-Ncut):(j1+Ncut))*pp(j2-1,(j1-Ncut):(j1+Ncut))*k((j1-Ncut):(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk1/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk3=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk2/2)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt/2)+&
			&sum(V(j1,(j1-Ncut):(j1+Ncut))*pp(j2-1,(j1-Ncut):(j1+Ncut))*k((j1-Ncut):(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk2/2)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			kk4=-ij*(omega(j1)-Eg/h-exce(j1)/h)*(pp(j2-1,j1)+stt*kk3)-ij*(ne(j1)-nh(j1))*&
			&(mu(j1)*E_field(t(j2-1)+stt)+&
			&sum(V(j1,(j1-Ncut):(j1+Ncut))*pp(j2-1,(j1-Ncut):(j1+Ncut))*k((j1-Ncut):(j1+Ncut))*stk))/h-damp*(pp(j2-1,j1)+stt*kk3)/h!+sum((Co((1+(j1-1)*size(k)):(j1*size(k)))
	!+ij*Co1((1+(j1-1)*size(k)):(j1*size(k))))*pp(j2-1,:)*k*stk)

			pp(j2,j1)=pp(j2-1,j1)+(stt/6)*(kk1+2*kk2+2*kk3+kk4)
	!        pp(j2,j1)=pp(j2-1,j1)+(stt)*kk1
			P(j2)=P(j2)+mu(j1)*k(j1)*pp(j2,j1)*stk
		
		endif

        

    enddo

	E_ft(j2)=E_field(t(j2-1))
	print *,'P(',t(j2),')=',P(j2)
!	print *,'E(',t(j2),')=',E_ft(j2)

! call METAFL('CONS')
! CALL DISINI()
!       CALL PAGERA()
!       CALL COMPLX()
 
!       CALL COLOR('RED')
!       CALL ELLIPS(t(j2)/0.15e-11, P(j2)/0.12366E-27, 2, 2)
 
!       CALL DISFIN()

	print *,n_charact,' P(',t(j2),')=',P(j2)

enddo

!-----------------Fourie transformation------------------

do j=1,l_f  
	ES(j)=sum(E_ft*exp(ij*fff(j)*t)*stt)
	PS(j)=sum(P*exp(ij*fff(j)*t)*stt)/(4.0*pi*eps0*eps)
enddo

!PS=(fff)*(PS/(ES*Vol))/(c*n_reff)
! PS=(PS/(4.0*eps0*eps*ES))
PSr=(fff+Eg/h)*imag(PS/ES)/(c*n_reff)

!------------------------Data output---------------------

!ne=ne-n

open(7,file='p.dat')
write (7,'(e10.3)'), real(p)
close(7)
open(7,file='pp.dat')
write (7,'(300e10.3)'), transpose(imag(pp))
close(7)
open(9,file='ps11r.dat')
write (9,'(e12.5)'), PSr
close(9)
open(10,file='E_ft.dat')
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

end subroutine polarization1

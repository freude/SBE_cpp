program gain

!use dislin

implicit none

!----------------Matrix sizes----------------------------

integer, parameter :: l_k=400		!length of k array
integer, parameter :: l_f=9000		!length of frequency array

!----------------Definition of constants-----------------

real, parameter :: h=1.054e-34		! Planck constant
real, parameter :: e=1.602e-19		! elementary charge
real, parameter :: m0=9.109e-31		! electron mass
real, parameter :: kb=1.38e-23		! electron mass
real, parameter :: pi=3.14			! pi	
real, parameter ::c=3.0e8			!Light speed
real,parameter :: n_reff=3.61		!refraction index

!----------------Material parameters---------------------

real :: me							! Electrons' effective mass
real :: mh							! Holes' effective mass
real :: Tempr						! temperature
real :: conc						! concentration of carriers


!---------------------------Counters----------------------------

integer :: j
integer :: j1
integer :: j2
integer :: j3
integer :: j4
integer :: j5
integer :: j6
integer :: j7
integer :: j8
integer :: j9
integer :: j10
integer :: j11
integer :: j12
integer :: j13
integer :: j14
integer :: j15
integer :: j16
integer :: j17

integer :: j18
integer :: j19
integer :: j20
integer :: j21
integer :: j22
integer :: j23

!---------------------------k-space-----------------------------

real :: k_min			
real :: k_max			
real :: stk				
real, dimension(l_k) :: k			! k-space (1st particle)
real, dimension(l_k-1) :: q			! wave vector difference for interacted particles

!---------------------Fourier frequencies------------------------

real :: fff_min
real :: fff_max
real :: stf

real, dimension(l_f) :: fff

!----------------Single-particles characteristics---------------

real :: Eg							! Band gap
real :: Eg11						! Band gap
real :: Eg21						! Band gap
real :: Eg31						! Band gap
real :: Eg41						! Band gap
real :: Eg51						! Band gap
real :: Eg61						! Band gap
real :: Eg71						! Band gap
real :: Eg81						! Band gap
real :: Eg91						! Band gap
real :: Eg101						! Band gap
real :: Eg111						! Band gap
real :: Eg121						! Band gap


real :: Eg12						! Band gap
real :: Eg22						! Band gap
real :: Eg32						! Band gap
real :: Eg42						! Band gap
real :: Eg52						! Band gap
real :: Eg62						! Band gap
real :: Eg72						! Band gap
real :: Eg82						! Band gap
real :: Eg92						! Band gap
real :: Eg102						! Band gap
real :: Eg112						! Band gap
real :: Eg122						! Band gap


real :: Eg13						! Band gap
real :: Eg23						! Band gap
real :: Eg33						! Band gap
real :: Eg43						! Band gap
real :: Eg53						! Band gap
real :: Eg63						! Band gap
real :: Eg73						! Band gap


real, dimension(l_k) :: mu11		! dipole matrix element
real, dimension(l_k) :: mu12		! dipole matrix element
real, dimension(l_k) :: mu13		! dipole matrix element
real, dimension(l_k) :: mu14		! dipole matrix element
real, dimension(l_k) :: mu15		! dipole matrix element
real, dimension(l_k) :: mu16		! dipole matrix element
real, dimension(l_k) :: mu17		! dipole matrix element
real, dimension(l_k) :: mu18		! dipole matrix element
real, dimension(l_k) :: mu19		! dipole matrix element
real, dimension(l_k) :: mu110		! dipole matrix element
real, dimension(l_k) :: mu111		! dipole matrix element
real, dimension(l_k) :: mu112		! dipole matrix element
real, dimension(l_k) :: mu21		! dipole matrix element
real, dimension(l_k) :: mu22		! dipole matrix element
real, dimension(l_k) :: mu23		! dipole matrix element
real, dimension(l_k) :: mu24		! dipole matrix element
real, dimension(l_k) :: mu25		! dipole matrix element
real, dimension(l_k) :: mu26		! dipole matrix element
real, dimension(l_k) :: mu27		! dipole matrix element
real, dimension(l_k) :: mu28		! dipole matrix element
real, dimension(l_k) :: mu29		! dipole matrix element
real, dimension(l_k) :: mu210		! dipole matrix element
real, dimension(l_k) :: mu211		! dipole matrix element
real, dimension(l_k) :: mu212		! dipole matrix element
real, dimension(l_k) :: mu31		! dipole matrix element
real, dimension(l_k) :: mu32		! dipole matrix element
real, dimension(l_k) :: mu33		! dipole matrix element
real, dimension(l_k) :: mu34		! dipole matrix element
real, dimension(l_k) :: mu35		! dipole matrix element
real, dimension(l_k) :: mu36		! dipole matrix element

real, dimension(l_k) :: Ee1			! Conduction band structure
real, dimension(l_k) :: Ee2			! Conduction band structure
real, dimension(l_k) :: Ee3			! Conduction band structure
real, dimension(l_k) :: Ee4			! Conduction band structure
real, dimension(l_k) :: Ee5			! Conduction band structure
real, dimension(l_k) :: Ee6			! Conduction band structure
real, dimension(l_k) :: Ee7			! Conduction band structure
real, dimension(l_k) :: Ee8			! Conduction band structure
real, dimension(l_k) :: Ee9			! Conduction band structure
real, dimension(l_k) :: Ee10			! Conduction band structure
real, dimension(l_k) :: Ee11			! Conduction band structure
real, dimension(l_k) :: Ee12			! Conduction band structure

real, dimension(l_k) :: Eh1			! Valence band structure
real, dimension(l_k) :: Eh2			! Valence band structure
real, dimension(l_k) :: Eh3			! Valence band structure
real, dimension(l_k) :: Eh4			! Valence band structure
real, dimension(l_k) :: Eh5			! Valence band structure
real, dimension(l_k) :: Eh6			! Valence band structure
real, dimension(l_k) :: Eh7			! Valence band structure
real, dimension(l_k) :: Eh8			! Valence band structure
real, dimension(l_k) :: Eh9			! Valence band structure
real, dimension(l_k) :: Eh10			! Valence band structure
real, dimension(l_k) :: Eh11			! Valence band structure
real, dimension(l_k) :: Eh12			! Valence band structure

real, dimension(l_f) :: res
real, dimension(l_f) :: res1
real, dimension(l_f) :: res2
real, dimension(l_f) :: res3
real, dimension(l_f) :: res4
real, dimension(l_f) :: res5
real, dimension(l_f) :: res6
real, dimension(l_f) :: res7
real, dimension(l_f) :: res8
real, dimension(l_f) :: res9
real, dimension(l_f) :: res10
real, dimension(l_f) :: res11
real, dimension(l_f) :: res12
real, dimension(l_f) :: res13
real, dimension(l_f) :: res14
real, dimension(l_f) :: res15
real, dimension(l_f) :: res16
real, dimension(l_f) :: res17
real, dimension(l_f) :: res18
real, dimension(l_f) :: res19
real, dimension(l_f) :: res20
real, dimension(l_f) :: res21
real, dimension(l_f) :: res22

real, dimension(l_f) :: ps1
real, dimension(l_f) :: ps2
real, dimension(l_f) :: ps3
real, dimension(l_f) :: ps4
real, dimension(l_f) :: ps5
real, dimension(l_f) :: ps6
real, dimension(l_f) :: ps7
real, dimension(l_f) :: ps8
real, dimension(l_f) :: ps9
real, dimension(l_f) :: ps10
real, dimension(l_f) :: ps11
real, dimension(l_f) :: ps12
real, dimension(l_f) :: ps13
real, dimension(l_f) :: ps14
real, dimension(l_f) :: ps15
real, dimension(l_f) :: ps16
real, dimension(l_f) :: ps17
real, dimension(l_f) :: ps18
real, dimension(l_f) :: ps19
real, dimension(l_f) :: ps20
real, dimension(l_f) :: ps21
real, dimension(l_f) :: ps22
real, dimension(l_f) :: ps23
real, dimension(l_f) :: ps24
real, dimension(2,l_f) :: ps

real :: Ef_e
real :: Ef_h

real, dimension(l_k,l_k) :: V			! interaction matrix
real :: eps
integer :: info

real, dimension(l_k) :: ne				! Electrons' effective mass
real, dimension(l_k) :: nh				! Holes' effective mass
real, dimension(l_k) :: exce
real, dimension(l_k) :: exce1
integer, dimension(30) :: allowed

character (len=100) :: gway

! call saxpy(k,k,k,k)


!---------------------Band structure----------------------

Eg=1.52*e				!nominal band gap
me=0.02*m0				!electrons effective mass
mh=0.5*m0				!holes effective mass
eps=12.3				!permitivity

!----------------Formation of arrays---------------------

k_min=0.0					! min kx/ky
k_max=2.0e9					! max kx/ky
stk=(k_max-k_min)/l_k				! step in k-grid

k=(/(j*stk+k_min,j=0,int(l_k-1))/)		! k array

fff_min=-1.0*Eg/h				! min fr
fff_max=1.0*Eg/h				! max fr			
				
stf=(fff_max-fff_min)/(l_f)			! step in fr-grid

fff=(/ ((j*stf+fff_min),j=0,(int(l_f)-1)) /)	! fr array

!----------------------Plasma parameters-----------------

Tempr=300.0
conc=8.85e11

gway="/home/mk/results/InGaN_4nm_0d1_2003/dependencies_s/0_0_0d7_8d85/"

open(7,file=trim(gway)//'Ee1.dat')
read (7,'(e15.5)'), Ee1
 close(7)

Ee1=Ee1*e/h
!Ee1=Eg/h+h*(k**2.0)/(2.0*me)

open(7,file=trim(gway)//'Ee2.dat')
read (7,'(e15.5)'), Ee2
 close(7)

Ee2=Ee2*e/h
!Ee2=Eg/h+h*(k**2.0)/(2.0*me)

! 
! open(7,file='/home/mk/matlab/work/results/segregation/without/Ee3.dat')
! read (7,'(e15.5)'), Ee3
!  close(7)

!Ee3=Ee3*e/h
Ee3=Ee2
Ee4=Ee2
Ee5=Ee2
Ee6=Ee2
Ee7=Ee2
Ee8=Ee2
Ee9=Ee2

open(7,file=trim(gway)//'E1.dat')
read (7,'(e15.5)'), Eh1
close(7)

Eh1=Eh1*e/h
!Eh1=-h*(k**2.0)/(2.0*mh)

 open(7,file=trim(gway)//'E2.dat')
 read (7,'(e15.5)'), Eh2
 close(7)
 
 Eh2=Eh2*e/h
 ! Eh2=-0.1*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E3.dat')
 read (7,'(e15.5)'), Eh3
 close(7)
 
 Eh3=Eh3*e/h
 ! Eh3=-0.15*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E4.dat')
 read (7,'(e15.5)'), Eh4
 close(7)
 
 Eh4=Eh4*e/h
! ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)
! 
 open(7,file=trim(gway)//'E5.dat')
 read (7,'(e15.5)'), Eh5
 close(7)
 
 Eh5=Eh5*e/h
 ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E6.dat')
 read (7,'(e15.5)'), Eh6
 close(7)
 
 Eh6=Eh6*e/h
 ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E7.dat')
 read (7,'(e15.5)'), Eh7
 close(7)
 
 Eh7=Eh7*e/h
 ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E8.dat')
 read (7,'(e15.5)'), Eh8
 close(7)
 
 Eh8=Eh8*e/h
 ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)
 
 open(7,file=trim(gway)//'E9.dat')
 read (7,'(e15.5)'), Eh9
 close(7)
! 
! Eh9=Eh9*e/h
! ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)

 open(7,file=trim(gway)//'E10.dat')
 read (7,'(e15.5)'), Eh10
 close(7)
! 
! Eh9=Eh9*e/h
! ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)

 open(7,file=trim(gway)//'E11.dat')
 read (7,'(e15.5)'), Eh11
 close(7)
! 
! Eh9=Eh9*e/h
! ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)

 open(7,file=trim(gway)//'E12.dat')
 read (7,'(e15.5)'), Eh12
 close(7)
! 
! Eh9=Eh9*e/h
! ! Eh4=-0.2*e/h-h*(k**2.0)/(2.0*mh)


open(7,file=trim(gway)//'me11te.dat')
read (7,'(e17.5)'), mu11
close(7)

!mu11=1.0

 open(7,file=trim(gway)//'me21te.dat')
 read (7,'(e17.5)'), mu12
  close(7)
 
 open(7,file=trim(gway)//'me31te.dat')
 read (7,'(e17.5)'), mu13
  close(7)
 
 open(7,file=trim(gway)//'me41te.dat')
 read (7,'(e17.5)'), mu14
  close(7)
 
 open(7,file=trim(gway)//'me51te.dat')
 read (7,'(e17.5)'), mu15
  close(7)
 
 open(7,file=trim(gway)//'me61te.dat')
 read (7,'(e17.5)'), mu16
  close(7)
 
 open(7,file=trim(gway)//'me71te.dat')
 read (7,'(e17.5)'), mu17
  close(7)
 
 open(7,file=trim(gway)//'me81te.dat')
 read (7,'(e17.5)'), mu18
  close(7)
 
 open(7,file=trim(gway)//'me91te.dat')
 read (7,'(e17.5)'), mu19
  close(7)
  
 open(7,file=trim(gway)//'me101te.dat')
 read (7,'(e17.5)'), mu110
  close(7)
 
 open(7,file=trim(gway)//'me111te.dat')
 read (7,'(e17.5)'), mu111
  close(7)
 
 open(7,file=trim(gway)//'me121te.dat')
 read (7,'(e17.5)'), mu112
  close(7)  

 open(7,file=trim(gway)//'me12te.dat')
 read (7,'(e15.5)'), mu21
  close(7)
 
 open(7,file=trim(gway)//'me22te.dat')
 read (7,'(e15.5)'), mu22
  close(7)

  open(7,file=trim(gway)//'me32te.dat')
 read (7,'(e15.5)'), mu23
  close(7)
 
 open(7,file=trim(gway)//'me42te.dat')
 read (7,'(e15.5)'), mu24
  close(7) 
 
  open(7,file=trim(gway)//'me52te.dat')
 read (7,'(e15.5)'), mu25
  close(7)  
  
   open(7,file=trim(gway)//'me62te.dat')
 read (7,'(e15.5)'), mu26
  close(7) 
 
  open(7,file=trim(gway)//'me72te.dat')
 read (7,'(e15.5)'), mu27
  close(7)  
  
     open(7,file=trim(gway)//'me82te.dat')
 read (7,'(e15.5)'), mu28
  close(7) 
 
  open(7,file=trim(gway)//'me92te.dat')
 read (7,'(e15.5)'), mu29
  close(7)  
  
    open(7,file=trim(gway)//'me102te.dat')
 read (7,'(e15.5)'), mu210
  close(7)  
  
     open(7,file=trim(gway)//'me112te.dat')
 read (7,'(e15.5)'), mu211
  close(7) 
 
  open(7,file=trim(gway)//'me122te.dat')
 read (7,'(e15.5)'), mu212
  close(7) 
  
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me32te.dat')
! read (7,'(e15.5)'), mu23
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me42te.dat')
! read (7,'(e15.5)'), mu24
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me52te.dat')
! read (7,'(e15.5)'), mu25
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me62te.dat')
! read (7,'(e15.5)'), mu26
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me13te.dat')
! read (7,'(e15.5)'), mu31
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me23te.dat')
! read (7,'(e15.5)'), mu32
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me33te.dat')
! read (7,'(e15.5)'), mu33
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me43te.dat')
! read (7,'(e15.5)'), mu34
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me53te.dat')
! read (7,'(e15.5)'), mu35
!  close(7)
! 
! open(7,file='/home/mk/matlab/work/Новая папка (2)/me63te.dat')
! read (7,'(e15.5)'), mu36
!  close(7)

!call spotri( 'U', 1, k, 1, info )

!   call concentration(k,l_k,ne,conc)
  call Fermi_level(k,l_k,3,'c',Tempr,conc,Ef_e,h*Ee1,h*Ee2,h*Ee3,h*Ee4,h*Ee5,h*Ee6,h*Ee7,h*Ee8,h*Ee9)
  call Fermi_level(k,l_k,7,'v',Tempr,conc,Ef_h,h*Eh1,h*Eh2,h*Eh3,h*Eh4,h*Eh5,h*Eh6,h*Eh7,h*Eh8,h*Eh9)

!do j=2,l_k

!call Lind_formula(k(j),k,l_k,Eg,me,mh,Tempr,Ef_h,Ef_e,12.0,mu12(j))

!enddo

! open(7,file='eps.dat')
! write (7,'(e12.3)'), mu12
!close(7)

ne=1.0/(1+exp((Ee1*h-Ef_e)/kb/Tempr))
nh=1.0/(1+exp((Eh1*h-Ef_h)/kb/Tempr))

print *,'concentration=',ne(1)
print *,'concentration=',nh(1)
!print *,'Fermi level e=',Ef_e/e
!print *,'Fermi level h=',Ef_h/e

call int_matrix(k,l_k,mh,me,nh(1),ne(1),conc,Tempr,eps,11,gway,V)

!open(7,file='/home/mk/science/Fortran_programm/gain/src/v.dat')
!read (7,'(e17.5)'), V
!close(7)

!call exchange(k,l_k,ne,1.0-nh,V,exce)

if ((mu11(10))<(20.0*mu12(10))) then
	allowed(2)=1
else
	allowed(2)=0
endif
if ((mu11(10))<(20.0*mu13(10))) then
	allowed(3)=1
else 
	allowed(3)=0
endif
if ((mu11(1))<(20.0*mu14(10))) then
	allowed(4)=1
else
	allowed(4)=0
endif
if ((mu11(1))<(20.0*mu15(10))) then
	allowed(5)=1
else
	allowed(5)=0
endif
if ((mu11(10))<(20.0*mu16(10))) then
	allowed(6)=1
else
	allowed(6)=0
endif
if ((mu11(10))<(20.0*mu17(10))) then
	allowed(7)=1
else
	allowed(7)=0
endif
if ((mu11(10))<(20.0*mu18(10))) then
	allowed(8)=1
else
	allowed(8)=0
endif
if ((mu11(10))<(20.0*mu19(10))) then
	allowed(9)=1
else
	allowed(9)=0
endif
if ((mu11(10))<(20.0*mu110(10))) then
	allowed(10)=1
else
	allowed(10)=0
endif
if ((mu11(10))<(20.0*mu111(10))) then
	allowed(11)=1
else
	allowed(11)=0
endif
if ((mu11(10))<(20.0*mu112(10))) then
	allowed(12)=1
else
	allowed(12)=0
endif
if ((mu11(10))<(20.0*mu21(10))) then
	allowed(13)=1
else
	allowed(13)=0
endif
if ((mu11(10))<(20.0*mu22(10))) then
	allowed(14)=1
else
	allowed(14)=0
endif
if ((mu11(10))<(20.0*mu23(10))) then
	allowed(15)=1
else
	allowed(15)=0
endif
if ((mu11(10))<(20.0*mu24(10))) then
	allowed(16)=1
else
	allowed(16)=0
endif
if ((mu11(10))<(10.0*mu25(10))) then
	allowed(17)=1
else
	allowed(17)=0
endif
if ((mu11(10))<(20.0*mu26(10))) then
	allowed(18)=1
else
	allowed(18)=0
endif
if ((mu11(10))<(10.0*mu27(10))) then
	allowed(19)=1
else
	allowed(19)=0
endif
if ((mu11(10))<(20.0*mu28(10))) then
	allowed(20)=1
else
	allowed(20)=0
endif
if ((mu11(10))<(10.0*mu29(10))) then
	allowed(21)=1
else
	allowed(21)=0
endif
if ((mu11(10))<(10.0*mu210(10))) then
	allowed(22)=1
else
	allowed(22)=0
endif
if ((mu11(10))<(20.0*mu211(10))) then
	allowed(23)=1
else
	allowed(23)=0
endif
if ((mu11(10))<(10.0*mu212(10))) then
	allowed(24)=1
else
	allowed(24)=0
endif

allowed(1)=1

print *,'allowed=',allowed
!pause

call polarization1(1,k,l_k,Eh1,Ee1,mu11,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps1)  

 
 if (allowed(2)==1) then
	call polarization1(2,k,l_k,Eh2,Ee1,mu12,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps2) 
 endif
 if (allowed(3)==1) then
	call polarization1(3,k,l_k,Eh3,Ee1,mu13,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps3) 
 endif
 if (allowed(4)==1) then
	call polarization1(4,k,l_k,Eh4,Ee1,mu14,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps4) 
 endif
 if (allowed(5)==1) then
	call polarization1(5,k,l_k,Eh5,Ee1,mu15,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps5) 
 endif
 if (allowed(6)==1) then
	call polarization1(6,k,l_k,Eh6,Ee1,mu16,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps6) 
 endif
 if (allowed(7)==1) then
	call polarization1(7,k,l_k,Eh7,Ee1,mu17,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps7) 
endif
! if (allowed(8)==1) then
!	call polarization1(8,k,l_k,Eh8,Ee1,mu18,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps8) 
! endif
! if (allowed(9)==1) then
!	call polarization1(9,k,l_k,Eh9,Ee1,mu19,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps9) 
! endif
! if (allowed(10)==1) then
!	call polarization1(10,k,l_k,Eh10,Ee1,mu110,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps10) 
! endif
! if (allowed(11)==1) then
!	call polarization1(11,k,l_k,Eh11,Ee1,mu111,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps11) 
! endif
! if (allowed(12)==1) then
!	call polarization1(12,k,l_k,Eh12,Ee1,mu112,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps12) 
! endif
 
 
! if (allowed(13)==1) then
!	call polarization1(13,k,l_k,Eh1,Ee2,mu21,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps13) 
! endif
! if (allowed(14)==1) then
!	call polarization1(14,k,l_k,Eh2,Ee2,mu22,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps14) 
! endif
! if (allowed(15)==1) then
!	call polarization1(15,k,l_k,Eh3,Ee2,mu23,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps15) 
! endif
! if (allowed(16)==1) then
!	call polarization1(16,k,l_k,Eh4,Ee2,mu24,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps16) 
! endif
! if (allowed(17)==1) then
!	call polarization1(17,k,l_k,Eh5,Ee2,mu25,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps17) 
! endif
! if (allowed(18)==1) then
!	call polarization1(18,k,l_k,Eh6,Ee2,mu26,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps18) 
! endif
! if (allowed(19)==1) then
!	call polarization1(19,k,l_k,Eh7,Ee2,mu27,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps19) 
! endif
!  if (allowed(20)==1) then
!	call polarization1(20,k,l_k,Eh8,Ee2,mu28,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps20) 
! endif
! if (allowed(21)==1) then
!	call polarization1(21,k,l_k,Eh9,Ee2,mu29,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps21) 
! endif
!  if (allowed(22)==1) then
!	call polarization1(22,k,l_k,Eh10,Ee2,mu210,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps22) 
! endif
!  if (allowed(23)==1) then
!	call polarization1(23,k,l_k,Eh11,Ee2,mu211,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps23) 
! endif
! if (allowed(24)==1) then
!	call polarization1(24,k,l_k,Eh12,Ee2,mu212,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps24) 
! endif
 
 
! call polarization_normal(1,k,l_k,Eh1,Ee1,mu11,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps1)  
! call polarization_normal(2,k,l_k,Eh2,Ee1,mu12,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps2) 
! call polarization_normal(3,k,l_k,Eh3,Ee1,mu13,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps3) 
! call polarization_normal(4,k,l_k,Eh4,Ee1,mu14,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps4) 
! call polarization_normal(5,k,l_k,Eh5,Ee1,mu15,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps5) 
 
! call polarization1(k,l_k,Eh6,Ee1,mu16,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps6) 
! call polarization1(k,l_k,Eh7,Ee1,mu17,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps7) 
! call polarization1(k,l_k,Eh8,Ee1,mu18,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps8) 
! call polarization1(k,l_k,Eh9,Ee1,mu19,Ef_h,Ef_e,Tempr,conc,l_f,fff,V,exce,ps9)

! call polarization1(k,l_k,Eh1,Ee2,mu21,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps7) 
! call polarization1(k,l_k,Eh2,Ee2,mu22,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps8) 
! call polarization1(k,l_k,Eh3,Ee2,mu23,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps9) 
! call polarization1(k,l_k,Eh4,Ee2,mu24,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps10) 
! call polarization1(k,l_k,Eh5,Ee2,mu25,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps11) 
! call polarization1(k,l_k,Eh6,Ee2,mu26,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps12) 
! call polarization1(k,l_k,Eh1,Ee3,mu31,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps13) 
! call polarization1(k,l_k,Eh2,Ee3,mu32,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps14) 
! call polarization1(k,l_k,Eh3,Ee3,mu33,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps15) 
! call polarization1(k,l_k,Eh4,Ee3,mu34,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps16) 
! call polarization1(k,l_k,Eh5,Ee3,mu35,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps17) 
! call polarization1(k,l_k,Eh6,Ee3,mu36,Ef_h,Ef_e,Tempr,conc,l_f,fff,ps18) 

Eg=(Ee1(1)-Eh1(1))*h

!ps(1,:)=fff*h/(0.0042*e)
ps(1,:)=(fff*h+Eg)/e
! ps(1,:)=(fff*h+2.33*e)/e
! ps(2,:)=ps1

Eg11=Ee1(1)-Eh1(1)
Eg21=Ee1(1)-Eh2(1)
Eg31=Ee1(1)-Eh3(1)
Eg41=Ee1(1)-Eh4(1)
Eg51=Ee1(1)-Eh5(1)
Eg61=Ee1(1)-Eh6(1)
Eg71=Ee1(1)-Eh7(1)
Eg81=Ee1(1)-Eh8(1)
Eg91=Ee1(1)-Eh9(1)
Eg101=Ee1(1)-Eh10(1)
Eg111=Ee1(1)-Eh11(1)
Eg121=Ee1(1)-Eh12(1)


Eg12=Ee2(1)-Eh1(1)
Eg22=Ee2(1)-Eh2(1)
Eg32=Ee2(1)-Eh3(1)
Eg42=Ee2(1)-Eh4(1)
Eg52=Ee2(1)-Eh5(1)
Eg62=Ee2(1)-Eh6(1)
Eg72=Ee2(1)-Eh7(1)
Eg82=Ee2(1)-Eh8(1)
Eg92=Ee2(1)-Eh9(1)
Eg102=Ee2(1)-Eh10(1)
Eg112=Ee2(1)-Eh11(1)
Eg122=Ee2(1)-Eh12(1)

Eg13=Ee3(1)-Eh1(1)
Eg23=Ee3(1)-Eh2(1)
Eg33=Ee3(1)-Eh3(1)
Eg43=Ee3(1)-Eh4(1)
Eg53=Ee3(1)-Eh5(1)
Eg63=Ee3(1)-Eh6(1)
Eg73=Ee3(1)-Eh7(1)

res=abs(fff-fff(1)-(Eg21-Eg11))
res1=abs(fff-fff(1)-(Eg31-Eg11))
res2=abs(fff-fff(1)-(Eg41-Eg11))
res3=abs(fff-fff(1)-(Eg51-Eg11))
res4=abs(fff-fff(1)-(Eg61-Eg11))
res5=abs(fff-fff(1)-(Eg71-Eg11))
res6=abs(fff-fff(1)-(Eg81-Eg11))
res7=abs(fff-fff(1)-(Eg91-Eg11))
res8=abs(fff-fff(1)-(Eg101-Eg11))
res9=abs(fff-fff(1)-(Eg111-Eg11))
res10=abs(fff-fff(1)-(Eg121-Eg11))

res11=abs(fff-fff(1)-(Eg12-Eg11))
res12=abs(fff-fff(1)-(Eg22-Eg11))
res13=abs(fff-fff(1)-(Eg32-Eg11))
res14=abs(fff-fff(1)-(Eg42-Eg11))
res15=abs(fff-fff(1)-(Eg52-Eg11))
res16=abs(fff-fff(1)-(Eg62-Eg11))
res17=abs(fff-fff(1)-(Eg72-Eg11))
res18=abs(fff-fff(1)-(Eg82-Eg11))
res19=abs(fff-fff(1)-(Eg92-Eg11))
res20=abs(fff-fff(1)-(Eg102-Eg11))
res21=abs(fff-fff(1)-(Eg112-Eg11))
res22=abs(fff-fff(1)-(Eg122-Eg11))

!res13=abs(fff-fff(1)-(Eg33-Eg11))
!res14=abs(fff-fff(1)-(Eg43-Eg11))
!res15=abs(fff-fff(1)-(Eg53-Eg11))
!res16=abs(fff-fff(1)-(Eg63-Eg11))
 
j1=1.0
j2=1.0
j3=1.0
j4=1.0
j5=1.0
j6=1.0
j7=1.0
j8=1.0
j9=1.0
j10=1.0
j11=1.0
j12=1.0
j13=1.0
j14=1.0
j15=1.0
j16=1.0
j17=1.0
j18=1.0
j19=1.0
j20=1.0
j21=1.0
j22=1.0
j23=1.0

do j=2,l_f
	if (res(j)<res(j-1)) then
		j1=j
	endif
	if (res1(j)<res1(j-1)) then
		j2=j
	endif
	if (res2(j)<res2(j-1)) then
		j3=j
	endif
	if (res3(j)<res3(j-1)) then
		j4=j
	endif
	if (res4(j)<res4(j-1)) then
		j5=j
	endif
	if (res5(j)<res5(j-1)) then
		j6=j
	endif
	if (res6(j)<res6(j-1)) then
		j7=j
	endif
	if (res7(j)<res7(j-1)) then
		j8=j
	endif
	if (res8(j)<res8(j-1)) then
		j9=j
	endif
	if (res9(j)<res9(j-1)) then
		j10=j
	endif

	if (res10(j)<res10(j-1)) then
		j11=j
	endif
	if (res11(j)<res11(j-1)) then
		j12=j
	endif
	if (res12(j)<res12(j-1)) then
		j13=j
	endif
	if (res13(j)<res13(j-1)) then
		j14=j
	endif
	if (res14(j)<res14(j-1)) then
		j15=j
	endif
	if (res15(j)<res15(j-1)) then
		j16=j
	endif
	if (res16(j)<res16(j-1)) then
		j17=j
	endif
	
	if (res17(j)<res17(j-1)) then
		j18=j
	endif
	if (res18(j)<res18(j-1)) then
		j19=j
	endif
	if (res19(j)<res19(j-1)) then
		j20=j
	endif
	if (res20(j)<res20(j-1)) then
		j21=j
	endif
	if (res21(j)<res21(j-1)) then
		j22=j
	endif
	if (res22(j)<res22(j-1)) then
		j23=j
	endif
enddo
 
 !1.376900 
 !1.385000   
 !1.615330    
 !1.623560    
 
print *,j1
print *,j2
print *,j3
print *,j4
print *,j5
print *,j6
print *,j7
print *,j8
print *,j9
print *,j10
print *,j11
print *,j12
print *,j13
print *,j14
print *,j15
print *,j16
print *,j17
print *,j18
print *,j19
print *,j20
print *,j21
print *,j22
print *,j23
 
ps(2,:)=ps1
ps(2,j1:l_f)=ps(2,j1:l_f)+ps2(1:(l_f-j1+1))
ps(2,j2:l_f)=ps(2,j2:l_f)+ps3(1:(l_f-j2+1))
ps(2,j3:l_f)=ps(2,j3:l_f)+ps4(1:(l_f-j3+1))
ps(2,j4:l_f)=ps(2,j4:l_f)+ps5(1:(l_f-j4+1))
ps(2,j5:l_f)=ps(2,j5:l_f)+ps6(1:(l_f-j5+1))
ps(2,j6:l_f)=ps(2,j6:l_f)+ps7(1:(l_f-j6+1))
ps(2,j7:l_f)=ps(2,j7:l_f)+ps8(1:(l_f-j7+1))
ps(2,j8:l_f)=ps(2,j8:l_f)+ps9(1:(l_f-j8+1))
ps(2,j9:l_f)=ps(2,j9:l_f)+ps10(1:(l_f-j9+1))
ps(2,j10:l_f)=ps(2,j10:l_f)+ps11(1:(l_f-j10+1))
ps(2,j11:l_f)=ps(2,j11:l_f)+ps12(1:(l_f-j11+1))
ps(2,j12:l_f)=ps(2,j12:l_f)+ps13(1:(l_f-j12+1))
ps(2,j13:l_f)=ps(2,j13:l_f)+ps14(1:(l_f-j13+1))
ps(2,j14:l_f)=ps(2,j14:l_f)+ps15(1:(l_f-j14+1))
ps(2,j15:l_f)=ps(2,j15:l_f)+ps16(1:(l_f-j15+1))
ps(2,j16:l_f)=ps(2,j16:l_f)+ps17(1:(l_f-j16+1))
ps(2,j17:l_f)=ps(2,j17:l_f)+ps18(1:(l_f-j17+1))

ps(2,j18:l_f)=ps(2,j18:l_f)+ps19(1:(l_f-j18+1))
ps(2,j19:l_f)=ps(2,j19:l_f)+ps20(1:(l_f-j19+1))
ps(2,j20:l_f)=ps(2,j20:l_f)+ps21(1:(l_f-j20+1))
ps(2,j21:l_f)=ps(2,j21:l_f)+ps22(1:(l_f-j21+1))
ps(2,j22:l_f)=ps(2,j22:l_f)+ps23(1:(l_f-j22+1))
ps(2,j23:l_f)=ps(2,j23:l_f)+ps24(1:(l_f-j23+1))

 !ps(2,:)=ps1
 !ps(3,j1:l_f)=ps2(1:(l_f-j1+1))
 !ps(4,j2:l_f)=ps3(1:(l_f-j2+1))
 !ps(5,j3:l_f)=ps4(1:(l_f-j3+1))

ps(2,:)=ps(1,:)*ps(2,:)/(c*n_reff)

 open(7,file=trim(gway)//'ps4_n5m1ws.dat')
 write (7,'(2e15.5)'), ps	
close(7)

 open(7,file=trim(gway)//'ps_1.dat')
 write (7,'(e15.4)'), ps1
close(7)

 open(7,file=trim(gway)//'ps_2.dat')
 write (7,'(e15.4)'), ps2
close(7)

 open(7,file=trim(gway)//'ps_3.dat')
 write (7,'(e15.4)'), ps3
close(7)

 open(7,file=trim(gway)//'ps_4.dat')
 write (7,'(e15.4)'), ps4
close(7)

 open(7,file=trim(gway)//'ps_5.dat')
 write (7,'(e15.4)'), ps5
close(7)

 open(7,file=trim(gway)//'ps_6.dat')
 write (7,'(e15.4)'), ps6
close(7)

 open(7,file=trim(gway)//'ps_7.dat')
 write (7,'(e15.4)'), ps7
close(7)

 open(7,file=trim(gway)//'ps_8.dat')
 write (7,'(e15.4)'), ps8
close(7)

 open(7,file=trim(gway)//'Eg1.dat')
 write (7,'(e15.4)'), Eg/e
close(7)

 open(7,file=trim(gway)//'exce.dat')
 write (7,'(e15.4)'), exce
close(7)
 open(7,file=trim(gway)//'exce1.dat')
 write (7,'(e15.4)'), exce1
close(7)
 open(7,file=trim(gway)//'allowed.dat')
 write (7,'(e15.4)'), allowed
close(7)

end program gain

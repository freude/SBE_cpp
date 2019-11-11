!--------------------------------------------------------
!---------------Interaction matrix subroutine------------
!--------------------------------------------------------

subroutine int_matrix(k,l_k,mh,me,nh,ne,conc,Tempr,eps,ind,gway,V)

implicit none

!----------------Constantes------------------------------

integer, parameter						:: dp = selected_real_kind(15, 307)
real :: h							! Planck constant
real :: e							! elementary charge
real :: m0							! electron mass
real :: pi=3.14						! pi	
real :: eps0						! dielectric constant
complex :: ij						! imaginary unit

!----------------Material parameters---------------------

real :: epss						! dielectric constant
real :: eps							! Permittivity

integer :: j1
integer :: j2
integer :: j

integer :: l_k						! length of k array
real, dimension(l_k) :: k			! k-space (1st particle)
real, dimension(l_k) :: k1			! k-space (2nd particle)
real, dimension(l_k,l_k) :: V		! Coulomb potential
real, dimension(10002) :: phi
double precision, dimension(9) :: p	! k-space (1st particle)
double precision :: q	
double precision :: f_f
real :: Tempr
real :: me
real :: mh
real :: ne
real :: nh
real :: conc

integer :: ind
character (len=100) :: gway
real(dp), dimension(9) :: ccc

!----------------Definition of constants-----------------

h=1.05e-34
e=1.602e-19
m0=9.1e-31
eps0=8.85e-12
ij=sqrt(CMPLX(-1.0))

if (ind==11) then
open(7,file=trim(gway)//'f_factor.dat')
read (7,'(e15.5)'), p
close(7)
endif

!call screening(l_k,k,conc,gway,ccc)

phi=(/(j1*2*pi/10001,j1=0,(10001))/)

k1=k

do j1=1,l_k
    do j2=1,l_k  
        if ((j1==j2).OR.(j2==j1)) then
            V(j1,j2)=0.0
	    do j=2,10001
	                    
!		call dielectrical_const(q,mh,me,nh,ne,conc,Tempr,epss)	
       		
                q=sqrt(k(j1)**2+k1(j2)**2-2.0*k(j1)*k1(j2)*cos(phi(j)))	
				f_f=p(1)*(q**8.0)+p(2)*(q**7.0)+p(3)*(q**6.0)+p(4)*(q**5.0)+p(5)*(q**4.0)+p(6)*(q**3.0)+p(7)*(q**2.0)+p(8)*q+p(9)
				!epss=ccc(1)+ccc(2)*q+ccc(3)*(q**2.0)+ccc(4)*(q**3.0)+ccc(5)*(q**4.0)+ccc(6)*(q**5.0)+ccc(7)*(q**6.0)+ccc(8)*(q**7.0)+ccc(9)*(q**8.0)
                !f_f=1.0
				call dielectrical_const(q,mh,me,nh,ne,conc,Tempr,epss)	
				V(j1,j2)=V(j1,j2)+f_f*(epss*(e**2.0)/(8.0*(pi**2.0)*eps*eps0*q))*abs(phi(2)-phi(1))
       		
	    enddo
        else
	    V(j1,j2)=0.0
	    do j=1,10001
                q=sqrt(k(j1)**2+k1(j2)**2-2.0*k(j1)*k1(j2)*cos(phi(j)))
				f_f=p(1)*(q**8.0)+p(2)*(q**7.0)+p(3)*(q**6.0)+p(4)*(q**5.0)+p(5)*(q**4.0)+p(6)*(q**3.0)+p(7)*(q**2.0)+p(8)*q+p(9)
				!epss=ccc(1)+ccc(2)*q+ccc(3)*(q**2.0)+ccc(4)*(q**3.0)+ccc(5)*(q**4.0)+ccc(6)*(q**5.0)+ccc(7)*(q**6.0)+ccc(8)*(q**7.0)+ccc(9)*(q**8.0)
                !f_f=1.0
				call dielectrical_const(q,mh,me,nh,ne,conc,Tempr,epss)	
				V(j1,j2)=V(j1,j2)+f_f*(epss*(e**2.0)/(8.0*(pi**2.0)*eps*eps0*q))*abs(phi(2)-phi(1))
	    enddo
        endif
    enddo
enddo

V(1,1)=0.0

open(7,file='v.dat')
write (7,'(400e15.7)'), V

return
end subroutine int_matrix

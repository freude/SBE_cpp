 subroutine appr(x,l_x,y,M,coef)
! 
!  integer :: l_x
!  integer :: M
!  integer :: info
! real, dimension(M+1) :: ipiv
!  real, dimension(l_x) :: x
!  real, dimension(l_x) :: y
!  real, dimension(M+1) :: b
!  real, dimension(M+1) :: coef
!  real, dimension(l_x,M+1) :: A
!  real, dimension(M+1,M+1) :: AA
! 
! do j1=1,(M+1)
!     do j1=1,(M+1)
! 	A(j,j1)=x(j1)**(M+1-j)
!     enddo
! enddo
! 
! do j1=1,(M+1)
!     do j=1,(M+1)
! 	AA(j,j1)=sum(A(j,:)*A(j1,:))
!     enddo
!     b(j1)=sum(A(j1,:)*y)
! enddo
! 
! call sgetrf( M+1, M+1, AA, M+1, ipiv, info )
! call sgetrs( 'N', 1, M+1, AA, M+1, ipiv, b, M+1, info )
! 
!  coef=b
! 
 end subroutine appr
 

!--------------------------------------------------------
!---------------Electric field subroutine----------------
!--------------------------------------------------------

complex function E_field(ttt)

real, intent(in) :: ttt
real :: h				! Planck constant
real :: Eg
real :: stt
real :: omega
complex :: ij

h=1.05e-34
Eg=1.519*1.6e-19
ij=sqrt(CMPLX(-1.0))

!------Electrical field ha5e form of delta function-----

!stt=

!if (ttt<(3*stt)) then
!E_field=1.0
!else
!E_field=0.0
!endif

!------Electrical field have form of Gaussian function-----

stt=0.3e-14
omega=Eg/h
E_field=5.0e-5*exp(-((ttt-30*stt)**2)/(2*(stt**2)))!*exp(ij*(omega-Eg/h)*ttt)

return
end function E_field


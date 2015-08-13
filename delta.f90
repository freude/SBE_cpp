real function delta(x,eps)

implicit none

real, parameter :: pi=3.14
real x, eps


delta=1.0/(2.0*sqrt(pi*eps))*exp((-x**2)/(4*eps))


return
end function delta
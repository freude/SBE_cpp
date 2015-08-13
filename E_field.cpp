//--------------------------------------------------------
//---------------Electric field subroutine----------------
//--------------------------------------------------------

#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <cstdio>
#include <cstdlib>

complex<double> E_field(double ttt){

double h;				// Planck constant
double Eg;
double stt;
double omega;
const complex<double> ij(0,1);				// imaginary unit
complex<double> E_f;

h=1.05e-34;
Eg=1.519*1.6e-19;

//------Electrical field ha5e form of delta function-----

//stt=

//if (ttt<(3*stt)) then
//E_field=1.0
//else
//E_field=0.0
//endif

//------Electrical field have form of Gaussian function-----

stt=0.5e-14;
omega=Eg/h;
E_f=exp(-pow((ttt-10*stt),2)/(2*stt*stt))*exp(ij*(omega-Eg/h)*ttt);

return E_f;
}


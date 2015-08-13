#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <fstream>

#define l_t 1500

using namespace std;

void exchange(Vec_DP &k, Vec_DP &ne, Vec_DP &nh, Mat_DP &V, Vec_DP &exce);
complex<double> E_field(double ttt);

void polarization1(Vec_DP &k, Vec_DP &Eh, Vec_DP &Ee, Vec_DP &mu, double &Ef_h, double &Ef_e, double &Tempr, double &conc, Vec_DP &fff, Mat_DP &V, Vec_DP &PSr)
{

int l_k=k.size();				//length of k array
int l_f=fff.size();

//----------------Constantes------------------------------

const double h=1.05e-34;				// Planck constant
const double e=1.6e-19;				// elementary charge
const double m0=9.1e-31;				// electron mass
const double pi=3.14;				// pi	
const double eps0=8.85e-12;				// dielectric constant
const double kb=1.38e-23;				// Boltzmann constant
const double c=3e8;				// light speed in vacuum
const complex<double> ij(0,1);				// imaginary unit

//----------------Material parameters---------------------

double eps;					// Permittivity
double Eg;					// Band gap
//double Ef_e;					// Fermi level
//double Ef_h;					// Fermi level
double me;					// Electrons' effective mass
double mh;					// Holes' effective mass
//double Tempr;				// temperature
//double conc;					// concentration of carriers
//double Vol;					// active volume
double n_reff;				// refractive index

//---------------------------Counters----------------------------

int j1;
int j2;
int j;

//----------------------------time-------------------------------

double t_min;
double t_max;
double stt;

Vec_DP t(l_t);

//---------------------------k-space-----------------------------

double stk;				//

//Vec_DP k(l_k);		// k-space (1st particle)
Vec_DP k1(l_k);		// k-space (2nd particle)


//---------------------Fourier frequencies------------------------

//Vec_DP fff(l_f);

//----------------Single-particles characteristics---------------

//Vec_DP mu(l_k);		// dipole matrix element
Vec_DP ne(l_k);		// Electrons distribution function
Vec_DP nh(l_k);		// Holes distribution function
Vec_DP omega(l_k);	// Single-paticle transition frequency
//Vec_DP Ee(l_k);		// Conduction band structure
//Vec_DP Eh(l_k);		// Valence band structure

//---------------------Additional variables----------------------

complex<double> RS;

complex<double> kk1;				// Runge-Kutta RS
complex<double> kk2;				// Runge-Kutta RS
complex<double> kk3;				// Runge-Kutta RS
complex<double> kk4;				// Runge-Kutta RS

//-----------------Many-particle characteristics-----------------

Vec_DP exce(0.0,l_k);		// Exchange energy
double damp;				// dephasing
double vol;
//Mat_SP V(l_k,l_k);		// interaction matrix

//-------------------------Polarization--------------------------

Vec_CPLX_DP A(0.0,l_k);

Mat_CPLX_DP pp(0.0,l_t,l_k);
Vec_CPLX_DP P(0.0,l_t);
Vec_CPLX_DP PS(0.0,l_f);

//------------------------Electric field-------------------------

Vec_CPLX_DP E_ft(l_t);
Vec_CPLX_DP ES(l_f);

//------------------Material parameters-------------------

//eps=12.3					//permitivity
//n_reff=3.61				//refraction index
//Eg=2.5679*e				//nominal band gap
//me=0.148*m0				//electrons effective mass
//mh=0.234*m0				//holes effective mass

eps=12.3;					//permitivity
n_reff=3.61;				//refraction index
//Eg=1.519*e;				//nominal band gap
me=0.067*m0;				//electrons effective mass
mh=0.377*m0;				//holes effective mass

////----------------Formation of arrays---------------------

t_min=0.0;						// min time
t_max=0.5e-12;					// max time
stt=(t_max-t_min)/l_t;				// step in t-grid

for (j1=0;j1<l_t;j1++){
	t[j1]=t_min+j1*stt;			// time array
	}

stk=k[4]-k[3];

//--------------------------------------------------------

damp=0.015*e;					//damping

//----------------------Plasma parameters-----------------

//Ef_e=Eg+0.02*e
//Ef_h=-0.003*e
//Tempr=300
//conc=1.0e14
//conc=0.0


for (j1=0;j1<l_k;j1++)
{
//---------------------Transition frequency----------------------
	k1[j1]=k[j1];						// k' array
	omega[j1]=Ee[j1]-Eh[j1];

//-----------------Distribution functions--------------------

	ne[j1]=1.0/(1+exp((Ee[j1]*h-Ef_e)/kb/Tempr));
	nh[j1]=1.0/(1+exp((Eh[j1]*h-Ef_h)/kb/Tempr));
	//ne[j1]=0.0;
	//nh[j1]=1.0;
	pp[0][j1]=0.0;
	
}

Eg=h*omega[0];

//concentration(k,l_k,1.0-nh,conc1);
//print *,'concentration1=',conc1

//-----------------------------------------------------------

// call dielectrical_const(k,l_k,mh,me,0.0*minval(Ee),nh[1],ne[1],Tempr,conc,eps)

//----------------------Interaction matrix--------------------

//	call int_matrix(k,k1,l_k,eps,V,mh,me,Tempr,conc,ne[1],nh[1])
exchange(k,ne,nh,V,exce);

//-----------Solving semiconductor Bloch equations------------

//mu=mu*1.0e-27
vol=1.0;
ofstream outdata;
outdata.open("/home/mk/science/Cpp_programms/SBE_project/src/Es.dat",ios::in|ios::trunc);

FILE* Gplt = popen("gnuplot -persist","w");

for (j2=1;j2<l_t;j2++){
	for (j1=0;j1<l_k;j1++){	
					
		A[j1]=0.0;
		
		if ((j1>30)&&(j1<=(l_k-30))){
			for (j=(j1-30);j<(j1+30);j++){
				A[j1]+=V[j1][j]*pp[j2-1][j]*k[j]*stk;
			}
		}
		else if (j1<=30){
			for (j=0;j<(j1+30);j++){
				A[j1]+=V[j1][j]*pp[j2-1][j]*k[j]*stk;
			}
		}
		else if (j1>(l_k-30)){
			for (j=(j1-30);j<l_k;j++){
				A[j1]+=V[j1][j]*pp[j2-1][j]*k[j]*stk;
			}
		}

		RS=-ij*(omega[j1]-Eg/h-exce[j1]/h)*pp[j2-1][j1]-ij*(ne[j1]-nh[j1])*(mu[j1]*E_field(t[j2-1])+A[j1])/h-damp*pp[j2-1][j1]/h;
		
		kk1=RS;

		kk2=-ij*(omega[j1]-Eg/h-exce[j1]/h)*(pp[j2-1][j1]+stt*kk1/2.0)-ij*(ne[j1]-nh[j1])*(mu[j1]*E_field(t[j2-1]+stt/2)+A[j1])/h-damp*(pp[j2-1][j1]+stt*kk1/2.0)/h;
	
		kk3=-ij*(omega[j1]-Eg/h-exce[j1]/h)*(pp[j2-1][j1]+stt*kk2/2.0)-ij*(ne[j1]-nh[j1])*(mu[j1]*E_field(t[j2-1]+stt/2)+A[j1])/h-damp*(pp[j2-1][j1]+stt*kk2/2.0)/h;

		kk4=-ij*(omega[j1]-Eg/h-exce[j1]/h)*(pp[j2-1][j1]+stt*kk3)-ij*(ne[j1]-nh[j1])*(mu[j1]*E_field(t[j2-1]+stt)+A[j1])/h-damp*(pp[j2-1][j1]+stt*kk3)/h;

        pp[j2][j1]=pp[j2-1][j1]+(stt/6)*(kk1+2.0*kk2+2.0*kk3+kk4);
//        pp(j2,j1)=pp(j2-1,j1)+(stt)*kk1
		P[j2]+=2.0*pi/vol*mu[j1]*k[j1]*pp[j2][j1]*stk;
		cout << j2 << " " << "pp=" << pp[j2][j1] << " " << "ne=" << ne[j1] << " " << "nh=" << nh[j1] << " " << "exce=" << exce[j1] << " " << "A=" << A[j1] << "\n";
//		cout << mu[j1]<<ij<<"\n";
    }
	E_ft[j2]=E_field(t[j2-1]);
//	cout <<P[j2]<<"\n";
	outdata<<t[j2]<<" "<<real(P[j2])<<endl;
	
	if(j2==1){
		fprintf(Gplt,"plot '/home/mk/science/Cpp_programms/SBE_project/src/Es.dat' using 1:2 with lines \n");
		fflush (Gplt);
	}
	else if ((j2>1)&((j2%100)==0.0)){
		fprintf(Gplt,"reread; replot \n");
		fflush (Gplt);		
	}

}

outdata.close();
pclose(Gplt);

//-----------------Fourie transformation------------------

for (j=0;j<l_f;j++){
	for (j1=0;j1<l_t;j1++){
		ES[j]=ES[j]+E_ft[j1]*exp(ij*(fff[j]-Eg/h)*t[j1])*stt;
		PS[j]=PS[j]+P[j1]*exp(ij*(fff[j]-Eg/h)*t[j1])*stt/(4.0*pi*eps0*eps);
	}
	PSr[j]=(fff[j]+Eg/h)*imag(PS[j]/ES[j])/(c*n_reff);
}


}

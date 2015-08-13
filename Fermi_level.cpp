//--------------------------------------------------------
//-----------Fermi level computation subroutine-----------
//--------------------------------------------------------

#include <complex>
#include <iostream>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <cmath>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#define l_E 300

double Fermi_level(Vec_DP &k, int N_bands, char flag, double &Tempr, double &conc, Mat_DP &Ee1)
{

int l_k=k.size();

float kb;
float el;
float pi;
float h;
float m0;
float me;
float mh;
int j;				// counter
int j1;				// counter

float stk;				// length of k array
Vec_DP ne(l_k);			// k-space (1st particle)

float conc1;				// concentration
Vec_DP diff(2);				// concentration

double Ef;				// concentration
float Ef_min;				// concentration
float Ef_max;				// concentration
float st_Ef;				// concentration
Vec_DP Ef_pr(1000);		// concentration
Vec_DP a(1000);		// concentration

float stE;
float E_max;
float E_min;
Vec_DP dos_el(l_E);
Vec_DP dos_h(l_E);
Vec_DP E(l_E);
float damp;
float lz;

kb=1.38e-23;
el=1.6e-19;
h=1.05e-34;				//Planck constant
m0=9.1e-31;				//Electon mass
pi=3.14;

me=0.0665*m0;				//electrons effective mass
mh=0.234*m0;				//holes effective mass

E_min=-0.5;					// min kx/ky
E_max=0.5;					// max kx/ky
stE=(E_max-E_min)/l_E;				// step in k-grid
damp=0.005;

for (j=0;j<l_E;j++)
{
	E[j]=E_min+j*stE;
}

//--------------------------------------------------------

stk=k[3]-k[2];
lz=1.0;

//--------------------------------------------------------

diff[0]=conc*conc;

if (flag=='c'){
	Ef_min=-0.5;						// min kx/ky
	Ef_max=0.5;						// max kx/ky
	st_Ef=(Ef_max-Ef_min)/1000;				// step in k-grid
	
	for (j=0;j<l_E;j++)
	{
		Ef_pr[j]=h*Ee1[0][0]/el+Ef_min+j*st_Ef;
	}
	
	for (j=0;j<1000;j++){
		conc1=0.0;
		for (j1=0;j1<l_k;j1++){
			switch(N_bands){
				case 1 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;
				case 2 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 3 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 4 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 5 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 6 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 7 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 8 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[7][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;
				case 9 : ne[j1]=1/(1+exp((h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[7][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((h*Ee1[8][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				default:
				cout<< "Undefined number of subbands";
				break;			
			}
			conc1=conc1+2/lz*ne[j1]*k[j1]*stk/2/pi;
		}
		diff[1]=abs(conc-conc1);
		a[j]=diff[1];

		if (diff[0]>diff[1]){
			diff[0]=diff[1];
			Ef=Ef_pr[j];
		}
	}
}
else{
	Ef_min=-0.5;					// min kx/ky
	Ef_max=0.5;					// max kx/ky
	st_Ef=(Ef_max-Ef_min)/1000;			// step in k-grid
	
	for (j=0;j<l_E;j++)
	{
		Ef_pr[j]=-h*Ee1[0][l_k-1]/el+Ef_min+j*st_Ef;
	}
	
	for (j=0;j<1000;j++){
		conc1=0.0;
		for (j1=0;j1<l_k;j1++){
			switch(N_bands){
				case 1 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;			
				case 2 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;
				case 3 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 4 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 5 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 6 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 7 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 8 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[7][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				case 9 : ne[j1]=1/(1+exp((-h*Ee1[0][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[1][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[2][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[3][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[4][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[5][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[6][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[7][j1]-Ef_pr[j]*el)/kb/Tempr))+
								1/(1+exp((-h*Ee1[8][j1]-Ef_pr[j]*el)/kb/Tempr));
				break;						
				default:
				cout<< "Undefined number of subbands";
				break;			
			}
			conc1=conc1+2/lz*ne[j1]*k[j1]*stk/2/pi;
		}

		diff[1]=abs(conc-conc1);
		a[j]=diff[1];

		if (diff[0]>diff[1]){
			diff[0]=diff[1];
			Ef=-Ef_pr[j];
		}
	}
}

Ef=Ef*el;

return Ef;
}

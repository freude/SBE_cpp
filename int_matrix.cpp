//--------------------------------------------------------
//---------------Interaction matrix subroutine------------
//--------------------------------------------------------

#include <complex>
#include <iostream>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <cmath>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>

using namespace std;
using std::endl;
using std::ofstream;

void int_matrix(Vec_DP &k, double &mh, double &me, Vec_DP &nh, Vec_DP &ne, double &conc, double &Tempr, double &eps, int ind, string gway, Mat_DP &V){

//----------------Constantes------------------------------

//double h=1.05e-34;				// Planck constant
double e=1.602e-19;				// elementary charge
//double m0=9.1e-31;				// electron mass
double pi=3.14;				// pi	
double eps0=8.85e-12;				// dielectric constant
complex<double> ij(0,1);				// imaginary unit

//----------------Material parameters---------------------

double epss;				// dielectric constant

int j1;
int j2;
int j;

int l_k=k.size();		//length of k array
Vec_DP k1(k);		// k-space (2nd particle)
double phi[10002];
double q;

//----------------Definition of constants-----------------

if (ind==11){
	
	Vec_DP p(9);
	
	std::stringstream sstm;
	sstm << gway<<"f_factor.dat";	//valence band data reading
	string flnm=sstm.str();
	std::ifstream in(flnm.c_str(),std::ios::out);
	if(!in){  
		std::cout<<"Could not open file"<<sstm.str()<<std::endl;
		exit(0); 
	}
	j=0;
	while (!in.eof() && j<l_k) {
		in>>p[j];
		j++;
	}
	in.close();
}

for (j1=0;j1<3001;j1++){
	phi[j1]=j1*2*pi/3001;			// time array
	}
	
	ofstream outdata;
	outdata.open("/home/mk/my_projects/SBE_cpp/V.dat",ios::in|ios::trunc);
	

for (j1=0;j1<l_k;j1++){
    for (j2=0;j2<l_k;j2++){
        if ((j1==j2)||(j2==j1)){
            V[j1][j2]=0.0;
			for (j=1;j<3001;j++){
                q=sqrt(pow(k[j1],2)+pow(k1[j2],2)-2.0*k[j1]*k1[j2]*cos(phi[j]));	
				epss=1.0;
				V[j1][j2]=V[j1][j2]+(epss*e*e/(8.0*pow(pi,2)*eps*eps0*q))*fabs(phi[2]-phi[1]);
			}
		}
        else{
			V[j1][j2]=0.0;
			for (j=0;j<3001;j++){
                q=sqrt(pow(k[j1],2)+pow(k1[j2],2)-2.0*k[j1]*k1[j2]*cos(phi[j]));	
				epss=1.0;
				V[j1][j2]=V[j1][j2]+(epss*e*e/(8.0*pow(pi,2)*eps*eps0*q))*fabs(phi[2]-phi[1]);
			}
		}
		V[0][0]=0.0;
		outdata<<V[j1][j2]<<endl;
    }
}

outdata.close();

///* Open the file for writing. If it exists, append to it;
//otherwise, create a new file. */
//int fd = open ('/home/mk/science/Cpp_programms/SBE_project/src/V1.dat', O_WRONLY | O_CREAT | O_APPEND, 0666);
///* Write the timestamp to the file. */
//write(fd, V, l_k*l_k*sizeof(double));
///* All done. */
//close (fd);

//return 0;

}

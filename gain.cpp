#include <cstring>
#include <complex>
#include <iostream>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <cmath>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "gain.h"
#include <unistd.h>


using namespace std;
//using std::endl;
//using std::ofstream;

//----------------Matrix sizes----------------------------

#define l_k 400					//length of k array
#define l_f 9000				//length of frequency array!

const int Nsb_e=1;			//number of conduction subbands
const int Nsb_h=1;			//number of valence subbands


//----------------Definition of constants-----------------

const double h=1.054e-34;		//Planck constant
const double e=1.602e-19;		//elementary charge
const double m0=9.109e-31;		//electron mass
const double kb=1.38e-23;		//electron mass
const double pi=3.14;			//pi	
const double c=3.0e8;			//Light speed
const double n_reff=3.61;		//refraction index

//----------------Material parameters---------------------

double me;							//Electrons' effective mass
double mh;							//Holes' effective mass
double Tempr;						//temperature
double conc;						//concentration of carriers

//---------------------------k-space-----------------------------

double k_min;
double k_max;
double stk;
Vec_DP k(l_k);		// k-space (1st particle)
Vec_DP q(l_k-1);		// k-space (1st particle)

//---------------------Fourier frequencies------------------------

double fff_min;
double fff_max;
double stf;

Vec_DP fff(l_f);

//----------------Single-particles characteristics---------------

Mat_DP Eg(Nsb_e,Nsb_h);
double mu[Nsb_e][Nsb_h][l_k];
Mat_DP Ee(Nsb_e,l_k);
Mat_DP Eh(Nsb_h,l_k);

Vec_DP mu11(l_k);
Vec_DP Ee1(l_k);
Vec_DP Eh1(l_k);
 
Vec_DP ps1(l_f);
Vec_DP ps(l_f);

double Ef_e, Ef_h;

Mat_DP V(l_k,l_k);			//interaction matrix
double eps;
int info;

Vec_DP ne(l_k);				//Electrons' effective mass
Vec_DP nh(l_k);				//Holes' effective mass
Vec_DP exce(l_k);
Vec_DP exce1(l_k);
int allowed[30];

std::string gway;
std::stringstream sstm;

void gain()
{

gway = get_current_dir_name();
gway.append("/");
cout << gway << endl;	


//---------------------------Counters----------------------------

int j, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17, j18, j19, j20, j21, j22, j23;

//---------------------Band structure----------------------


Eg[0][0]=1.7*e;		//nominal band gap
me=0.02*m0;				//electrons effective mass
mh=0.5*m0;				//holes effective mass
eps=12.3;				//permitivity

//----------------Formation of arrays---------------------

k_min=0.0;					//min kx/ky
k_max=2.0e9;					//max kx/ky
stk=(k_max-k_min)/l_k;		//step in k-grid

for (j1=0;j1<l_k;j1++){
	k[j1]=k_min+j1*stk;   //k array
}

fff_min=0.001;				//min fr
fff_max=1.0*Eg[0][0]/h;				//max fr							
stf=(fff_max-fff_min)/(l_f);	//step in fr-grid

for (j1=0;j1<l_f;j1++){
	fff[j1]=fff_min+j1*stf+Eg[0][0]/h;			// fr array
	}

//--------------------------------------------------------------
//-----------------Reading data from files---------------------
//-----------------------BS and ME-----------------------------
//--------------------------------------------------------------

std::stringstream sstm;
string flnm;



for (j1=0;j1<Nsb_e;j1++){						//loop for conduction subbands
	//sstm << gway<<"Ee"<<(j1+1)<<".dat";			//conduction band data reading
	sstm << gway<<"Ee"<<".dat";			//conduction band data reading
	flnm=sstm.str();
	sstm.str(std::string());
	std::ifstream in(flnm.c_str(),std::ios::out);
	if(!in)
	{  
		std::cout<<"Could not open file1"<<sstm.str()<<std::endl;
		exit(0); 
	}
	j=0;
	while (!in.eof() && j<l_k) {
			in>>Ee[j1][j];
			Ee[j1][j]=Ee[j1][j];
			j++;
			}
    in.close();

	for (j2=0;j2<Nsb_h;j2++){					//loop for valence subbands
	
		if (j1==0){ 
			//sstm << gway<<"E"<<(j2+1)<<".dat";	//valence band data reading
			sstm << gway<<"Eh"<<".dat";	//valence band data reading
			flnm=sstm.str();
			sstm.str(std::string());
			std::ifstream in(flnm.c_str(),std::ios::out);
			if(!in)
			{  
				std::cout<<"Could not open file2"<<sstm.str()<<std::endl;
				exit(0); 
			}
			j=0;
			while (!in.eof() && j<l_k) {
				in>>Eh[j2][j];
				Eh[j2][j]=Eh[j2][j];
				j++;
			}
			in.close();
		}
		
		
		//sstm << gway<<"me"<<(j2+1)<<(j1+1)<<"te.dat";	//matrix element data reading
		sstm << gway<<"me"<<".dat";	//matrix element data reading
		flnm=sstm.str();
		sstm.str(std::string());
		cout<<flnm<<"\n";
		std::ifstream in(flnm.c_str(),std::ios::out);
		if(!in)
		{  
			std::cout<<"Could not open file3"<<sstm.str()<<std::endl;
			exit(0); 
		}
		j=0;
		while (!in.eof() && j<l_k) {
			in>>mu[j1][j2][j];
			j++;
		}
		in.close();
		
	}
}

ofstream outdata;
outdata.open(gway.append("ps.dat"),ios::in|ios::trunc);

for (j=0;j<l_k;j++){
outdata<<mu[0][0][j]<<endl;
}

outdata.close();

FILE* Gplt = popen("gnuplot -persist","w");
fprintf(Gplt,"plot '/home/mk/my_projects/SBE_cpp/ps.dat' with lines \n");	
fflush(Gplt);
pclose(Gplt);



//--------------------------------------------------------------
//----------------Computing many-body contributions-------------
//--------------------------------------------------------------

//------------------------Plasma parameters---------------------

Tempr=300.0;
conc=8.85e14;


//  call concentration(k,l_k,ne,conc)
Ef_e=Fermi_level(k, 1,'c',Tempr,conc,Ee);
Ef_h=Fermi_level(k, 1,'v',Tempr,conc,Eh);

cout<<Ef_e<<"\n";
cout<<Ef_h<<"\n";

//ne=1.0/(1+exp((Ee1*h-Ef_e)/kb/Tempr))
//nh=1.0/(1+exp((Eh1*h-Ef_h)/kb/Tempr))

//call int_matrix(k,l_k,mh,me,nh(1),ne(1),conc,Tempr,eps,11,gway,V)

int_matrix(k,mh,me,nh,ne,conc,Tempr,eps,10,gway,V);

//--------------------------------------------------------------
//----------------Computing the polarization--------------------
//--------------------------------------------------------------

for (j=0;j<l_f;j++) ps[j]=0.0;

for (j1=0;j1<Nsb_e;j1++){
	for (j2=0;j2<Nsb_h;j2++){
		
		if ((mu[1][1][10])<=(20.0*mu[j1][j2][10])){

			for (j=0;j<l_k;j++){
				Eh1[j]=Eh[j2][j];
				Ee1[j]=Ee[j1][j];
				mu11[j]=mu[j1][j2][j];
			}
			
			cout<<j1<<""<<j2<<"\n";
			polarization1(k,Eh1,Ee1,mu11,Ef_h,Ef_e,Tempr,conc,fff,V,ps1);
			
			for (j=0;j<l_f;j++){
				ps[j]=ps[j]+ps1[j];
			}	
											
		}		
	}
}

outdata.open("/home/mk/my_projects/SBE_cpp/ps.dat",ios::in|ios::trunc);

for (j=0;j<l_f;j++) outdata<<fff[j]*h/e<<" "<<ps[j]<<endl;

outdata.close();

Gplt = popen("gnuplot -persist","w");
fprintf(Gplt,"plot '/home/mk/my_projects/SBE_cpp/ps.dat' with lines \n");	
fflush(Gplt);
pclose(Gplt);

return;

}

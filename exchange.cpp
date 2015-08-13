//--------------------------------------------------------
//---------Exchange energy computation subroutine---------
//--------------------------------------------------------

#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <cmath>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>

using namespace std;

void exchange(Vec_DP &k, Vec_DP &ne, Vec_DP &nh, Mat_DP &V, Vec_DP &exce)
{

int j,j1;							// counter
int l_k=k.size();				// length of k array
double stk;						// length of k array
//const double pi=3.14;			// length of k array

stk=k[3]-k[2];

for (j=1;j<l_k;j++){
	exce[j]=0;
	for (j1=1;j1<l_k;j1++){
	exce[j]=exce[j]+V[j][j1]*(ne[j1]+nh[j1])*k[j1]*stk;
	}
}

}

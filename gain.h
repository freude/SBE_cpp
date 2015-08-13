#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <string>

#ifndef ADD_H_GUARD
#define ADD_H_GUARD

void polarization1(Vec_DP &k, Vec_DP &Eh, Vec_DP &Ee, Vec_DP &mu, double &Ef_h, double &Ef_e, double &Tempr, double &conc, Vec_DP &fff, Mat_DP &V, Vec_DP &ps1);

void int_matrix(Vec_DP &k, double &mh, double &me, Vec_DP &nh, Vec_DP &ne, double &conc, double &Tempr, double &eps, int ind, string gway, Mat_DP &V);

double Fermi_level(Vec_DP &k, int N_bands, char flag, double &Tempr, double &conc, Mat_DP &Ee);

#endif

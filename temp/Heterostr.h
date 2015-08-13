/*
 * Heterostr.h
 *
 *  Created on: Mar 12, 2012
 *      Author: mk
 */

#ifndef HETEROSTR_H_
#define HETEROSTR_H_

#include "TernaryAlloy.h"
#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include "BS.h"
#include "compar.h"

namespace std {

class Heterostr {
    string heterostr;
    TernaryAlloy* matarray;
    double* widtharray;
    size_t NumLayers;
    double Temper;
    double rat;
	int Het2Materials(const string);
public:
	Heterostr();
	Heterostr(string, Vec_DP&);
	Heterostr(string, Vec_DP&, double);
	virtual ~Heterostr();
	void PrintStr();
	Vec_DP getParam(Vec_DP&, string);
	Vec_DP pos_dep(Vec_DP&, Vec_DP&);
	void computeBS(BS&, string, Compar);
	void qwzb4x4ps(BS&, Compar);
};

} /* namespace std */
#endif /* HETEROSTR_H_ */

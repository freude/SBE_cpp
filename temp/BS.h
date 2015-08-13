/*
 * BS.h
 *
 *  Created on: Mar 18, 2012
 *      Author: mk
 */

#ifndef BS_H_
#define BS_H_

#include "nrutil_nr.h"
#include "nrtypes_nr.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>

namespace std {

class BS {
	string config;
	bool is_wfs_alloc;
public:
	Vec_DP kx;
	Vec_DP ky;
	Vec_DP x;
	Mat_DP bs;
	Mat_DP* wfs;
	size_t NumBands;
	BS();
	BS(string, double, double, double);
	BS(string, double, double, double, size_t);
	BS(string, double, double, double, size_t, Vec_DP&);
	void initia(size_t);
	void initia(size_t, Vec_DP&);
	virtual ~BS();

};

} /* namespace std */
#endif /* BS_H_ */

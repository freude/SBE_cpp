/*
 * BS.cpp
 *
 *  Created on: Mar 18, 2012
 *      Author: mk
 */

#include "BS.h"
#include <cmath>

namespace std {

BS::BS() {

	double k_min=0;
	double k_max=1;
	double step=1;

	size_t size=floor((k_max-k_min)/step);
	kx.vector=gsl_vector_alloc(size);
	kx.makeVar(k_min,k_max,step);
	config="one point";
	is_wfs_alloc=false;

	ky.vector=gsl_vector_alloc(1);
	ky(0)=0;

}

BS::BS(string appr, double k_min, double k_max, double step) {

	size_t size=floor((k_max-k_min)/step);
	kx.vector=gsl_vector_alloc(size);
	kx.makeVar(k_min,k_max,step);
	config=appr;
	is_wfs_alloc=false;

	ky.vector=gsl_vector_alloc(1);
	ky(0)=0;

}

BS::BS(string appr, double k_min, double k_max, double step, size_t numb) {

	size_t size=floor((k_max-k_min)/step);
	kx.vector=gsl_vector_alloc(size);
	kx.makeVar(k_min,k_max,step);
	config=appr;
	initia(numb);

	ky.vector=gsl_vector_alloc(1);
	ky(0)=0;

}

BS::BS(string appr, double k_min, double k_max, double step, size_t numb, Vector& coord) {

	size_t size=floor((k_max-k_min)/step);
	kx.vector=gsl_vector_alloc(size);
	kx.makeVar(k_min,k_max,step);
	config=appr;
	initia(numb,coord);

	ky.vector=gsl_vector_alloc(1);
	ky(0)=0;

}

BS::~BS() {
	if (is_wfs_alloc){
		delete [] wfs;
		gsl_vector_free(x.vector);
	}
	gsl_vector_free(kx.vector);
	gsl_matrix_free(bs.matrix);
}

void BS::initia(size_t numb){

	NumBands=numb;
	bs.matrix=gsl_matrix_calloc(NumBands,kx.vector->size);
	is_wfs_alloc=false;
    return;
}

void BS::initia(size_t numb, Vector& coord){

	NumBands=numb;
	bs.matrix=gsl_matrix_calloc(NumBands,kx.vector->size);
	x=coord;
	wfs=new Matrix [NumBands];
	is_wfs_alloc=true;
	for (size_t j=0; j<NumBands; j++) wfs[j].matrix=gsl_matrix_calloc(x.vector->size,kx.vector->size);
    return;
}

} /* namespace std */

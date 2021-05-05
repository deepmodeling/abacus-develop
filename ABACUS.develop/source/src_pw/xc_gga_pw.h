//==========================================================
// AUTHOR : mohan, Lixin He
// DATE : 2009-12-16
// UPDATE : Peize Lin change all double** to Vector3<double>* 2017-12-20
//==========================================================
#ifndef GGA_PW_H 
#define GGA_PW_H

#include "src_pw/tools.h"
#include "src_global/vector3.h"
#include "src_parallel/parallel_global.h"

//==========================================================
// Calculate the charge gradient using plane wave basis. 
//==========================================================

namespace GGA_PW
{
	void gradcorr(double &etxc, double &vtxc, matrix &v);
	void grad_rho( const complex<double> *rhog, Vector3<double> *gdr );
	void grad_dot( const Vector3<double> *h, double *dh);
	void noncolin_rho(double *rhoout1,double *rhoout2,double *seg);
}

#endif

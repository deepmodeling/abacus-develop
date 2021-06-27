//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-09-08
//==========================================================

#ifndef COMPLEXMATRIX_TEST_H
#define COMPLEXMATRIX_TEST_H

#include "module_base/complexmatrix.h"
#include<limits>

static std::ostream & operator<<( std::ostream & os, const ComplexMatrix & m )
{
	for( int ir=0; ir!=m.nr; ++ir )
	{
		for( int ic=0; ic!=m.nc; ++ic )
			if(std::norm(m(ir,ic))>1E-10)
			{
				if(std::imag(m(ir,ic))>1E-10)
					os<<m(ir,ic)<<"\t";
				else
					os<<std::real(m(ir,ic))<<"\t";
			}
			else
				os<<0<<"\t";
		os<<std::endl;
	}	
	return os;
}

#endif
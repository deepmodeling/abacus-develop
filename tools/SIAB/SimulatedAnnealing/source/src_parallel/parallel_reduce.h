#ifndef PARALLEL_REDUCE_H
#define PARALLEL_REDUCE_H

#ifdef __MPI
#include <mpi.h>
#endif

using namespace std;
#include <complex>

namespace Parallel_Reduce
{
	void reduce_int_all(int &object);
	void reduce_int_pool(int &object);
	
	// reduce double in all process
	void reduce_double_all(double &object);
	void reduce_double_all(double *object, const int n);

	// reduce double only in this pool
	// (each pool contain different k points)
	void reduce_double_pool(double &object);
	void reduce_double_pool(double *object, const int n);
	void reduce_double_allpool(double &object);

	void reduce_complex_double_all(complex<double> &object);
	void reduce_complex_double_pool(complex<double> &object);
	void reduce_complex_double_pool(complex<double> *object,const int n);

	void gather_min_int_all(int &v);
	void gather_max_double_all(double &v);
	void gather_min_double_all(double &v);
	void gather_max_double_pool(double &v);
	void gather_min_double_pool(double &v);

	bool check_if_equal(double &v);//mohan add 2009-11-11
}

#endif

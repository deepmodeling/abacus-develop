#ifndef PARALLEL_REDUCE_H
#define PARALLEL_REDUCE_H

#ifdef __MPI
#include <mpi.h>
#endif

#include <complex>
using std::complex;

namespace Parallel_Reduce
{
    /// reduce in all process
    template<typename T>
    void reduce_all(T& object);
    template<typename T>
    void reduce_all(T* object, const int n);
    template<typename T>
    void reduce_pool(T& object);
    template<typename T>
    void reduce_pool(T* object, const int n);

    void reduce_int_diag(int& object);//mohan add 2012-01-12

    void reduce_int_grid(int* object, const int n);//mohan add 2012-01-12

	// reduce double only in this pool
	// (each pool contain different k points)
	void reduce_double_grid(double *object, const int n);
	void reduce_double_diag(double *object, const int n);

	void reduce_double_allpool(double &object);
	void reduce_double_allpool(double *object, const int n);

	void gather_min_int_all(int &v);
	void gather_max_double_all(double &v);
	void gather_min_double_all(double &v);
	void gather_max_double_pool(double &v);
	void gather_min_double_pool(double &v);

	// mohan add 2011-04-21
	void gather_int_all(int &v, int *all);

	bool check_if_equal(double &v);//mohan add 2009-11-11
}

#endif

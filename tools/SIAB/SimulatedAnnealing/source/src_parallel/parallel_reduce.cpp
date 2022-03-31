#include "parallel_reduce.h"
#include "../src_spillage/tools.h"
#include "parallel_global.h"

void Parallel_Reduce::reduce_int_all(int &object)
{
#ifdef __MPI
	int swap = object;
	MPI_Allreduce(&swap , &object , 1, MPI_INT , MPI_SUM , MPI_COMM_WORLD);
#endif
	return;
}

void Parallel_Reduce::reduce_int_pool(int &object)
{
#ifdef __MPI
	int swap = object;
	MPI_Allreduce(&swap , &object , 1, MPI_INT , MPI_SUM , POOL_WORLD);
#endif
	return;
}



void Parallel_Reduce::reduce_double_all(double &object)
{
#ifdef __MPI
	double swap = object;
	MPI_Allreduce(&swap , &object , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
#endif
	return;
}

void Parallel_Reduce::reduce_double_all(double *object, const int n)
{
#ifdef __MPI
	double *swap = new double[n];
	for(int i=0;i<n;i++) swap[i] = object[i];
	MPI_Allreduce(swap, object, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	delete[] swap;
#endif
	return;
}


void Parallel_Reduce::reduce_double_pool(double &object)
{
#ifdef __MPI
	double swap = object;
	MPI_Allreduce(&swap , &object , 1, MPI_DOUBLE , MPI_SUM , POOL_WORLD);
#endif
	return;
}

void Parallel_Reduce::reduce_double_pool(double *object, const int n)
{
#ifdef __MPI
	double *swap = new double[n];
	for(int i=0;i<n;i++) swap[i] = object[i];
	MPI_Allreduce(swap, object, n, MPI_DOUBLE, MPI_SUM, POOL_WORLD);
	delete[] swap;
#endif
	return;
}


// (1) the value in same in each pool.
// (2) we need to reduce the value from different pool.
void Parallel_Reduce::reduce_double_allpool(double &object)
{
	if(KPAR==1) return;
#ifdef __MPI
	double swap = object / NPROC_IN_POOL;
	MPI_Allreduce(&swap , &object , 1, MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
#endif
}

void Parallel_Reduce::reduce_complex_double_all(complex<double> &object)
{
#ifdef __MPI
	complex<double> swap = object;
    MPI_Allreduce(&swap, &object, 1, mpicomplex, myOp, MPI_COMM_WORLD);
#endif
    return;
}

void Parallel_Reduce::reduce_complex_double_pool(complex<double> &object)
{
#ifdef __MPI
	complex<double> swap = object;
	MPI_Allreduce(&swap, &object, 1, mpicomplex, myOp, POOL_WORLD);
#endif
	return;
}

void Parallel_Reduce::reduce_complex_double_pool(complex <double> *object, const int n)
{
#ifdef __MPI
	complex<double> *swap = new complex<double>[n];
	for(int i=0;i<n;i++) swap[i] = object[i];
	MPI_Allreduce(swap, object, n, mpicomplex, myOp, POOL_WORLD);
	delete[] swap;
#endif
	return;
}

void Parallel_Reduce::gather_min_int_all(int &v)
{
#ifdef __MPI
	int *all = new int[NPROC];
	MPI_Allgather(&v, 1, MPI_INT, all, 1, MPI_INT, MPI_COMM_WORLD);
	for(int i=0; i<NPROC; i++)
	{
		if(v>all[i])
		{
			v = all[i];		
		}
	}
	delete[] all;
#endif
}

void Parallel_Reduce::gather_max_double_all(double &v)
{
#ifdef __MPI
	double *value=new double[NPROC];
	ZEROS(value, NPROC);
	MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	for(int i=0; i<NPROC; i++)
	{
		if(v<value[i])
		{
			v = value[i];
		}
	}
	delete[] value;
#endif
}

void Parallel_Reduce::gather_max_double_pool(double &v)
{
#ifdef __MPI
	if(NPROC_IN_POOL == 1) return;
	double *value=new double[NPROC_IN_POOL];
	ZEROS(value, NPROC_IN_POOL);
	MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, POOL_WORLD);
	for(int i=0; i<NPROC_IN_POOL; i++)
	{
		if(v<value[i])
		{
			v = value[i];
		}
	}
	delete[] value;
#endif
}

void Parallel_Reduce::gather_min_double_pool(double &v)
{
#ifdef __MPI
	if(NPROC_IN_POOL == 1) return;
	double *value=new double[NPROC_IN_POOL];
	ZEROS(value, NPROC_IN_POOL);
	MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, POOL_WORLD);
	for(int i=0; i<NPROC_IN_POOL; i++)
	{
		if(v>value[i])
		{
			v = value[i];
		}
	}
	delete[] value;
#endif
}

void Parallel_Reduce::gather_min_double_all(double &v)
{
#ifdef __MPI
	double *value=new double[NPROC];
	ZEROS(value, NPROC);
	MPI_Allgather(&v, 1, MPI_DOUBLE, value, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	for(int i=0; i<NPROC; i++)
	{
		if(v>value[i])
		{
			v = value[i];
		}
	}
	delete[] value;
#endif
}

bool Parallel_Reduce::check_if_equal(double &v)
{
#ifdef __MPI
	double *all=new double[NPROC];
	MPI_Allgather(&v, 1, MPI_DOUBLE, all, 1, MPI_DOUBLE, MPI_COMM_WORLD);
	for(int i=0; i<NPROC; i++)
	{
		if( abs(all[i] - all[0]) > 1.0e-9 )
		{
			for(int j=0; j<NPROC; j++)
			{
				cout << "\n processor = " << j << " value = " << all[j];
			}
			delete[] all;
			return 0;
		}
	}
	delete[] all;
	return 1;
#endif
}

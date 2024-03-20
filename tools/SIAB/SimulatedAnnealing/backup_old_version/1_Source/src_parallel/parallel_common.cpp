#include "parallel_common.h"
#include <string.h> //this is needed because strcpy

#ifdef __MPI
void Parallel_Common::bcast_string(string &object)
{
	char swap[100];
	if(MY_RANK == 0) strcpy(swap, object.c_str() );
	MPI_Bcast(swap, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
	if(MY_RANK != 0) object = static_cast<string>( swap );
//	cout<<"\n"<<object;
	return;
}

void Parallel_Common::bcast_string(string *object,const int n)
{
	for(int i=0;i<n;i++)
	{
		char swap[100];
		if(MY_RANK == 0) strcpy(swap, object[i].c_str() );
		MPI_Bcast(swap, 100, MPI_CHAR, 0, MPI_COMM_WORLD);
		if(MY_RANK != 0) object[i] = static_cast<string>( swap );
	}
	return;
}

void Parallel_Common::bcast_complex_double(complex<double> &object)
{
	double a = object.real();
	double b = object.imag();
	Parallel_Common::bcast_double(a);
	Parallel_Common::bcast_double(b);
	object = complex<double>( a, b);
	return;
}

void Parallel_Common::bcast_complex_double(complex<double> *object, const int n)
{
	double *a = new double[n];
	double *b = new double[n];

	for(int i=0; i<n; i++)
	{
		a[i] = object[i].real();
		b[i] = object[i].imag();
	}

	Parallel_Common::bcast_double( a, n);
	Parallel_Common::bcast_double( b, n);

	for(int i=0; i<n; i++)
	{
		object[i] = complex<double>( a[i], b[i]);
	}
	delete[] a;
	delete[] b;
}

void Parallel_Common::bcast_double(double &object)
{
	MPI_Bcast(&object, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	cout<<"\n"<<object;
	return;
}

void Parallel_Common::bcast_double(double *object,const int n)
{
	MPI_Bcast(object, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Parallel_Common::bcast_int(int &object)
{
	MPI_Bcast(&object, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	cout<<"\n"<<object;
	return;
}

void Parallel_Common::bcast_int(int *object,const int n)
{
	MPI_Bcast(object,n,MPI_INT,0,MPI_COMM_WORLD);
}

void Parallel_Common::bcast_bool(bool &object)
{
	int swap = object;
	if(MY_RANK == 0)swap = object;
    MPI_Bcast(&swap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(MY_RANK != 0)object = static_cast<bool>( swap );
//	cout<<"\n"<<object;
	return;
}

void Parallel_Common::bcast_char(char *object,const int n)
{
	MPI_Bcast(object, n, MPI_CHAR, 0, MPI_COMM_WORLD);
	return;
}



#endif

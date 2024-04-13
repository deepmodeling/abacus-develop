#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/complexarray.h"
#include "module_base/realarray.h"

class Parallel_Kpoints
{
	public:
	
	Parallel_Kpoints();
	~Parallel_Kpoints();

	void kinfo(int &nkstot);
	
	// collect value from each pool to wk.
	void pool_collection(double &value, const double *wk, const int &ik);

	// collect value from each pool to overlap.
	void pool_collection(double *valuea, double *valueb, const ModuleBase::realArray &a, const ModuleBase::realArray &b, const int &ik);
	void pool_collection(std::complex<double> *value, const ModuleBase::ComplexArray &w, const int &ik);

	// information about pool
	int *nproc_pool;
	int *startpro_pool;

	// inforamation about kpoints //qianrui add comment
	int* nks_pool; //number of k-points in each pool
	int* startk_pool; //the first k-point in each pool
	int* whichpool; //whichpool[k] : the pool which k belongs to

	private:

#ifdef __MPI
	void get_nks_pool(const int &nkstot);
	void get_startk_pool(const int &nkstot);
	void get_whichpool(const int &nkstot);
#endif


};

#endif

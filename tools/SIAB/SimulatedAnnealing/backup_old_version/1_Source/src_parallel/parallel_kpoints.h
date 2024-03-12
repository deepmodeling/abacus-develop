#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "../src_spillage/common.h"

class Parallel_Kpoints
{
	public:
	
	Parallel_Kpoints();
	~Parallel_Kpoints();	

	void init();
	void kinfo(int &nkstot);
	
	void pool_collection(double &value, const double *wk, const int &ik);
	void pool_collection(double *value, const realArray& overlap, const int &ik);

	// information about pool
	int *nproc_pool;
	int *startpro_pool;

	// inforamation about kpoints
	int* nks_pool;
	int* startk_pool;
	int* whichpool;

	private:

	void divide_pools();

	void get_nks_pool(const int &nkstot);
	void get_startk_pool(const int &nkstot);
	void get_whichpool(const int &nkstot);


};

#endif

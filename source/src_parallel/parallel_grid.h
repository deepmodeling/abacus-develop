#ifndef PARALLEL_GRID_H
#define PARALLEL_GRID_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

class Parallel_Grid
{
	public:

	friend class Efield;

	Parallel_Grid();
	~Parallel_Grid();
	
	void init(const int &ncx, const int &ncy, const int &ncz, 
		const int &nczp, const int &nrxx, const int &nbz, const int &bz);

	void init_final_scf(const int &ncx, const int &ncy, const int &ncz, 
		const int &nczp, const int &nrxx, const int &nbz, const int &bz); //LiuXh add 20180606

#ifdef __MPI	
	void zpiece_to_all(double *zpiece, const int &iz, double *rho);
	
	void reduce_to_fullrho(double *rhotot, double *rhoin);
#endif
	
	private:

	void z_distribution(void);
	
	int *nproc_in_pool;
	int **numz;
	int **startz;
	int **whichpro;
	int **numdata;
	int **startdata;

	int ncx;
	int ncy;
	int ncz;
	int ncxy;
	int ncxyz;
	int nczp; // different processors have different values.
	int nrxx;
	int nbz;
	int bz;

	bool allocate;
    bool allocate_final_scf; //LiuXh add 20180619
};

#endif

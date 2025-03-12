#ifndef PARALLEL_GRID_H
#define PARALLEL_GRID_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"

class Parallel_Grid
{
	public:

	friend class Efield;
	friend class Symmetry_rho;

    Parallel_Grid();
    Parallel_Grid(const int ncx_in, const int ncy_in, const int ncz_in, const int nczp_in, const int nrxx_in, const int nbz_in, const int bz_in)
        :ncx(ncx_in), ncy(ncy_in), ncz(ncz_in), nczp(nczp_in), nrxx(nrxx_in), nbz(nbz_in), bz(bz_in),
        ncxy(ncx_in* ncy_in), ncxyz(ncx_in* ncy_in* ncz_in)
    {
        assert(ncx > 0 && ncy > 0 && ncz > 0 && nczp >= 0 && nrxx > 0 && nbz > 0 && bz > 0);
    }
    ~Parallel_Grid();
	
	void init(const int &ncx, const int &ncy, const int &ncz, 
		const int &nczp, const int &nrxx, const int &nbz, const int &bz);

	void init_final_scf(const int &ncx, const int &ncy, const int &ncz, 
		const int &nczp, const int &nrxx, const int &nbz, const int &bz); //LiuXh add 20180606

#ifdef __MPI	
    void zpiece_to_all(double* zpiece, const int& iz, double* rho) const;
    void zpiece_to_stogroup(double* zpiece, const int& iz, double* rho) const; //qainrui add for sto-dft 2021-7-21

    /// @brief  Broadcast data from root to all processors. The index order is [x][y][z].
    void bcast(const double* const data_global, double* data_local, const int& rank)const;
    /// @brief  Reduce data from all processors to root. The index order is [x][y][z].
    void reduce(double* rhotot, const double* constrhoin)const;
#endif

    const int& nx = this->ncx;
    const int& ny = this->ncy;
    const int& nz = this->ncz;

	private:

	void z_distribution(void);
	
    int *nproc_in_pool = nullptr;
    int **numz = nullptr;
    int **startz = nullptr;
    int **whichpro = nullptr;

	int ncx=0;
	int ncy=0;
	int ncz=0;
	int ncxy=0;
	int ncxyz=0;
    int nczp=0; // number of z-layers (xy-planes) in each processor
	int nrxx=0;
	int nbz=0;
	int bz=0;

    bool allocate = false;
    bool allocate_final_scf = false; //LiuXh add 20180619
};

#endif

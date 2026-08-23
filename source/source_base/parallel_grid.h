#ifndef PARALLEL_GRID_H
#define PARALLEL_GRID_H

#include <cassert>
#include <vector>

class Parallel_Grid
{
    public:

    Parallel_Grid();
    Parallel_Grid(const int ncx_in, const int ncy_in, const int ncz_in, const int nczp_in, const int nrxx_in, const int nbz_in, const int bz_in)
        :ncx(ncx_in), ncy(ncy_in), ncz(ncz_in), nczp(nczp_in), nrxx(nrxx_in), nbz(nbz_in), bz(bz_in),
        ncxy(ncx_in* ncy_in), ncxyz(ncx_in* ncy_in* ncz_in)
    {
        assert(ncx > 0 && ncy > 0 && ncz > 0 && nczp >= 0 && nrxx > 0 && nbz > 0 && bz > 0);
    }
    ~Parallel_Grid();

    void init(const int &ncx, const int &ncy, const int &ncz,
        const int &nczp, const int &nrxx, const int &nbz, const int &bz,
        const int nprocgroup);

    /**
     * @brief Sum a distributed real-space grid across k-point pools.
     *
     * Uses a direct local-slab reduction for equal-sized pools. For uneven
     * pools, reconstructs a common global layout before the cross-pool sum.
     *
     * @param data Local real-space grid data ordered as [x][y][z].
     */
    void reduce_across_pools(double* data) const;

#ifdef __MPI
    /// @brief  Broadcast data from root to all processors. The index order is [x][y][z].
    void bcast(const double* const data_global, double* data_local, const int& rank, const bool is_sdft) const;
    /// @brief  Reduce data from all processors to root. The index order is [x][y][z].
    void reduce(double* rhotot, const double* constrhoin, const bool reduce_all_pool) const;
#endif

    int get_nx() const { return ncx; }
    int get_ny() const { return ncy; }
    int get_nz() const { return ncz; }

    private:

    void z_distribution(void);

#ifdef __MPI
    void zpiece_to_all(double* zpiece, const int& iz, double* rho) const;
    void zpiece_to_stogroup(double* zpiece, const int& iz, double* rho) const; //qainrui add for sto-dft 2021-7-21
#endif

    std::vector<int> nproc_in_pool;
    std::vector<std::vector<int>> numz;
    std::vector<std::vector<int>> startz;
    std::vector<std::vector<int>> whichpro;
    std::vector<std::vector<int>> whichpro_loc;

    int ncx=0;
    int ncy=0;
    int ncz=0;
    int ncxy=0;
    int ncxyz=0;
    int nczp=0; // number of z-layers (xy-planes) in each processor
    int nrxx=0;
    int nbz=0;
    int bz=0;
};

#endif

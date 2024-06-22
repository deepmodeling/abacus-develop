#include "local_orbital_wfc.h"

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/read_wfc_nao.h"

Local_Orbital_wfc::Local_Orbital_wfc()
{
    allocate_flag = false;
    wfck_flag = false;
    complex_flag = false;
    nks = 0;
}

Local_Orbital_wfc::~Local_Orbital_wfc()
{

    // used for k-points.
    if (this->complex_flag)
    {
        delete[] this->wfc_k_grid2;
    }
    if (this->wfck_flag)
    {
        for (int i = 0; i < nks; i++)
        {
            delete[] this->wfc_k_grid[i];
        }
        delete[] this->wfc_k_grid;
    }
}

int Local_Orbital_wfc::globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

int Local_Orbital_wfc::localIndex(int globalindex, int nblk, int nprocs, int& myproc)
{
    myproc = int((globalindex % (nblk * nprocs)) / nblk);
    return int(globalindex / (nblk * nprocs)) * nblk + globalindex % nblk;
}

#ifdef __MPI
void Local_Orbital_wfc::wfc_2d_to_grid(const double* wfc_2d,
                                       double** wfc_grid,
                                       const int ik,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg)
{
    ModuleBase::TITLE(" Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int myid = 0;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info = 0;

    // calculate maxnloc for bcasting 2d-wfc
    long maxnloc; // maximum number of elements in local matrix
    info = MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info = MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    std::vector<double> work(maxnloc); // work/buffer matrix

    int naroc[2]; // maximum number of row or column
    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work.data(), inc);
                naroc[0] = pv->nrow;
                naroc[1] = pv->ncol_bands;
            }
            info = MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info = MPI_Bcast(work.data(), maxnloc, MPI_DOUBLE, src_rank, pv->comm_2D);

            info = this->set_wfc_grid(naroc, pv->nb, pv->dim0, pv->dim1, iprow, ipcol, work.data(), wfc_grid);

        } // loop ipcol
    }     // loop iprow
    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");
}

void Local_Orbital_wfc::wfc_2d_to_grid(const std::complex<double>* wfc_2d,
                                       std::complex<double>** wfc_grid,
                                       int ik,
                                       const ModuleBase::matrix& ekb,
                                       const ModuleBase::matrix& wg,
                                       const std::vector<ModuleBase::Vector3<double>>& kvec_c)
{
    ModuleBase::TITLE("Local_Orbital_wfc", "wfc_2d_to_grid");
    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");

    const Parallel_Orbitals* pv = this->ParaV;
    const int inc = 1;
    int myid = 0;
    MPI_Comm_rank(pv->comm_2D, &myid);
    int info = 0;

    // calculate maxnloc for bcasting 2d-wfc
    long maxnloc = 0; // maximum number of elements in local matrix
    info = MPI_Reduce(&pv->nloc_wfc, &maxnloc, 1, MPI_LONG, MPI_MAX, 0, pv->comm_2D);
    info = MPI_Bcast(&maxnloc, 1, MPI_LONG, 0, pv->comm_2D);
    std::vector<std::complex<double>> work(maxnloc); // work/buffer matrix

    int naroc[2] = {0}; // maximum number of row or column
    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            const int coord[2] = {iprow, ipcol};
            int src_rank;
            info = MPI_Cart_rank(pv->comm_2D, coord, &src_rank);
            if (myid == src_rank)
            {
                BlasConnector::copy(pv->nloc_wfc, wfc_2d, inc, work.data(), inc);
                naroc[0] = pv->nrow;
                naroc[1] = pv->ncol_bands;
            }
            info = MPI_Bcast(naroc, 2, MPI_INT, src_rank, pv->comm_2D);
            info = MPI_Bcast(work.data(), maxnloc, MPI_DOUBLE_COMPLEX, src_rank, pv->comm_2D);
            // mohan update 2021-02-12, delte BFIELD option
            info = this->set_wfc_grid(naroc, pv->nb, pv->dim0, pv->dim1, iprow, ipcol, work.data(), wfc_grid);
        } // loop ipcol
    }     // loop iprow

    ModuleBase::timer::tick("Local_Orbital_wfc", "wfc_2d_to_grid");
}
#endif

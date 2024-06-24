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

void Local_Orbital_wfc::gamma_file(psi::Psi<double>* psid, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "gamma_file");
    std::cout << " Read in gamma point wave function files " << std::endl;

    double** ctot;

    // allocate psi
    int ncol = this->ParaV->ncol_bands;

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "lapack_gvx" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "cg_in_lcao"
#ifdef __CUDA
        || GlobalV::KS_SOLVER == "cusolver"
#endif
#ifdef __CUSOLVERMP
        || GlobalV::KS_SOLVER == "cusolvermp"
#endif
    )
    {
        ncol = this->ParaV->ncol;
    }

    if (psid == nullptr)
    {
        ModuleBase::WARNING_QUIT("gamma_file", "psid should be allocated first!");
    }
    else
    {
        psid->resize(GlobalV::NSPIN, ncol, this->ParaV->nrow);
    }
    ModuleBase::GlobalFunc::ZEROS(psid->get_pointer(), psid->size());

    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        this->error = ModuleIO::read_wfc_nao(ctot,
                                             is,
                                             GlobalV::GAMMA_ONLY_LOCAL,
                                             GlobalV::NB2D,
                                             GlobalV::NBANDS,
                                             GlobalV::NLOCAL,
                                             GlobalV::global_readin_dir,
                                             this->ParaV,
                                             psid,
                                             pelec);
#ifdef __MPI
        Parallel_Common::bcast_int(this->error);
#endif
        switch (this->error)
        {
        case 1:
            std::cout << "Can't find the wave function file: WFC_NAO_GAMMA" << is + 1 << ".txt" << std::endl;
            break;
        case 2:
            std::cout << "In wave function file, band number doesn't match" << std::endl;
            break;
        case 3:
            std::cout << "In wave function file, nlocal doesn't match" << std::endl;
            break;
        case 4:
            std::cout << "In k-dependent wave function file, k point is not correct" << std::endl;
            break;
        default:
            std::cout << " Successfully read in wave functions " << is << std::endl;
        }
        if (this->error)
        {
            ModuleBase::WARNING_QUIT("Local_Orbital_wfc::gamma_file", "Failed to read in wavefunction.");
        }
    } // loop ispin
}

void Local_Orbital_wfc::allocate_k(const int& lgd,
                                   psi::Psi<std::complex<double>>* psi,
                                   elecstate::ElecState* pelec,
                                   const int& nks,
                                   const int& nkstot,
                                   const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                                   const int& istep)
{
    this->nks = nks;

    ModuleBase::TITLE("Local_Orbital_wfc", "allocate_k");
    if (GlobalV::NLOCAL < GlobalV::NBANDS)
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate", "NLOCAL<NBANDS");
    }

    // mohan add the flag 2011-03-02
    // allocate the first part (only once!).
    if (this->wfck_flag == false)
    {
        this->wfc_k_grid = new std::complex<double>**[nks];
        for (int ik = 0; ik < nks; ik++)
        {
            this->wfc_k_grid[ik] = new std::complex<double>*[GlobalV::NBANDS];
        }
        this->wfck_flag = true;
    }

    if (this->complex_flag)
    {
        delete[] this->wfc_k_grid2;
        this->complex_flag = false;
    }
    // allocate the second part and initialize value as zero.
    // if(lgd != 0) xiaohui modify 2015-02-04, fixed memory bug
    if (lgd != 0)
    {
        const int page = GlobalV::NBANDS * lgd; // lgd: local grid dimension
        this->wfc_k_grid2
            = new std::complex<double>[nks
                                       * page]; // wfc_k_grid2 stores nks * nbands * lgd number of basis coefficients
        ModuleBase::GlobalFunc::ZEROS(wfc_k_grid2, nks * page);
        for (int ik = 0; ik < nks; ik++)
        {
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
            {
                this->wfc_k_grid[ik][ib]
                    = &wfc_k_grid2[ik * page + ib * lgd
                                   + 0]; // then wfc_k_grid stores the starting address of each band
                                         // but now there are less number of coefficients stored, if lgd < nbasis.
                // this is because in grid intergration, for each grid point, only few grid points (near some
                // neighboring atoms) would be involved in calculation, therefore only few basis coefficients are
                // needed.
            }
            ModuleBase::Memory::record("LOWF::wfc_k_grid",
                                       sizeof(std::complex<double>) * GlobalV::NBANDS * GlobalV::NLOCAL);
            this->complex_flag = true;
        }
    }

    // read wavefunction from file, then divide, distribute and broadcast
    if (INPUT.init_wfc == "file") // init_wfc can also be "atomic" but actually do nothing.
    {
        // confusing, seems istep to be the index of scf step. if not the first scf, why call this function?
        if (istep > 0)
        {
            return;
        }
        std::cout << " Read in wave functions files: " << nkstot << std::endl;
        if (psi == nullptr)
        {
            ModuleBase::WARNING_QUIT("allocate_k", "psi should be allocated first!");
        }
        else
        {
            psi->resize(nkstot, this->ParaV->ncol_bands, this->ParaV->nrow);
        }
        for (int ik = 0; ik < nkstot; ++ik)
        {
            std::complex<double>** ctot;
            this->error = ModuleIO::read_wfc_nao_complex(ctot,
                                                         ik,
                                                         GlobalV::NB2D,
                                                         GlobalV::NBANDS,
                                                         GlobalV::NLOCAL,
                                                         GlobalV::global_readin_dir,
                                                         kvec_c[ik],
                                                         this->ParaV,
                                                         psi,
                                                         pelec);
#ifdef __MPI
            Parallel_Common::bcast_int(this->error);
#endif
            switch (this->error)
            {
            case 1:
                std::cout << "Can't find the wave function file: WFC_NAO_K" << ik + 1 << ".txt" << std::endl;
                break;
            case 2:
                std::cout << "In wave function file, band number doesn't match" << std::endl;
                break;
            case 3:
                std::cout << "In wave function file, nlocal doesn't match" << std::endl;
                break;
            case 4:
                std::cout << "In k-dependent wave function file, k point is not correct" << std::endl;
                break;
            default:
                std::cout << " Successfully read in wave functions " << ik + 1 << std::endl;
            }
            if (this->error)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate_k", "Failed to read in wavefunction.");
            }
        }
    }

    return;
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

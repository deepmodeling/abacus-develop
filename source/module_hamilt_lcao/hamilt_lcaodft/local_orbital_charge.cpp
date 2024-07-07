#include "local_orbital_charge.h"

#include "module_base/blas_connector.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/read_wfc_nao.h"

// Shen Yu add 2019/5/9
extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int* npcol, int* myprow, int* mypcol);
    void Cblacs_pinfo(int* myid, int* nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int* prow, int* pcol);
    int Cblacs_pnum(int icontxt, int prow, int pcol);
}

int Local_Orbital_Charge::out_dm = 0;
int Local_Orbital_Charge::out_dm1 = 0;

Local_Orbital_Charge::Local_Orbital_Charge()
{
    // for gamma algorithms.
    this->init_DM = false;
    this->lgd_now = 0;
    this->lgd_last = 0;

    // for k-dependent algorithms.
    this->init_DM_R = false;

    // xiaohui add 2014-06-19
    // band_local = nullptr;
    // Z_wg = nullptr;
    // Z_LOC = nullptr;

    // move from Local_Orbital_Wfc
    allocate_flag = false;
    wfck_flag = false;
    complex_flag = false;
    nks = 0;
}

Local_Orbital_Charge::~Local_Orbital_Charge()
{
    // with gamma point only
    if (this->init_DM)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] DM[is];
            delete[] DM_pool[is];
        }
        delete[] DM;
        delete[] DM_pool;
    }

    // with k points
    if (this->init_DM_R)
    {
        for (int is = 0; is < GlobalV::NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
    }
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

void Local_Orbital_Charge::allocate_dm_wfc(const Grid_Technique& gt,
                                           elecstate::ElecState* pelec,
                                           psi::Psi<double>* psi,
                                           const K_Vectors& kv,
                                           const int& istep)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "allocate_dm_wfc");
    // here we reset the density matrix dimension.
    this->allocate_gamma(gt.lgd, gt.trace_lo, psi, pelec, kv.get_nks(), istep);
    return;
}

void Local_Orbital_Charge::allocate_dm_wfc(const Grid_Technique& gt,
                                           elecstate::ElecState* pelec,
                                           psi::Psi<std::complex<double>>* psi,
                                           const K_Vectors& kv,
                                           const int& istep)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "allocate_dm_wfc");
    // here we reset the density matrix dimension.
    this->allocate_k(gt.lgd, psi, pelec, kv.get_nks(), kv.get_nkstot(), kv.kvec_c, istep);
    this->allocate_DM_k(kv.get_nks(), gt.nnrg);
    return;
}

void Local_Orbital_Charge::set_dm_k(int ik, std::complex<double>* dm_k_in)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "set_dm_k");
    dm_k[ik].create(ParaV->ncol, ParaV->nrow);
    for (int i = 0; i < ParaV->ncol; ++i)
    {
        for (int j = 0; j < ParaV->nrow; ++j)
        {
            dm_k[ik](i, j) = dm_k_in[i * ParaV->nrow + j];
        }
    }
    return;
}

void Local_Orbital_Charge::set_dm_gamma(int is, double* dm_gamma_in)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "set_dm_gamma");
    dm_gamma[is].create(ParaV->ncol, ParaV->nrow);
    for (int i = 0; i < ParaV->ncol; ++i)
    {
        for (int j = 0; j < ParaV->nrow; ++j)
        {
            dm_gamma[is](i, j) = dm_gamma_in[i * ParaV->nrow + j];
        }
    }
    return;
}

void Local_Orbital_Charge::gamma_file(psi::Psi<double>* psid, elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("Local_Orbital_Charge", "gamma_file");
    std::cout << " Read in gamma point wave function files " << std::endl;

    double** ctot;

    // allocate psi
    int ncol = this->ParaV->ncol_bands;

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "lapack" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "cg_in_lcao"
#ifdef __CUDA
        || GlobalV::KS_SOLVER == "cusolver"
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

    int error;
    for (int is = 0; is < GlobalV::NSPIN; ++is)
    {
        error = ModuleIO::read_wfc_nao(ctot,
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
        Parallel_Common::bcast_int(error);
#endif
        switch (error)
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
        if (error)
        {
            ModuleBase::WARNING_QUIT("Local_Orbital_wfc::gamma_file", "Failed to read in wavefunction.");
        }
    } // loop ispin
}

void Local_Orbital_Charge::allocate_k(const int& lgd,
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
                this->wfc_k_grid[ik][ib] = &wfc_k_grid2[ik * page + ib * lgd + 0];
            }
            ModuleBase::Memory::record("LOC::wfc_k_grid",
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
        int error;
        for (int ik = 0; ik < nkstot; ++ik)
        {
            std::complex<double>** ctot;
            error = ModuleIO::read_wfc_nao_complex(ctot,
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
            Parallel_Common::bcast_int(error);
#endif
            switch (error)
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
            if (error)
            {
                ModuleBase::WARNING_QUIT("Local_Orbital_wfc::allocate_k", "Failed to read in wavefunction.");
            }
        }
    }

    return;
}
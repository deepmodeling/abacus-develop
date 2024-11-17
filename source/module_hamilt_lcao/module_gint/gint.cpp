#include "gint.h"

#include "module_parameter/parameter.h"
#if ((defined __CUDA))
#include "gint_force_gpu.h"
#include "gint_rho_gpu.h"
#include "gint_vl_gpu.h"
#endif

#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint::~Gint() {

    delete this->hRGint;
    delete this->hRGintCd;
    // in gamma_only case, DMRGint.size()=0, 
    // in multi-k case, DMRGint.size()=nspin
    for (int is = 0; is < this->DMRGint.size(); is++) {
        delete this->DMRGint[is];
    }
    for(int is = 0; is < this->hRGint_tmp.size(); is++) {
        delete this->hRGint_tmp[is];
    }
#ifdef __MPI
    delete this->DMRGint_full;
#endif
}

void Gint::cal_gint(Gint_inout* inout) {
    ModuleBase::TITLE("Gint_interface", "cal_gint");
    ModuleBase::timer::tick("Gint_interface", "cal_gint");
    // In multi-process environments,
    // some processes may not be allocated any data.
    if (this->gridt->get_init_malloced() == false) {
        ModuleBase::WARNING_QUIT("Gint_interface::cal_gint",
                                 "gridt has not been allocated yet!");
    }
    if (this->gridt->max_atom > 0) {
#ifdef __CUDA
        if (PARAM.inp.device == "gpu"
            && (inout->job == Gint_Tools::job_type::vlocal
                || inout->job == Gint_Tools::job_type::rho
                || inout->job == Gint_Tools::job_type::force)) {
            if (inout->job == Gint_Tools::job_type::vlocal) {
                gpu_vlocal_interface(inout);
            } else if (inout->job == Gint_Tools::job_type::rho) {
                gpu_rho_interface(inout);
            } else if (inout->job == Gint_Tools::job_type::force) {
                gpu_force_interface(inout);
            }
        } else
#endif
        {
#ifdef __MKL
            const int mkl_threads = mkl_get_max_threads();
            mkl_set_num_threads(mkl_threads);
#endif
            {
                if (inout->job == Gint_Tools::job_type::vlocal) {
                    gint_kernel_vlocal(inout);
                } else if (inout->job == Gint_Tools::job_type::dvlocal) {
                    gint_kernel_dvlocal(inout);
                } else if (inout->job == Gint_Tools::job_type::vlocal_meta) {
                    gint_kernel_vlocal_meta(inout);
                } else if (inout->job == Gint_Tools::job_type::rho) {
                    gint_kernel_rho(inout);
                } else if (inout->job == Gint_Tools::job_type::tau) {
                    gint_kernel_tau(inout);
                } else if (inout->job == Gint_Tools::job_type::force) {
                    gint_kernel_force(inout);
                } else if (inout->job == Gint_Tools::job_type::force_meta) {
                    gint_kernel_force_meta(inout);
                }
            }
        }
    }
    ModuleBase::timer::tick("Gint_interface", "cal_gint");
    return;
}
void Gint::prep_grid(const Grid_Technique& gt,
                     const int& nbx_in,
                     const int& nby_in,
                     const int& nbz_in,
                     const int& nbz_start_in,
                     const int& ncxyz_in,
                     const int& bx_in,
                     const int& by_in,
                     const int& bz_in,
                     const int& bxyz_in,
                     const int& nbxx_in,
                     const int& ny_in,
                     const int& nplane_in,
                     const int& startz_current_in,
                     const UnitCell* ucell_in,
                     const LCAO_Orbitals* orb_in) {
    ModuleBase::TITLE(GlobalV::ofs_running, "Gint_k", "prep_grid");

    this->gridt = &gt;
    this->nbx = nbx_in;
    this->nby = nby_in;
    this->nbz = nbz_in;
    this->ncxyz = ncxyz_in;
    this->nbz_start = nbz_start_in;
    this->bx = bx_in;
    this->by = by_in;
    this->bz = bz_in;
    this->bxyz = bxyz_in;
    this->nbxx = nbxx_in;
    this->ny = ny_in;
    this->nplane = nplane_in;
    this->startz_current = startz_current_in;
    this->ucell = ucell_in;
    assert(nbx > 0);
    assert(nby > 0);
    assert(nbz >= 0);
    assert(ncxyz > 0);
    assert(bx > 0);
    assert(by > 0);
    assert(bz > 0);
    assert(bxyz > 0);
    assert(nbxx >= 0);
    assert(ny > 0);
    assert(nplane >= 0);
    assert(startz_current >= 0);
    assert(this->ucell->omega > 0.0);

    return;
}

void Gint::initialize_pvpR(const UnitCell& ucell_in, Grid_Driver* gd, const int& nspin) {
    ModuleBase::TITLE("Gint", "initialize_pvpR");

    int npol = 1;
    // there is the only resize code of DMRGint
    if (this->DMRGint.size() == 0) {
        this->DMRGint.resize(nspin);
    }
    hRGint_tmp.resize(nspin);
    if (nspin != 4) {
        if (this->hRGint != nullptr) {
            delete this->hRGint;
        }
        this->hRGint = new hamilt::HContainer<double>(ucell_in.nat);
    } else {
        npol = 2;
        if (this->hRGintCd != nullptr) {
            delete this->hRGintCd;
        }
        this->hRGintCd
            = new hamilt::HContainer<std::complex<double>>(ucell_in.nat);
        for (int is = 0; is < nspin; is++) {
            if (this->DMRGint[is] != nullptr) {
                delete this->DMRGint[is];
            }
            if (this->hRGint_tmp[is] != nullptr) {
                delete this->hRGint_tmp[is];
            }
            this->DMRGint[is] = new hamilt::HContainer<double>(ucell_in.nat);
            this->hRGint_tmp[is] = new hamilt::HContainer<double>(ucell_in.nat);
        }
#ifdef __MPI
        if (this->DMRGint_full != nullptr) {
            delete this->DMRGint_full;
        }
        this->DMRGint_full = new hamilt::HContainer<double>(ucell_in.nat);
#endif
    }

    if (PARAM.globalv.gamma_only_local && nspin != 4) {
        this->hRGint->fix_gamma();
    }
    if (npol == 1) {
        this->hRGint->insert_ijrs(this->gridt->get_ijr_info(), ucell_in);
        this->hRGint->allocate(nullptr, true);
        ModuleBase::Memory::record("Gint::hRGint",
                            this->hRGint->get_memory_size());
        // initialize DMRGint with hRGint when NSPIN != 4
        for (int is = 0; is < this->DMRGint.size(); is++) {
            if (this->DMRGint[is] != nullptr) {
                delete this->DMRGint[is];
            }
            this->DMRGint[is] = new hamilt::HContainer<double>(*this->hRGint);
        }
        ModuleBase::Memory::record("Gint::DMRGint",
                                   this->DMRGint[0]->get_memory_size()
                                       * this->DMRGint.size());
    } else {
        this->hRGintCd->insert_ijrs(this->gridt->get_ijr_info(), ucell_in, npol);
        this->hRGintCd->allocate(nullptr, true);
        for(int is = 0; is < nspin; is++) {
            this->hRGint_tmp[is]->insert_ijrs(this->gridt->get_ijr_info(), ucell_in);
            this->DMRGint[is]->insert_ijrs(this->gridt->get_ijr_info(), ucell_in);
            this->hRGint_tmp[is]->allocate(nullptr, true);
            this->DMRGint[is]->allocate(nullptr, true);
        }
        ModuleBase::Memory::record("Gint::hRGint_tmp",
                                       this->hRGint_tmp[0]->get_memory_size()*nspin);
        ModuleBase::Memory::record("Gint::DMRGint",
                                       this->DMRGint[0]->get_memory_size()
                                           * this->DMRGint.size()*nspin);
#ifdef __MPI
        this->DMRGint_full->insert_ijrs(this->gridt->get_ijr_info(), ucell_in, npol);
        this->DMRGint_full->allocate(nullptr, true);
        ModuleBase::Memory::record("Gint::DMRGint_full",
                                   this->DMRGint_full->get_memory_size());
#endif
    }
}

void Gint::reset_DMRGint(const int& nspin)
{
    if (this->hRGint)
    {
        for (auto& d : this->DMRGint) { delete d; }
        this->DMRGint.resize(nspin);
        this->DMRGint.shrink_to_fit();
        for (auto& d : this->DMRGint) { d = new hamilt::HContainer<double>(*this->hRGint); }
        if (nspin == 4)
        {
            for (auto& d : this->DMRGint) { d->allocate(nullptr, false); }
#ifdef __MPI
            delete this->DMRGint_full;
            this->DMRGint_full = new hamilt::HContainer<double>(*this->hRGint);
            this->DMRGint_full->allocate(nullptr, false);
#endif
        }
    }
}

void Gint::transfer_DM2DtoGrid(std::vector<hamilt::HContainer<double>*> DM2D) {
    ModuleBase::TITLE("Gint", "transfer_DMR");

    // To check whether input parameter DM2D has been initialized
#ifdef __DEBUG
    assert(!DM2D.empty()
           && "Input parameter DM2D has not been initialized while calling "
              "function transfer_DM2DtoGrid!");
#endif

    ModuleBase::timer::tick("Gint", "transfer_DMR");
    if (PARAM.inp.nspin != 4) {
        for (int is = 0; is < this->DMRGint.size(); is++) {
#ifdef __MPI
            hamilt::transferParallels2Serials(*DM2D[is], DMRGint[is]);
#else
            this->DMRGint[is]->set_zero();
            this->DMRGint[is]->add(*DM2D[is]);
#endif
        }
    } else // NSPIN=4 case
    {
#ifdef __MPI
        hamilt::transferParallels2Serials(*DM2D[0], this->DMRGint_full);
#else
        this->DMRGint_full = DM2D[0];
#endif
        std::vector<double*> tmp_pointer(4, nullptr);
        for (int iap = 0; iap < this->DMRGint_full->size_atom_pairs(); ++iap) {
            auto& ap = this->DMRGint_full->get_atom_pair(iap);
            int iat1 = ap.get_atom_i();
            int iat2 = ap.get_atom_j();
            for (int ir = 0; ir < ap.get_R_size(); ++ir) {
                const ModuleBase::Vector3<int> r_index = ap.get_R_index(ir);
                for (int is = 0; is < 4; is++) {
                    tmp_pointer[is] = this->DMRGint[is]
                                          ->find_matrix(iat1, iat2, r_index)
                                          ->get_pointer();
                }
                double* data_full = ap.get_pointer(ir);
                for (int irow = 0; irow < ap.get_row_size(); irow += 2) {
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2) {
                        *(tmp_pointer[0])++ = data_full[icol];
                        *(tmp_pointer[1])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                    for (int icol = 0; icol < ap.get_col_size(); icol += 2) {
                        *(tmp_pointer[2])++ = data_full[icol];
                        *(tmp_pointer[3])++ = data_full[icol + 1];
                    }
                    data_full += ap.get_col_size();
                }
            }
        }
    }
    ModuleBase::timer::tick("Gint", "transfer_DMR");
}
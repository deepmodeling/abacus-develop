#include "gint.h"

#include "source_io/module_parameter/parameter.h"
#if ((defined __CUDA))
#include "gint_force_gpu.h"
#include "gint_rho_gpu.h"
#include "gint_vl_gpu.h"
#endif

#include "source_base/memory.h"
#include "source_base/timer.h"
#include "source_basis/module_ao/ORB_read.h"
#include "source_lcao/module_hcontainer/hcontainer_funcs.h"
#include "source_pw/module_pwdft/global.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint::~Gint() {

    delete this->hRGint;
    delete this->hRGintCd;
    // in gamma_only case, dmr_gint.size()=0, 
    // in multi-k case, dmr_gint.size()=nspin
    for (int is = 0; is < this->dmr_gint.size(); is++) {
        delete this->dmr_gint[is];
    }
    for(int is = 0; is < this->hr_gint_tmp .size(); is++) {
        delete this->hr_gint_tmp [is];
    }
#ifdef __MPI
    delete this->dm2d_tmp;
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

void Gint::initialize_pvpR(const UnitCell& ucell_in, const Grid_Driver* gd, const int& nspin)
{
    ModuleBase::TITLE("Gint", "initialize_pvpR");
    int npol = 1;
    // there is the only resize code of dmr_gint
    if (this->dmr_gint.size() == 0) {
        this->dmr_gint.resize(nspin);
    }
    hr_gint_tmp.resize(nspin);
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
            if (this->dmr_gint[is] != nullptr) {
                delete this->dmr_gint[is];
            }
            if (this->hr_gint_tmp[is] != nullptr) {
                delete this->hr_gint_tmp[is];
            }
            this->dmr_gint[is] = new hamilt::HContainer<double>(ucell_in.nat);
            this->hr_gint_tmp[is] = new hamilt::HContainer<double>(ucell_in.nat);
        }
#ifdef __MPI
        if (this->dm2d_tmp != nullptr) {
            delete this->dm2d_tmp;
        }
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
        // initialize dmr_gint with hRGint when NSPIN != 4
        for (int is = 0; is < this->dmr_gint.size(); is++) {
            if (this->dmr_gint[is] != nullptr) {
                delete this->dmr_gint[is];
            }
            this->dmr_gint[is] = new hamilt::HContainer<double>(*this->hRGint);
        }
        ModuleBase::Memory::record("Gint::dmr_gint",
                                   this->dmr_gint[0]->get_memory_size()
                                       * this->dmr_gint.size());
    } else {
        this->hRGintCd->insert_ijrs(this->gridt->get_ijr_info(), ucell_in, npol);
        this->hRGintCd->allocate(nullptr, true);
        for(int is = 0; is < nspin; is++) {
            this->hr_gint_tmp[is]->insert_ijrs(this->gridt->get_ijr_info(), ucell_in);
            this->dmr_gint[is]->insert_ijrs(this->gridt->get_ijr_info(), ucell_in);
            this->hr_gint_tmp[is]->allocate(nullptr, true);
            this->dmr_gint[is]->allocate(nullptr, true);
        }
        ModuleBase::Memory::record("Gint::hr_gint_tmp",
                                       this->hr_gint_tmp[0]->get_memory_size()*nspin);
        ModuleBase::Memory::record("Gint::dmr_gint",
                                       this->dmr_gint[0]->get_memory_size()
                                           * this->dmr_gint.size()*nspin);
    }
}

void Gint::reset_DMRGint(const int& nspin)
{
    if (this->hRGint)
    {
        for (auto& d : this->dmr_gint) { delete d; }
        this->dmr_gint.resize(nspin);
        this->dmr_gint.shrink_to_fit();
        for (auto& d : this->dmr_gint) { d = new hamilt::HContainer<double>(*this->hRGint); }
        if (nspin == 4)
        {
            for (auto& d : this->dmr_gint) { d->allocate(nullptr, false); }
#ifdef __MPI
            delete this->dm2d_tmp;
#endif
        }
    }
}

void Gint::transfer_DM2DtoGrid(std::vector<hamilt::HContainer<double>*> dm2d) {
    ModuleBase::TITLE("Gint", "transfer_DMR");
    // To check whether input parameter dm2d has been initialized
#ifdef __DEBUG
    assert(!dm2d.empty()
           && "Input parameter dm2d has not been initialized while calling "
              "function transfer_DM2DtoGrid!");
#endif
    ModuleBase::timer::tick("Gint", "transfer_DMR");
    if (PARAM.inp.nspin != 4) {
        for (int is = 0; is < this->dmr_gint.size(); is++) {
#ifdef __MPI
            hamilt::transferParallels2Serials(*dm2d[is], dmr_gint[is]);
#else
            this->dmr_gint[is]->set_zero();
            this->dmr_gint[is]->add(*dm2d[is]);
#endif
        }
    } else // NSPIN=4 case
    {
        // is=0:↑↑, 1:↑↓, 2:↓↑, 3:↓↓
        const int row_set[4] = {0, 0, 1, 1};
        const int col_set[4] = {0, 1, 0, 1};
        int mg = dm2d[0]->get_paraV()->get_global_row_size()/2;
        int ng = dm2d[0]->get_paraV()->get_global_col_size()/2;
        int nb = dm2d[0]->get_paraV()->get_block_size()/2;
        auto ijr_info = dm2d[0]->get_ijr_info();
#ifdef __MPI
        int blacs_ctxt = dm2d[0]->get_paraV()->blacs_ctxt;
        std::vector<int> iat2iwt(ucell->nat);
        for (int iat = 0; iat < ucell->nat; iat++) {
            iat2iwt[iat] = ucell->get_iat2iwt()[iat]/2;
        }
        Parallel_Orbitals pv{};
        pv.set(mg, ng, nb, blacs_ctxt);
        pv.set_atomic_trace(iat2iwt.data(), ucell->nat, mg);
        this-> dm2d_tmp = new hamilt::HContainer<double>(&pv, nullptr, &ijr_info);
#else
        if (this->dm2d_tmp != nullptr) {
            delete this->dm2d_tmp;
        }
        this-> dm2d_tmp = new hamilt::HContainer<double>(*this->hRGint);
        this-> dm2d_tmp -> insert_ijrs(this->gridt->get_ijr_info(), *(this->ucell));
        this-> dm2d_tmp -> allocate(nullptr, true);
#endif
        ModuleBase::Memory::record("Gint::dm2d_tmp", this->dm2d_tmp->get_memory_size());
        for (int is = 0; is < 4; is++){
            for (int iap = 0; iap < dm2d[0]->size_atom_pairs(); ++iap) {
                auto& ap = dm2d[0]->get_atom_pair(iap);
                int iat1 = ap.get_atom_i();
                int iat2 = ap.get_atom_j();
                for (int ir = 0; ir < ap.get_R_size(); ++ir) {
                    const ModuleBase::Vector3<int> r_index = ap.get_R_index(ir);
                    double* matrix_out = this -> dm2d_tmp -> find_matrix(iat1, iat2, r_index)->get_pointer();
                    double* matrix_in = ap.get_pointer(ir);
                    for (int irow = 0; irow < ap.get_row_size()/2; irow ++) {
                        for (int icol = 0; icol < ap.get_col_size()/2; icol++){
                            int index_i = irow* ap.get_col_size()/2 + icol;
                            int index_j = (irow*2+row_set[is]) * ap.get_col_size() + icol*2+col_set[is];
                            matrix_out[index_i] = matrix_in[index_j];
                        }
                    }
                }
            }
#ifdef __MPI
            hamilt::transferParallels2Serials( *(this->dm2d_tmp), this->dmr_gint[is]);
#else
            this->dmr_gint[is]->set_zero();
            this->dmr_gint[is]->add(*(this->dm2d_tmp));
#endif
        }//is=4
        delete this->dm2d_tmp;
        this->dm2d_tmp = nullptr;
    }
    ModuleBase::timer::tick("Gint", "transfer_DMR");
}
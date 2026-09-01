#include "dftu_lcao.h"

#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"

#include <complex>
#include <vector>

Plus_U::Plus_U()
{}

Plus_U::~Plus_U()
{}

void Plus_U::init(UnitCell& cell,
                const Parallel_Orbitals* pv,
                const int npol,
                const int nspin,
                const std::vector<int>& orbital_corr,
                const bool yukawa_potential,
                const double yukawa_lambda,
                const std::string& global_readin_dir,
                const std::string& global_out_dir,
                const std::string& init_chg,
                const int nlocal,
                const std::string& ks_solver,
                const std::string& device,
                const int kpar,
                const std::vector<double>& hubbard_u,
                const double uramping,
                const int occ_mat_ctrl,
                const int mixing_dftu
#ifdef __LCAO
                , const LCAO_Orbitals* orb
#endif
                )
{
    ModuleBase::TITLE("Plus_U", "init");

    this->yukawa_lambda = yukawa_lambda;

#ifdef __LCAO
    ptr_orb_ = orb;
    if(ptr_orb_ != nullptr)
    {
        orb_cutoff_ = orb->cutoffs();
    }
#endif

    if (pv != nullptr)
    {
        const int global_rows = pv->get_global_row_size();
        const int global_cols = pv->get_global_col_size();
        if (global_rows != global_cols)
        {
            ModuleBase::WARNING_QUIT("Plus_U::init", "Global row and column dimensions do not match");
        }
        if (nlocal != global_rows)
        {
            ModuleBase::WARNING_QUIT("Plus_U::init", "nlocal does not match global matrix dimension");
        }
    }

    this->init_base(cell,
                    npol,
                    nspin,
                    orbital_corr,
                    yukawa_potential,
                    global_readin_dir,
                    global_out_dir,
                    init_chg,
                    device,
                    kpar,
                    hubbard_u,
                    uramping,
                    occ_mat_ctrl,
                    mixing_dftu);
    return;
}

// uramping_update() and u_converged() are now implemented in
// dftu_base.cpp as Plus_U_Base methods (inherited by Plus_U).

#ifdef __LCAO

void Plus_U::set_dmr(const elecstate::DensityMatrix<std::complex<double>, double>* dmr)
{
    this->dm_in_dftu_cd = dmr;
    return;
}

void Plus_U::set_dmr(const elecstate::DensityMatrix<double, double>* dmr)
{
    this->dm_in_dftu_d = dmr;
    return;
}

const hamilt::HContainer<double>* Plus_U::get_dmr(int ispin) const
{
    if (this->dm_in_dftu_d != nullptr)
    {
        return this->dm_in_dftu_d->get_DMR_pointer(ispin + 1);
    }
    else if (this->dm_in_dftu_cd != nullptr)
    {
        return this->dm_in_dftu_cd->get_DMR_pointer(ispin + 1);
    }
    else
    {
        return nullptr;
    }
}

#endif

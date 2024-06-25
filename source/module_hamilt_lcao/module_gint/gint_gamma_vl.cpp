//=========================================================
// REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "module_base/blas_connector.h"
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

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int* npcol, int* myprow, int* mypcol);
    void Cblacs_pinfo(int* myid, int* nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int* prow, int* pcol);
}

void Gint_Gamma::cal_vlocal(Gint_inout* inout, bool new_e_iteration)
{
    const int max_size = this->gridt->max_atom;
    const int lgd = this->gridt->lgd;

    if (inout->job == Gint_Tools::job_type::vlocal || inout->job == Gint_Tools::job_type::vlocal_meta)
    {
        if (max_size > 0 && lgd > 0)
        {
            this->hRGint->set_zero();
        }

        this->cal_gint(inout);
    }
}

#ifdef __MPI
#include "module_hamilt_lcao/module_hcontainer/hcontainer_funcs.h"
#endif
void Gint_Gamma::transfer_pvpR(hamilt::HContainer<double>* hR, const UnitCell* ucell)
{
    ModuleBase::TITLE("Gint_Gamma", "transfer_pvpR");
    ModuleBase::timer::tick("Gint_Gamma", "transfer_pvpR");

    for (int iap = 0; iap < this->hRGint->size_atom_pairs(); iap++)
    {
        auto& ap = this->hRGint->get_atom_pair(iap);
        const int iat1 = ap.get_atom_i();
        const int iat2 = ap.get_atom_j();
        if (iat1 > iat2)
        {
            // fill lower triangle matrix with upper triangle matrix
            // gamma_only case, only 1 R_index in each AtomPair
            // the upper <IJR> is <iat2, iat1, 0>
            const hamilt::AtomPair<double>* upper_ap = this->hRGint->find_pair(iat2, iat1);
#ifdef __DEBUG
            assert(upper_ap != nullptr);
#endif
            double* lower_matrix = ap.get_pointer(0);
            for (int irow = 0; irow < ap.get_row_size(); ++irow)
            {
                for (int icol = 0; icol < ap.get_col_size(); ++icol)
                {
                    *lower_matrix++ = upper_ap->get_value(icol, irow);
                }
            }
        }
    }

#ifdef __MPI
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1)
    {
        hR->add(*this->hRGint);
    }
    else
    {
        hamilt::transferSerials2Parallels(*this->hRGint, hR);
    }
#else
    hR->add(*this->hRGint);
#endif

    ModuleBase::timer::tick("Gint_Gamma", "transfer_pvpR");

    return;
}

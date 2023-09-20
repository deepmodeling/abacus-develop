//Obsolete code, not compiled
//Please fix it by removing globalc::hm

#include "diago_lapack.h"

#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/global_variable.h"
#include "module_base/lapack_connector.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"

typedef hamilt::MatrixBlock<double> matd;
typedef hamilt::MatrixBlock<std::complex<double>> matcd;

namespace hsolver
{

void DiagoLapack::diag(hamilt::Hamilt<std::complex<double>> *phm_in, psi::Psi<std::complex<double>> &psi, double *eigenvalue_in)
{
    ModuleBase::TITLE("DiagoLapack", "diag");
    assert(GlobalV::NPROC == 1);

    ModuleBase::ComplexMatrix Htmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    ModuleBase::ComplexMatrix Stmp(GlobalV::NLOCAL, GlobalV::NLOCAL);
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
    for (int i = 0; i < GlobalV::NLOCAL; i++)
    {
        for (int j = 0; j < GlobalV::NLOCAL; j++)
        {
            Htmp(i, j) = h_mat.p[i * GlobalV::NLOCAL + j];
            Stmp(i, j) = s_mat.p[i * GlobalV::NLOCAL + j];
        }
    }

    //----------------------------
    // keep this for tests
    //    out.printcm_norm("Lapack_H", Htmp, 1.0e-5);
    //    out.printcm_norm("Lapack_S", Stmp, 1.0e-5);
    //----------------------------

    double *en = new double[GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(en, GlobalV::NLOCAL);

    ModuleBase::ComplexMatrix hvec(GlobalV::NLOCAL, GlobalV::NBANDS);
    GlobalC::hm.diagH_LAPACK(GlobalV::NLOCAL, GlobalV::NBANDS, Htmp, Stmp, GlobalV::NLOCAL, en, hvec);

    if (GlobalV::NSPIN != 4)
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
            {
                psi.get_pointer()[ib * GlobalV::NLOCAL + iw] = hvec(iw, ib);
            }
        }
    }
    /*
    else
    {
        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {
            for (int iw = 0; iw < GlobalV::NLOCAL / GlobalV::NPOL; iw++)
            {
                wfc_k_grid[ib][iw] = hvec(iw * GlobalV::NPOL, ib);
                wfc_k_grid[ib][iw + GlobalV::NLOCAL / GlobalV::NPOL] = hvec(iw * GlobalV::NPOL + 1, ib);
            }
        }
    }
    */

    // energy for k-point ik
    for (int ib = 0; ib < GlobalV::NBANDS; ib++)
    {
        eigenvalue_in[ib] = en[ib];
    }
}

void DiagoLapack::diag(hamilt::Hamilt<std::complex<double>> *phm_in, psi::Psi<double> &psi, double *eigenvalue_in)
{
    ModuleBase::TITLE("DiagoLapack", "diag");
}

} // namespace hsolver
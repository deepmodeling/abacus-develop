#include "boundary_fix.h"
#include "source_base/libm/libm.h"
#include "source_base/constants.h"
#include "source_base/vector3.h"

namespace module_rt{

void reset_matrix_boundary(const UnitCell& ucell,
                           const K_Vectors& kv,
                           const Parallel_Orbitals* pv,
                           ct::Tensor& hk_last,
                           ct::Tensor& sk_last,
                           psi::Psi<std::complex<double>>* psi_last,
                           const int len_hs)
{
    const ModuleBase::Vector3<int> zero = {0, 0, 0};
    for(int iat = 0; iat < ucell.nat; iat++)
    {
        const int it = ucell.iat2it[iat];
        const int ia = ucell.iat2ia[iat];
        if(ucell.atoms[it].boundary_shift[ia]!=zero)
        {
            const auto& rshift = ucell.atoms[it].boundary_shift[ia];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for(int ik = 0; ik < kv.get_nks(); ik++)
            {
                const ModuleBase::Vector3<double> tmp_rshift(rshift.x, rshift.y, rshift.z);
                const double arg = -kv.kvec_d[ik] * tmp_rshift * ModuleBase::TWO_PI;
                //skip unrelevent ik
                if(arg==0)continue;
                //calculate correction phase
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                const std::complex<double> phase = std::complex<double>(cosp, sinp);
                //phase correction for Hamiltionian, overlap matrix and c vec.
                module_rt::boundary_shift_mat(phase, hk_last.template data<std::complex<double>>() + ik * len_hs, pv, iat);
                module_rt::boundary_shift_mat(phase, sk_last.template data<std::complex<double>>() + ik * len_hs, pv, iat);
                psi_last->fix_k(ik);
                module_rt::boundary_shift_c(phase, psi_last[0].get_pointer(), pv, iat);
            }
        }
    }
    return;
}

void boundary_shift_mat(const std::complex<double>& phase,
                        std::complex<double>* matk,
                        const Parallel_Orbitals* pv,
                        const int iat)
{
    const std::complex<double> phase_conj = std::conj(phase);
    int row0 = pv->atom_begin_row[iat];
    int col0 = pv->atom_begin_col[iat];
    std::complex<double>* p_matkc = matk + col0 * pv->get_row_size();
    for(int nu = 0; nu < pv->get_col_size(iat); ++nu)
    {
        
        BlasConnector::scal(pv->get_row_size(),
                            phase,
                            p_matkc,
                            1);
        p_matkc += pv->get_row_size();
    }
    std::complex<double>* p_matkr = matk + row0;
    for(int mu = 0; mu < pv->get_row_size(iat); ++mu)
    {
        BlasConnector::scal(pv->get_col_size(),
                            phase_conj,
                            p_matkr,
                            pv->get_row_size());
        p_matkr += 1;
    }
    return;
}

void boundary_shift_c(const std::complex<double>& phase,
                      std::complex<double>* psi_k_last,
                      const Parallel_Orbitals* pv,
                      const int iat)
{
    const std::complex<double> phase_conj = std::conj(phase);
    int row0 = pv->atom_begin_row[iat];
    std::complex<double>* p_ck = psi_k_last + row0;
    for(int nu = 0; nu < pv->get_row_size(iat); ++nu)
    {
        BlasConnector::scal(pv->ncol_bands,
                            phase_conj,
                            p_ck,
                            pv->get_row_size());
        p_ck+=1;
    }
    return;
}
} //namespace module_rt
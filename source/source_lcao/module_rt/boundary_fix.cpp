#include "boundary_fix.h"
#include "source_base/libm/libm.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/vector3.h"

namespace module_rt{

void reset_matrix_boundary(const UnitCell& ucell,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const int nks,
                           const ModuleBase::Vector3<double>& cart_At_last,
                           const Parallel_Orbitals* pv,
                           ct::Tensor& hk_last,
                           ct::Tensor& sk_last,
                           psi::Psi<std::complex<double>>* psi_last,
                           const size_t len_hs)
{
    ModuleBase::TITLE("module_rt", "reset_matrix_boundary");
    ModuleBase::timer::start("module_rt", "reset_matrix_boundary");
    // Tensor-LAPACK keeps full H/S matrices on rank zero, while psi stays distributed.
    const bool gathered_matrices = len_hs != pv->para_2d().size();
    Parallel_Orbitals full_distribution;
    const Parallel_Orbitals* matrix_distribution = pv;
    bool update_matrices = true;
    if (gathered_matrices)
    {
        const int basis_size = pv->para_2d().nrglb();
        full_distribution.para_2d().set_serial(basis_size, basis_size);
        full_distribution.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, basis_size);
        matrix_distribution = &full_distribution;
        const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
        update_matrices = domain.rank() == 0;
    }
    const ModuleBase::Vector3<int> zero = {0, 0, 0};
    for(size_t iat = 0; iat < ucell.nat; iat++)
    {
        const size_t it = ucell.iat2it[iat];
        const size_t ia = ucell.iat2ia[iat];
        if(ucell.atoms[it].boundary_shift[ia]!=zero)
        {
            const auto& rshift = ucell.atoms[it].boundary_shift[ia];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for(int ik = 0; ik < nks; ik++)
            {
                const ModuleBase::Vector3<double> tmp_rshift(rshift.x, rshift.y, rshift.z);
                // Hybrid folding uses exp(i * (2*pi*k.R + A.R_cart)).
                // Transform the saved state with its own field endpoint.
                const ModuleBase::Vector3<double> cart_shift = tmp_rshift * ucell.latvec * ucell.lat0;
                const double arg = -kvec_d[ik] * tmp_rshift * ModuleBase::TWO_PI - cart_At_last * cart_shift;
                //skip unrelevent ik
                if(arg==0)continue;
                //calculate correction phase
                double sinp = 0.0, cosp = 0.0;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                const std::complex<double> phase = std::complex<double>(cosp, sinp);
                //phase correction for Hamiltionian, overlap matrix and c vec.
                if (update_matrices)
                {
                    std::complex<double>* h = hk_last.template data<std::complex<double>>() + ik * len_hs;
                    std::complex<double>* s = sk_last.template data<std::complex<double>>() + ik * len_hs;
                    module_rt::boundary_shift_mat(phase, h, matrix_distribution, iat);
                    module_rt::boundary_shift_mat(phase, s, matrix_distribution, iat);
                }
                if (pv->para_2d().nrloc() > 0 && pv->ncol_bands > 0)
                {
                    // Saved LCAO coefficients use k-first storage. Index it without
                    // changing the shared Psi cursor in this OpenMP loop.
                    std::complex<double>* psi_k = &(*psi_last)(ik, 0, 0);
                    module_rt::boundary_shift_c(phase, psi_k, pv, iat);
                }
            }
        }
    }
    ModuleBase::timer::end("module_rt", "reset_matrix_boundary");
    return;
}

void boundary_shift_mat(const std::complex<double>& phase,
                        std::complex<double>* matk,
                        const Parallel_Orbitals* pv,
                        const size_t iat)
{
    const std::complex<double> phase_conj = std::conj(phase);
    size_t row0 = pv->atom_begin_row[iat];
    size_t col0 = pv->atom_begin_col[iat];
    std::complex<double>* p_matkc = matk + col0 * pv->para_2d().nrloc();
    for(size_t nu = 0; nu < pv->get_ncol_atom(iat); ++nu)
    {
        
        BlasConnector::scal(pv->para_2d().nrloc(),
                            phase,
                            p_matkc,
                            1);
        p_matkc += pv->para_2d().nrloc();
    }
    std::complex<double>* p_matkr = matk + row0;
    for(size_t mu = 0; mu < pv->get_nrow_atom(iat); ++mu)
    {
        BlasConnector::scal(pv->para_2d().ncloc(),
                            phase_conj,
                            p_matkr,
                            pv->para_2d().nrloc());
        p_matkr += 1;
    }
    return;
}

void boundary_shift_c(const std::complex<double>& phase,
                      std::complex<double>* psi_k_last,
                      const Parallel_Orbitals* pv,
                      const size_t iat)
{
    const std::complex<double> phase_conj = std::conj(phase);
    size_t row0 = pv->atom_begin_row[iat];
    std::complex<double>* p_ck = psi_k_last + row0;
    for(size_t nu = 0; nu < pv->get_nrow_atom(iat); ++nu)
    {
        BlasConnector::scal(pv->ncol_bands,
                            phase_conj,
                            p_ck,
                            pv->para_2d().nrloc());
        p_ck+=1;
    }
    return;
}
} //namespace module_rt

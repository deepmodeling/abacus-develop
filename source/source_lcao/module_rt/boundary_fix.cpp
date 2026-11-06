#include "boundary_fix.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/block_cyclic_distribution.h"
#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/parallel_cell.h"
#include "source_base/vector3.h"

#include <cassert>

namespace module_rt{

namespace
{
struct LocalBasisRange
{
    int ibas0;
    int nbasis;
};

enum class BasisAxis
{
    Row,
    Column
};

LocalBasisRange local_basis_range(const BasisFlattenedIndexing& basis,
                                  const BlockCyclicDistribution& layout,
                                  const int iat,
                                  const BasisAxis axis)
{
    const BasisRange global_range = basis.range(iat);
    int ibas0 = -1;
    int nbasis = 0;
    for (std::size_t ibas = global_range.ibas0; ibas < global_range.ibas1; ++ibas)
    {
        const int iloc = axis == BasisAxis::Row
                             ? layout.irloc(static_cast<int>(ibas))
                             : layout.icloc(static_cast<int>(ibas));
        if (iloc >= 0)
        {
            if (ibas0 < 0)
            {
                ibas0 = iloc;
            }
            ++nbasis;
        }
    }
    const LocalBasisRange local_range = {ibas0, nbasis};
    return local_range;
}

void boundary_shift_mat(const std::complex<double>& phase,
                        std::complex<double>* matk,
                        const LocalBasisRange& rows,
                        const LocalBasisRange& columns,
                        const int nrloc,
                        const int ncloc)
{
    const std::complex<double> phase_conj = std::conj(phase);
    if (columns.nbasis > 0)
    {
        std::complex<double>* p_matkc = matk + columns.ibas0 * nrloc;
        for(int nu = 0; nu < columns.nbasis; ++nu)
        {
            BlasConnector::scal(nrloc, phase, p_matkc, 1);
            p_matkc += nrloc;
        }
    }
    if (rows.nbasis > 0)
    {
        std::complex<double>* p_matkr = matk + rows.ibas0;
        for(int mu = 0; mu < rows.nbasis; ++mu)
        {
            BlasConnector::scal(ncloc, phase_conj, p_matkr, nrloc);
            p_matkr += 1;
        }
    }
}

void boundary_shift_c(const std::complex<double>& phase,
                      std::complex<double>* psi_k,
                      const int nbands,
                      const int nbasis,
                      const LocalBasisRange& rows)
{
    const std::complex<double> phase_conj = std::conj(phase);
    std::complex<double>* p_ck = psi_k + rows.ibas0;
    for(int nu = 0; nu < rows.nbasis; ++nu)
    {
        BlasConnector::scal(nbands, phase_conj, p_ck, nbasis);
        p_ck += 1;
    }
}
} // namespace

void reset_matrix_boundary(const UnitCell& ucell,
                           const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                           const ModuleBase::Vector3<double>& cart_At_last,
                           const BasisFlattenedIndexing& basis,
                           const BlockCyclicDistribution& layout,
                           ct::Tensor& hk_last,
                           ct::Tensor& sk_last,
                           psi::Psi<std::complex<double>>* psi_last)
{
    ModuleBase::TITLE("module_rt", "reset_matrix_boundary");
    ModuleBase::timer::start("module_rt", "reset_matrix_boundary");
    const int nks = static_cast<int>(kvec_d.size());
    assert(nks > 0);
    if (psi_last == nullptr)
    {
        ModuleBase::timer::end("module_rt", "reset_matrix_boundary");
        return;
    }
    assert(hk_last.NumElements() % nks == 0);
    const std::size_t len_hs =
        static_cast<std::size_t>(hk_last.NumElements()) / static_cast<std::size_t>(nks);
    // Tensor-LAPACK keeps full H/S matrices on rank zero, while psi stays distributed.
    const bool gathered_matrices = len_hs != layout.size();
    BlockCyclicDistribution full_layout;
    const BlockCyclicDistribution* matrix_layout = &layout;
    bool update_matrices = true;
    if (gathered_matrices)
    {
        const int basis_size = layout.nrglb();
        full_layout.set_serial(basis_size, basis_size);
        matrix_layout = &full_layout;
        const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
        update_matrices = domain.rank() == 0;
    }
    const int matrix_nrloc = matrix_layout->nrloc();
    const int matrix_ncloc = matrix_layout->ncloc();
    const int psi_nbands = psi_last->get_nbands();
    const int psi_nbasis = psi_last->get_nbasis();
    const ModuleBase::Vector3<int> zero = {0, 0, 0};
    for(size_t iat = 0; iat < ucell.nat; iat++)
    {
        const size_t it = ucell.iat2it[iat];
        const size_t ia = ucell.iat2ia[iat];
        if(ucell.atoms[it].boundary_shift[ia]!=zero)
        {
            const auto& rshift = ucell.atoms[it].boundary_shift[ia];
            const LocalBasisRange matrix_rows =
                local_basis_range(basis, *matrix_layout, static_cast<int>(iat), BasisAxis::Row);
            const LocalBasisRange matrix_columns =
                local_basis_range(basis, *matrix_layout, static_cast<int>(iat), BasisAxis::Column);
            const LocalBasisRange psi_rows =
                local_basis_range(basis, layout, static_cast<int>(iat), BasisAxis::Row);
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
                    boundary_shift_mat(phase,
                                       h,
                                       matrix_rows,
                                       matrix_columns,
                                       matrix_nrloc,
                                       matrix_ncloc);
                    boundary_shift_mat(phase,
                                       s,
                                       matrix_rows,
                                       matrix_columns,
                                       matrix_nrloc,
                                       matrix_ncloc);
                }
                if (psi_rows.nbasis > 0 && psi_nbands > 0)
                {
                    // Saved LCAO coefficients use k-first storage. Index it without
                    // changing the shared Psi cursor in this OpenMP loop.
                    std::complex<double>* psi_k = &(*psi_last)(ik, 0, 0);
                    boundary_shift_c(phase, psi_k, psi_nbands, psi_nbasis, psi_rows);
                }
            }
        }
    }
    ModuleBase::timer::end("module_rt", "reset_matrix_boundary");
    return;
}

} //namespace module_rt

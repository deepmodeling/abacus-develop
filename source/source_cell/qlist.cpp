// ============================================================
// QList: q-point mesh generation and star reduction.
// ============================================================

#include "qlist.h"

#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/tool_quit.h"

namespace ModuleCell {

QList::QList() {}

QList::~QList() {}

void QList::generate_mesh(UnitCell& ucell, ModuleSymmetry::Symmetry& symm,
                          const std::vector<int>& mp_grid, bool use_irreps) {
    (void)use_irreps;

    if (mp_grid.size() != 3)
    {
        ModuleBase::WARNING_QUIT("QList::generate_mesh", "mp_grid must have three components.");
    }

    this->is_mp = true;
    this->nmp[0] = mp_grid[0];
    this->nmp[1] = mp_grid[1];
    this->nmp[2] = mp_grid[2];

    // Gamma-centered Monkhorst-Pack q mesh (k_type = 0), zero offset.
    const double offset[3] = {0.0, 0.0, 0.0};
    this->Monkhorst_Pack(this->nmp, offset, 0);

    this->nkstot_full = this->nkstot;
    this->nks = this->nkstot;

    // Star reduction: always use symmetry, always include the -q partner.
    bool match = true;
    std::string skpt;
    this->reduce_by_symmetry(ucell, symm, true, skpt, match);
    if (!match)
    {
        ModuleBase::WARNING("QList::generate_mesh",
                            "Reciprocal lattice is incompatible with the real-space lattice. "
                            "Falling back to the unreduced q-point mesh.");
        this->nkstot = this->nks = this->nkstot_full;
    }

    // weights sum to 1 (average over the full Brillouin zone)
    this->normalize_wk(1);

    // little-group irreducible-representation data
    this->get_irreps(ucell, symm);
}

void QList::read_from_file(const std::string& filename, UnitCell& ucell) {
    (void)filename;
    (void)ucell;
}

std::vector<int> QList::get_irrep_modes(int q_idx, int irrep_idx) const {
    if (q_idx < 0 || q_idx >= this->nkstot || irrep_idx < 0 || irrep_idx >= (int)this->nirr_[q_idx])
    {
        return std::vector<int>();
    }
    return this->irrep_modes_[q_idx][irrep_idx];
}

void QList::reduce_by_symmetry(const UnitCell& ucell,
                               const ModuleSymmetry::Symmetry& symm,
                               bool use_symm,
                               std::string& skpt,
                               bool& match) {
    (void)skpt;
    // q-points are spin-free: build the point-group operations and always
    // double them by the time-reversal operation -q (no magnetic group).
    std::vector<ModuleBase::Matrix3> kgmatrix(48 * 2);
    ModuleBase::Matrix3 inv(-1, 0, 0, 0, -1, 0, 0, 0, -1);

    ModuleBase::Matrix3 q_vec; // k-lattice basis of the q mesh
    int nrotkm = 0;
    if (!this->build_star_ops(ucell, symm, use_symm, q_vec, kgmatrix, nrotkm))
    {
        match = false;
        return;
    }
    if (nrotkm == 0)
    {
        // no operations to apply: the mesh stays unreduced
        match = true;
        return;
    }

    bool include_inv = false;
    for (int i = 0; i < nrotkm; ++i)
    {
        if (kgmatrix[i] == inv)
        {
            include_inv = true;
            break;
        }
    }
    if (!include_inv)
    {
        for (int i = 0; i < nrotkm; ++i)
        {
            kgmatrix[i + nrotkm] = inv * kgmatrix[i];
        }
        nrotkm *= 2;
    }

    ModuleBase::Matrix3* kkmatrix = new ModuleBase::Matrix3[nrotkm];
    symm.gmatrix_convert(kgmatrix.data(), kkmatrix, nrotkm, ucell.G, q_vec);

    std::vector<ModuleBase::Vector3<double>> qvec_ibz;
    std::vector<double> wk_ibz;
    std::vector<int> ibz_index;
    std::vector<int> ibz2bz;
    this->reduce_ibz(kgmatrix.data(), nrotkm, ucell.G, q_vec, kkmatrix, symm.epsilon, qvec_ibz, wk_ibz, ibz_index, ibz2bz);

    delete[] kkmatrix;

    // update the reduced q-point list (no spin expansion)
    const int nq_ibz = qvec_ibz.size();
    this->nkstot = this->nks = nq_ibz;
    this->kvec_d.resize(this->nkstot);
    this->wk.resize(this->nkstot);
    for (int i = 0; i < this->nkstot; ++i)
    {
        this->kvec_d[i] = qvec_ibz[i];
        this->wk[i] = wk_ibz[i];
    }
    this->kd_done = true;
    this->kc_done = false;

    match = true;
    return;
}

void QList::get_irreps(const UnitCell& ucell, const ModuleSymmetry::Symmetry& symm) {
    (void)ucell;

    // Decompose each q-point via its little group (placeholder: one
    // fully-symmetric A1 irrep per q-point; the LittleGroup projection-operator
    // basis is filled in a later iteration).
    this->nirr_.assign(this->nkstot, 0);
    this->irrep_modes_.assign(this->nkstot, std::vector<std::vector<int>>());
    for (int iq = 0; iq < this->nkstot; ++iq)
    {
        this->little_group_.set_q(this->kvec_d[iq], symm);
        this->nirr_[iq] = this->little_group_.get_nirr();
        this->irrep_modes_[iq].resize(this->nirr_[iq]);
        for (int iirr = 0; iirr < this->nirr_[iq]; ++iirr)
        {
            this->irrep_modes_[iq][iirr] = this->little_group_.get_mode_basis(iirr);
        }
    }
}

} // namespace ModuleCell
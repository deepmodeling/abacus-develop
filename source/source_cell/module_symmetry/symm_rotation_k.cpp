#include "symm_rotation_k.h"
#include "source_base/constants.h"
#include <cmath>
#include "source_base/parallel_reduce.h"
#include "source_base/parallel_global.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/tool_title.h"
#include "source_base/timer.h"

namespace ModuleSymmetry
{
    std::vector<TC> Symmetry_rotation_k::get_bvk_cells(const TC& period)
    {
        std::vector<TC> cells;
        cells.reserve(static_cast<size_t>(period[0]) * period[1] * period[2]);
        for (int ix = 0; ix < period[0]; ++ix) {
            for (int iy = 0; iy < period[1]; ++iy) {
                for (int iz = 0; iz < period[2]; ++iz) {
                    cells.push_back({ix, iy, iz});
        } } }
        return cells;
    }

    void Symmetry_rotation_k::cal_Ms(const K_Vectors& kv,
        const UnitCell& ucell, const Parallel_2D& pv, const int nspin)
    {
        ModuleBase::TITLE("Symmetry_rotation_k", "cal_Ms");
        ModuleBase::timer::start("Symmetry_rotation_k", "cal_Ms");

        this->nspin_ = nspin;
        this->nsym_ = ucell.symm.nrotk;
        this->nanti_ = ucell.symm.nrotk_anti;
        this->magnetic_nspin4_ = ucell.symm.magnetic_nspin4;
        this->eps_ = ucell.symm.epsilon;
        if (this->irs_.invmap_.empty())
        {
            this->irs_.invmap_.resize(ucell.symm.nrotk);
            ucell.symm.gmatrix_invmap(ucell.symm.gmatrix, ucell.symm.nrotk, this->irs_.invmap_.data());
        }
        // 1. calculate the rotation matrix in real spherical harmonics representation for each symmetry operation: [T_l (isym)]_mm'
        const int nop_tot = this->nsym_ + this->nanti_;
        std::vector<ModuleBase::Matrix3> gmatc(nop_tot);
        for (int i = 0;i < nsym_;++i) { gmatc[i] = this->irs_.direct_to_cartesian(ucell.symm.gmatrix[i], ucell.latvec); }
        for (int j = 0;j < this->nanti_;++j)
        { gmatc[nsym_ + j] = this->irs_.direct_to_cartesian(ucell.symm.gmatrix_anti[j], ucell.latvec); }
        this->cal_rotmat_Slm(gmatc.data(), std::max(this->abfs_Lmax_, ucell.lmax), nop_tot);

        // 1.5 (nspin=4) the SU(2) spin-1/2 rotation U(isym) for each symmetry operation. The AO
        // rotation matrix M becomes the spinor operator T(isym) (x) U(isym) so that the same
        // gemm D(k)=M^dagger D(k_ibz) M rotates both the orbital and the spin part at once.
        // For an antiunitary element Theta*g only the spatial part g enters M here; the Theta
        // (sigma_y (.)^* sigma_y) is applied afterwards in restore_dm.
        std::vector<SpinRotation::Su2> spin_U(nop_tot, SpinRotation::Su2{ 1.0, 0.0, 0.0, 1.0 });
        if (this->nspin_ == 4)
        {
            for (int i = 0;i < nop_tot;++i) { spin_U[i] = SpinRotation::so3_to_su2(gmatc[i]); }
        }
        this->spin_U_ = spin_U;  // keep for restore_HR_nspin4 (real-space EXX H(R) spin mixing)

        // 2. calculate the rotation matrix in AO-representation for each ibz_kpoint and symmetry operation: M(k, isym)
        int nks_ibz = kv.kstars.size(); // kv.nks = 2 * kv.nks_ibz when nspin=2
        this->Ms_.assign(nks_ibz, {});
        this->little_groups_.assign(nks_ibz, {});

        // (k-point pools, KPAR>1) kv.kvec_d only holds the k-points owned by this pool, so
        // kv.kvec_d[ik_ibz] is only valid for ik_ibz < kv.para_k.nks_np and is otherwise either
        // out of range or (after a caller's spin-doubling resize) a meaningless zero placeholder.
        // Gather the (small, size nks_ibz) global ibz-representative k-vector list once so every
        // pool builds the correct rotation matrix for every ibz k, not just the ones it owns.
#ifdef __MPI
        // inlined equivalent of Parallel_Kpoints::gatherkvec (avoided as a direct call so this
        // class does not pull in a link dependency on parallel_kpoints.cpp for every target that
        // links the "symmetry" library): every rank in the owning pool holds the same local
        // k-vectors, so only the pool root contributes to the MPI_Allreduce, matching gatherkvec.
        int world_rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        const bool is_pool_root = (world_rank == kv.para_k.get_startpro_pool(kv.para_k.my_pool));
        std::vector<ModuleBase::Vector3<double>> kvec_d_ibz_global(nks_ibz, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
        for (int i = 0; i < kv.para_k.nks_np; ++i)
        {
            if (is_pool_root) { kvec_d_ibz_global[i + kv.para_k.startk_pool[kv.para_k.my_pool]] = kv.kvec_d[i]; }
        }
        MPI_Allreduce(MPI_IN_PLACE, kvec_d_ibz_global.data(), 3 * nks_ibz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz_global = kv.kvec_d;
#endif

        // A k-star contains only one operation per distinct k point. The other
        // operations fixing k (modulo a reciprocal lattice vector) must still
        // be averaged: a finite-grid SCF density need not respect this little group.
        for (int ik_ibz = 0; ik_ibz < nks_ibz; ++ik_ibz)
        {
            std::set<int> needed;
            for (const auto& member : kv.kstars[ik_ibz])
            {
                const int op = (!this->magnetic_nspin4_ && member.first >= nsym_)
                                   ? member.first - nsym_ : member.first;
                needed.insert(op);
            }
            for (int op = 0; op < nsym_; ++op)
            {
                const auto delta = kvec_d_ibz_global[ik_ibz] * ucell.symm.kgmatrix[op] - kvec_d_ibz_global[ik_ibz];
                if (std::abs(delta.x - std::round(delta.x)) < this->eps_
                    && std::abs(delta.y - std::round(delta.y)) < this->eps_
                    && std::abs(delta.z - std::round(delta.z)) < this->eps_)
                {
                    this->little_groups_[ik_ibz].push_back(op);
                    needed.insert(op);
                }
            }
            for (const int op : needed)
            {
                this->Ms_[ik_ibz][op] = this->contruct_2d_rot_mat_ao(
                    ucell.symm, ucell.atoms, ucell.st, kvec_d_ibz_global[ik_ibz], op, pv, spin_U[op]);
            }
        }

        ModuleBase::timer::end("Symmetry_rotation_k", "cal_Ms");
    }

    std::vector<std::vector<std::complex<double>>> Symmetry_rotation_k::restore_dm(const K_Vectors& kv,
        const std::vector<std::vector<std::complex<double>>>& dm_k_ibz, const Parallel_2D& pv)const
    {
        ModuleBase::TITLE("Symmetry_rotation_k", "restore_dm");
        ModuleBase::timer::start("Symmetry_rotation_k", "restore_dm");
        std::vector<std::vector<std::complex<double>>> dm_k_full;
        int nspin0 = this->nspin_ == 2 ? 2 : 1;
        dm_k_full.reserve(kv.get_nkstot_nospin() * nspin0); //nkstot_nospin didn't doubled by spin
        int nk = kv.get_nkstot() / nspin0;

        // (nspin=4) Sigma_y = I (x) sigma_y for the time-reversal spin flip; k-independent, build once.
        std::vector<std::complex<double>> sigma_y;
        if (this->nspin_ == 4) { sigma_y = this->set_sigma_y_2d(pv); }

        for (int is = 0;is < nspin0;++is)
        {
            for (int ik_ibz = 0;ik_ibz < nk;++ik_ibz)
            {
                // P_k D = |G_k|^{-1} sum_g M_g^T D M_g^*. This preserves
                // Hermiticity and makes restoration independent of the chosen
                // star representative; rotating just one arbitrary D does not.
                const auto& little_group = this->little_groups_.at(ik_ibz);
                assert(!little_group.empty());
                std::vector<std::complex<double>> projected = dm_k_ibz[ik_ibz + is * nk];
                if (little_group.size() > 1)
                {
                    std::fill(projected.begin(), projected.end(), 0.0);
                    for (const int op : little_group)
                    {
                        const auto rotated = this->rot_matrix_ao(
                            dm_k_ibz[ik_ibz + is * nk], ik_ibz, little_group.size(), op, pv);
                        for (size_t i = 0; i < projected.size(); ++i)
                        {
                            projected[i] += rotated[i];
                        }
                    }
                }
                for (auto& isym_kvd : kv.kstars[ik_ibz])
                {
                    if (isym_kvd.first == 0)
                    {
                        double factor = 1.0 / static_cast<double>(kv.kstars[ik_ibz].size());
                        std::vector<std::complex<double>> dm_scaled(pv.get_local_size());
                        for (int i = 0;i < pv.get_local_size();++i) { dm_scaled[i] = factor * projected[i]; }
                        dm_k_full.push_back(dm_scaled);
                    }
                    else if (isym_kvd.first < nsym_)
                    { //space group operations
                        dm_k_full.push_back(this->rot_matrix_ao(projected, ik_ibz, kv.kstars[ik_ibz].size(), isym_kvd.first, pv));
                    }
                    else
                    {    // antiunitary elements: Theta * (spatial operation)
                        // D(Theta*g k_ibz) = sigma_y [D(g k_ibz)]^* sigma_y with D(g k_ibz) = M^dagger D M.
                        // For nspin=4, first do the (non-conjugated) spatial rotation, then the spin flip;
                        // for nspin<4 (Theta=K) the original TRS_conj path already gives the conjugate.
                        //
                        // Which spatial operation the index denotes depends on the regime, matching
                        // how the k-reduction filled kgmatrix[] (see K_Vectors::reduce_by_symmetry):
                        //  - nspin=4 magnetic (Shubnikov): index j+nsym_ is the antiunitary element
                        //    Theta*gmatrix_anti[j]; its Ms is stored under the RAW key j+nsym_.
                        //  - otherwise (grey group / nspin<4): index i+nsym_ is Theta*gmatrix[i],
                        //    i.e. the unitary operation i, whose Ms is stored under key i.
                        const int isym_M = this->magnetic_nspin4_ ? isym_kvd.first : (isym_kvd.first - nsym_);
                        if (this->nspin_ == 4)
                        {
                            // m=0: gray group: the space-group part of anti-unitary elements are the same of the unitary elements, isym_M < nsym_
                            // m!=0: Shubnikov group: using different space-group part of anti-unitary elements stored in gmatrix_anti with isym_M >= nsym_
                            dm_k_full.push_back(this->trs_spin_rotate(
                                this->rot_matrix_ao(projected, ik_ibz, kv.kstars[ik_ibz].size(), isym_M, pv, false),
                                sigma_y, pv, 1.0));
                        }
                        else
                        {
                            dm_k_full.push_back(this->rot_matrix_ao(projected, ik_ibz, kv.kstars[ik_ibz].size(), isym_M, pv, true));
                        }
                    }
                }
            }
        }
        ModuleBase::timer::end("Symmetry_rotation_k", "restore_dm");
        return dm_k_full;
    }
    std::vector<std::vector<double>> Symmetry_rotation_k::restore_dm(const K_Vectors& kv,
        const std::vector<std::vector<double>>& dm_k_ibz, const Parallel_2D& pv)const
    {
        return dm_k_ibz;// do nothing for gamma_only
    }

    // calculate Wigner D matrix
    double Symmetry_rotation_k::wigner_d(const double beta, const int l, const int m1, const int m2) const
    {
        auto factorial = [](int n) -> int {
            int result = 1;
            for (int i = 1;i <= n;++i) { result *= i;
}
            return result;
            };
        double result = 0.0;
        for (int i = std::max(0, m2 - m1);i <= std::min(l - m1, l + m2);++i) {
            result += std::pow(-1, i) * std::sqrt(factorial(l + m1) * factorial(l - m1) * factorial(l + m2) * factorial(l - m2))
            * std::pow(std::cos(beta / 2), 2 * l + m2 - m1 - 2 * i) * std::pow(-std::sin(beta / 2), m1 - m2 + 2 * i)
            / (factorial(i) * factorial(l - m1 - i) * factorial(l + m2 - i) * factorial(i - m2 + m1));
}
        return result;
    }

    std::complex<double> Symmetry_rotation_k::wigner_D(const TCdouble& euler_angle, const int l, const int m1, const int m2, const bool inv) const
    {
        std::complex<double> prefac(inv ? std::pow(-1, l) : 1, 0);
        return std::exp(-ModuleBase::IMAG_UNIT * static_cast<double>(m1) * euler_angle.x)
            * std::exp(-ModuleBase::IMAG_UNIT * static_cast<double>(m2) * euler_angle.z)
            * wigner_d(euler_angle.y, l, m1, m2) * prefac;
    }

    // c^l_{m1, m2}=<Y_l^m1|S_l^m2>
    std::complex<double> Symmetry_rotation_k::ovlp_Ylm_Slm(const int l, const int m1, const int m2) const
    {
        if (m1 == m2)
        {
            if (m1 == 0) { return 1.0;
}
            if (m1 > 0) { return 1 / std::sqrt(2);
}
            if (m1 < 0) { return std::pow(-1, m1) * ModuleBase::IMAG_UNIT / std::sqrt(2);
}
        }
        else if (m1 == -m2)
        {
            if (m1 > 0) { return -ModuleBase::IMAG_UNIT / std::sqrt(2);
}
            if (m1 < 0) { return std::pow(-1, m1) / std::sqrt(2);
}
        }
        return 0.0;
    }

    // reference: https://github.com/minyez/abf_trans/blob/f9e68e68069a94610d89e077bfe6e8ffac0b097d/src/rotate.cpp#L118
    // because the atom position here is row vector, the original gmatrix(eular angle) is transposed.
    // gmatc: the rotation matrix under the basis of cartesian coordinates
    // gmatc should be a rotation matrix, i.e. det(gmatc)=1
    TCdouble Symmetry_rotation_k::get_euler_angle(const ModuleBase::Matrix3& gmatc) const
    {
        double threshold = this->eps_;
        double alpha = 0.0, beta = 0.0, gamma = 0.0;
        if (std::fabs(gmatc.e32) > threshold || std::fabs(gmatc.e31) > threshold) // sin(beta) is not zero
        {
            // use the 2-angle elements to get alpha and gamma
            alpha = std::atan2(gmatc.e32, gmatc.e31);
            if (alpha < 0) { alpha += 2 * ModuleBase::PI;
}
            gamma = std::atan2(gmatc.e23, -gmatc.e13);
            if (gamma < 0) { gamma += 2 * ModuleBase::PI;
}
            // use the larger one of 2-angle elements to calculate beta
            if (std::fabs(gmatc.e32) > std::fabs(gmatc.e31)) {
                beta = std::atan2(gmatc.e32 / std::sin(alpha), gmatc.e33);
            } else {
                beta = std::atan2(gmatc.e31 / std::cos(alpha), gmatc.e33);
}
        }
        else
        {//sin(beta)=0, beta = 0 or pi, only (alpha+gamma) or (alpha-gamma) is important. now assign this to alpha.
            alpha = std::atan2(gmatc.e12, gmatc.e11);
            if (alpha < 0) { alpha += 2 * ModuleBase::PI;
}
            // if beta=0, gmatc.e11=cos(alpha+gamma), gmatc.e21=sin(alpha+gamma)
            // if beta=pi, gmatc.e11=cos(pi+alpha-gamma), gmatc.e21=sin(pi+alpha-gamma)
            if (gmatc.e33 > 0)
            {
                beta = 0;
                gamma = 0;  //alpha+gamma=alpha => gamma=0
            }
            else
            {
                beta = ModuleBase::PI;
                gamma = ModuleBase::PI;// pi+alpha-gamma=alpha  => gamma=pi
            }
        }
        return TCdouble(alpha, beta, gamma);
    }

    // in: the real value of m in range {-l, -l+1, ..., 0, ..., l-1, l}
    // out: the index of the orbital in a fixed {n， l}, i.e. the index in array [0, 1, -1, 2, -2, ...]
    inline int m2im_k(int m)
    {
        return (m > 0 ? 2 * m - 1 : -2 * m);
    }

    /// T_mm' = [c^\dagger D c]_mm'
    void Symmetry_rotation_k::cal_rotmat_Slm(const ModuleBase::Matrix3* gmatc, const int lmax, const int nop)
    {
        ++this->rotmat_Slm_version_;
        const int nop_tot = (nop < 0) ? this->nsym_ : nop;
        this->rotmat_Slm_.resize(nop_tot);
        // c matrix is independent on isym
        std::vector<ModuleBase::ComplexMatrix> c_mm(lmax + 1);
        for (int l = 0;l <= lmax;++l) {
            c_mm[l].create(2 * l + 1, 2 * l + 1);
}
        for (int l = 0;l <= lmax;++l) {
            for (int m1 = -l;m1 <= l;++m1) {
                for (int m2 = -l;m2 <= l;++m2) {
                    c_mm[l](m2im_k(m1), m2im_k(m2)) = ovlp_Ylm_Slm(l, m1, m2);
}
}
}

        for (int isym = 0;isym < nop_tot;++isym)
        {
            // if R is a reflection operation, calculate D^l(R)=(-1)^l*D^l(IR), so the euler angle of (IR) is needed.
            TCdouble euler_angle = get_euler_angle(gmatc[isym].Det() > 0 ?
                gmatc[isym] : gmatc[isym] * ModuleBase::Matrix3(-1, 0, 0, 0, -1, 0, 0, 0, -1));

            this->rotmat_Slm_[isym].resize(lmax + 1);
            for (int l = 0;l <= lmax;++l)
            {// wigner D matrix
                ModuleBase::ComplexMatrix D_mm(2 * l + 1, 2 * l + 1);
                for (int m1 = -l;m1 <= l;++m1) {
                    for (int m2 = -l;m2 <= l;++m2) {
                        D_mm(m2im_k(m1), m2im_k(m2)) = wigner_D(euler_angle, l, m1, m2, (gmatc[isym].Det() < 0));
}
}
                this->rotmat_Slm_[isym][l] = transpose(c_mm[l], true) * D_mm * c_mm[l];
            }
        }
    }

    void Symmetry_rotation_k::set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
        std::vector<std::complex<double>>& obj_mat, const Parallel_2D& pv, const bool trans) const
    {   // caution: ComplaxMatrix is row-major(col-continuous), but obj_mat is col-major(row-continuous)
        for (int j = 0;j < block.nr;++j) {//outside dimension
            for (int i = 0;i < block.nc;++i) { //inside dimension
                if (pv.in_this_processor(starti + i, startj + j))
                {
                    int index = pv.global2local_col(startj + j) * pv.get_row_size() + pv.global2local_row(starti + i);
                    obj_mat[index] = trans ? block(i, j) : block(j, i);
                }
}
}
    }

    void Symmetry_rotation_k::set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
        std::vector<double>& obj_mat, const Parallel_2D& pv, const bool trans) const
    {   // caution: ComplaxMatrix is row-major(col-continuous), but obj_mat is col-major(row-continuous)
        for (int j = 0;j < block.nr;++j) {//outside dimension
            for (int i = 0;i < block.nc;++i) { //inside dimension
                if (pv.in_this_processor(starti + i, startj + j))
                {
                    int index = pv.global2local_col(startj + j) * pv.get_row_size() + pv.global2local_row(starti + i);
                    obj_mat[index] = trans ? block(i, j).real() : block(j, i).real();
                }
}
}
    }

    // 2d-block parallized rotation matrix in AO-representation, denoted as M.
    // finally we will use D(k)=M(R, k)^\dagger*D(Rk)*M(R, k) to   D(k) from D(Rk) in cal_Ms.
    std::vector<std::complex<double>> Symmetry_rotation_k::contruct_2d_rot_mat_ao(const Symmetry& symm, const Atom* atoms, const Statistics& cell_st,
        const TCdouble& kvec_d_ibz, int isym, const Parallel_2D& pv, const SpinRotation::Su2& spin_U) const
    {
        const bool soc = (this->nspin_ == 4);
        const int npol = soc ? 2 : 1;  // spinor: global AO index is spin-fast interleaved, I = npol*iw_orb + s
        std::vector<std::complex<double>> M_isym(pv.get_local_size(), 0.0);
        // isym >= symm.nrotk addresses the antiunitary coset (spatial part gmatrix_anti[isym-nrotk]),
        // whose atom map lives in a separate table.
        const int nrotk_u = symm.nrotk;
        auto rotated_atom = [&symm, nrotk_u](const int is, const int iat) -> int
            {
                return (is < nrotk_u) ? symm.get_rotated_atom(is, iat)
                                      : symm.get_rotated_atom_anti(is - nrotk_u, iat);
            };
        for (int iat1 = 0;iat1 < cell_st.nat;++iat1)
        {
            int it = cell_st.iat2it[iat1];  // it1=it2
            int ia1 = cell_st.iat2ia[iat1];
            int iat2 = rotated_atom(isym, iat1); //iat2=rot(iat1)
            int ia2 = cell_st.iat2ia[iat2];
            // cal phase factor from return lattice:     exp(-ik_ibz*O)
            double arg = -2 * ModuleBase::PI * kvec_d_ibz * this->irs_.return_lattice_[iat1][isym];
            std::complex<double>phase_factor = std::complex<double>(std::cos(arg), std::sin(arg));
            int iw1start = atoms[it].stapos_wf + ia1 * atoms[it].nw;
            int iw2start = atoms[it].stapos_wf + ia2 * atoms[it].nw;
            int iw = 0;
            while (iw < atoms[it].nw)
            {
                int l = atoms[it].iw2l[iw];
                int nm = 2 * l + 1;
                //caution: the order of m in orbitals may be different from increasing
                if (!soc)
                {
                    set_block_to_mat2d(iw2start + iw, iw1start + iw,
                        phase_factor * this->rotmat_Slm_[isym][l], M_isym, pv, true);
                }
                else
                {
                    // M = T(isym) (x) U(isym): scatter phase * T_l(m,m') * U(a,b) to the interleaved
                    // spinor positions (row = rotated atom/spin, col = original atom/spin). For nspin=4
                    // stapos_wf already carries the npol factor, so the per-atom offset is ia*nw*npol
                    // and the within-atom spinor index is (iw_orb)*npol + spin (spin is the fast index).
                    const int base2 = atoms[it].stapos_wf + ia2 * atoms[it].nw * npol;
                    const int base1 = atoms[it].stapos_wf + ia1 * atoms[it].nw * npol;
                    const ModuleBase::ComplexMatrix& Tl = this->rotmat_Slm_[isym][l];
                    for (int m = 0;m < nm;++m)
                    {
                        for (int mp = 0;mp < nm;++mp)
                        {
                            const std::complex<double> t = phase_factor * Tl(m, mp);
                            for (int a = 0;a < npol;++a)
                            {
                                for (int b = 0;b < npol;++b)
                                {
                                    const int gi = base2 + (iw + m) * npol + a;
                                    const int gj = base1 + (iw + mp) * npol + b;
                                    if (pv.in_this_processor(gi, gj))
                                    {
                                        const int index = pv.global2local_col(gj) * pv.get_row_size() + pv.global2local_row(gi);
                                        // M(isym) = T_l (x) U is the spinor rep, with U = so3_to_su2 placed as-is:
                                        //   M[(m,a),(m',b)] = phase * T_l(m,m') * U_{ab},   U_{ab} = spin_U[a*npol + b].
                                        // Both T_l (rotmat_Slm) and U are ANTI-homomorphisms here (row-vector / R^T convention:
                                        // rotmat_Slm(g)=R_orb(g)^{-1}, so3_to_su2 likewise), so this M is a consistent rep
                                        //  and rot_matrix_ao's stored-DM rotation M^T D M^* is exact for ALL ops.
                                        M_isym[index] = t * spin_U[a * npol + b];
                                    }
                                }
                            }
                        }
                    }
                }
                iw += nm;
            }
        }
        return M_isym;
    }

    // D(k) = M^T(R, k) D(k_ibz) M^*(R, k), if D(k) is col-maj
    // D^T(k) = M^\dagger(R, k) D^T(k_ibz) M(R, k), if D(k) is row-maj
    // Ds from RI_2D_Comm are row-maj
    // the link  ik_ibz-isym-ik can be found in kstars.
    std::vector<std::complex<double>> Symmetry_rotation_k::rot_matrix_ao(const std::vector<std::complex<double>>& DMkibz,
        const int ik_ibz, const int kstar_size, const int isym, const Parallel_2D& pv, const bool TRS_conj) const
    {
        std::vector<std::complex<double>> DMk(pv.nloc, 0.0);
        std::vector<std::complex<double>> DMkibz_M(pv.nloc, 0.0);    // intermediate result
        const char dagger = 'C';
        const char transpose = 'T';
        const char notrans = 'N';
        std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        const int nbasis = pv.get_global_row_size();
        const int i1 = 1;
        if (TRS_conj)
        {
            // D^T* = M^T [M^T (D^T)^T]^\dagger
            ScalapackConnector::gemm(transpose, transpose, nbasis, nbasis, nbasis,
                alpha, this->Ms_[ik_ibz].at(isym).data(), i1, i1, pv.desc, DMkibz.data(), i1, i1, pv.desc,
                beta, DMkibz_M.data(), i1, i1, pv.desc);
            alpha.real(1.0 / static_cast<double>(kstar_size));
            ScalapackConnector::gemm(transpose, dagger, nbasis, nbasis, nbasis,
                alpha, this->Ms_[ik_ibz].at(isym).data(), i1, i1, pv.desc, DMkibz_M.data(), i1, i1, pv.desc,
                beta, DMk.data(), i1, i1, pv.desc);
        }
        else
        {
            // Physical DM rotation D(k) = M^dagger D(k_ibz) M, with M = T (x) U is the anti-homomorphism rep in row-major convention.
            // ABACUS stores the DM transposed (S = D^T), for which this becomes S(gk) = M^T S(k_ibz) M^* = (conj M)^dagger S (conj M)
            // For nspin<4 the orbital-only M is real, so Mc = M and this is bit-identical to the old M^dagger D M.
            const std::vector<std::complex<double>>& Mref = this->Ms_[ik_ibz].at(isym);
            std::vector<std::complex<double>> Mc(Mref.size());
            for (size_t i = 0; i < Mref.size(); ++i) { Mc[i] = std::conj(Mref[i]); }
            ScalapackConnector::gemm(dagger, notrans, nbasis, nbasis, nbasis,
                alpha, Mc.data(), i1, i1, pv.desc, DMkibz.data(), i1, i1, pv.desc,
                beta, DMkibz_M.data(), i1, i1, pv.desc);
            alpha.real(1.0 / static_cast<double>(kstar_size));
            ScalapackConnector::gemm(notrans, notrans, nbasis, nbasis, nbasis,
                alpha, DMkibz_M.data(), i1, i1, pv.desc, Mc.data(), i1, i1, pv.desc,
                beta, DMk.data(), i1, i1, pv.desc);
        }
        return DMk;
    }

    std::vector<std::complex<double>> Symmetry_rotation_k::set_sigma_y_2d(const Parallel_2D& pv) const
    {
        std::vector<std::complex<double>> sigma_y(pv.get_local_size(), 0.0);
        const int nlocal = pv.get_global_row_size();    // = 2*nao for nspin=4
        // sigma_y = [[0, -i], [i, 0]] on the interleaved spin index (I = 2*iorb + spin)
        const std::complex<double> sy[2][2] = { {std::complex<double>(0.0, 0.0), std::complex<double>(0.0, -1.0)},
                                                {std::complex<double>(0.0, 1.0), std::complex<double>(0.0, 0.0)} };
        for (int iorb = 0; 2 * iorb < nlocal; ++iorb)
        {
            for (int a = 0; a < 2; ++a)
            {
                const int b = 1 - a;    // only the off-diagonal spin entries are non-zero
                const int gi = 2 * iorb + a;
                const int gj = 2 * iorb + b;
                if (pv.in_this_processor(gi, gj))
                {
                    const int index = pv.global2local_col(gj) * pv.get_row_size() + pv.global2local_row(gi);
                    sigma_y[index] = sy[a][b];
                }
            }
        }
        return sigma_y;
    }

    std::vector<std::complex<double>> Symmetry_rotation_k::trs_spin_rotate(const std::vector<std::complex<double>>& X,
        const std::vector<std::complex<double>>& sigma_y, const Parallel_2D& pv, const double scale) const
    {
        // stored (transposed 2d-block) form of  D_new = sigma_y * conj(D) * sigma_y  is
        // Sigma_y * conj(X) * Sigma_y  (Sigma_y^T = -Sigma_y, the two minus signs cancel).
        const char notrans = 'N';
        const int nbasis = pv.get_global_row_size();
        const int i1 = 1;
        const std::complex<double> one(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        std::vector<std::complex<double>> Xc(X.size());
        for (size_t i = 0; i < X.size(); ++i) { Xc[i] = std::conj(X[i]); }
        std::vector<std::complex<double>> tmp(pv.get_local_size(), 0.0);
        std::vector<std::complex<double>> out(pv.get_local_size(), 0.0);
        // tmp = Sigma_y * conj(X)
        ScalapackConnector::gemm(notrans, notrans, nbasis, nbasis, nbasis,
            one, sigma_y.data(), i1, i1, pv.desc, Xc.data(), i1, i1, pv.desc,
            beta, tmp.data(), i1, i1, pv.desc);
        // out = scale * tmp * Sigma_y
        ScalapackConnector::gemm(notrans, notrans, nbasis, nbasis, nbasis,
            std::complex<double>(scale, 0.0), tmp.data(), i1, i1, pv.desc, sigma_y.data(), i1, i1, pv.desc,
            beta, out.data(), i1, i1, pv.desc);
        return out;
    }
}

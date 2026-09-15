#pragma once
#include "irreducible_sector.h"
#include "source_base/parallel_2d.h"
#include "source_base/complexmatrix.h"
#include "source_cell/module_symmetry/symm_rot_spin.h"

namespace ModuleSymmetry
{
    /// @brief k-space AO-representation symmetry restoration: reconstructs D(k) at every
    /// k-star member from D(k_ibz), for crystal-symmetry-reduced BZ sampling.
    ///
    /// This is the LibRI-independent subset of what used to be a single
    /// source_lcao/module_ri/module_exx_symmetry/Symmetry_rotation class: everything needed
    /// to go from an irreducible-k-point density matrix to the full-BZ one, so callers that
    /// only need that (e.g. DFT+U's occupation-matrix/DMR restoration) do not have to depend
    /// on module_ri or LibRI. EXX/RPA's own real-space H(R)/RI-coefficient restoration (which
    /// does need RI::Tensor) is built on top of this class in
    /// source_lcao/module_ri/module_exx_symmetry/symm_rotation.h (ModuleSymmetry::Symmetry_rotation,
    /// which inherits from this one).
    class Symmetry_rotation_k
    {
    public:
        Symmetry_rotation_k() {};
        virtual ~Symmetry_rotation_k() {};

        //--------------------------------------------------------------------------------
        // getters
        const std::map<Tap, std::set<TC>>& get_irreducible_sector()const { return this->irs_.get_irreducible_sector(); }
        TCdouble get_return_lattice(const Symmetry& symm,
            const ModuleBase::Matrix3& gmatd, const TCdouble gtransd,
            const TCdouble& posd_a1, const TCdouble& posd_a2)const
        {
            return this->irs_.get_return_lattice(symm, gmatd, gtransd, posd_a1, posd_a2);
        }
        TCdouble get_return_lattice(const int iat, const int isym) const
        {
            return this->irs_.get_return_lattice(iat, isym);
        }
        /// the rotation matrix under the basis of S_l^m. size: [nsym][lmax][nm*nm]
        const std::vector<std::vector<ModuleBase::ComplexMatrix>>& rotmat_Slm = this->rotmat_Slm_;
        const int& abfs_Lmax = this->abfs_Lmax_;
        //--------------------------------------------------------------------------------
        // setters
        void find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::vector<TC>& Rs, const TC& period, const Lattice& lat, const std::string& output_dir = "")
        {
            this->irs_.find_irreducible_sector(symm, atoms, st, Rs, period, lat, output_dir);
        }
        void set_abfs_Lmax(const int l) { this->abfs_Lmax_ = l; }
        //--------------------------------------------------------------------------------
        /// functions  to contruct rotation matrix in AO-representation

        /// The top-level calculation interface of this class. calculate the rotation matrix in AO representation: M
        /// only need once call in each ion step (decided by the configuration)
        /// @param kstars  equal k points to each ibz-kpont, corresponding to a certain symmetry operations.
        /// @param nspin  stored as a member so restore_dm()/contruct_2d_rot_mat_ao() do not each
        ///               need PARAM.inp.nspin (keeps this LibRI-free class free of a module_parameter
        ///               link dependency; every existing caller already has nspin in scope).
        void cal_Ms(const K_Vectors& kv,
            const UnitCell& ucell, const Parallel_2D& pv, const int nspin);

        /// Use calculated M matrix to recover D(k) from D(k_ibz): D(k) = M(R, k)^\dagger D(k_ibz) M(R, k)
        /// the link "ik_ibz-isym-ik" can be found in kstars: k_bz = gmat[isym](k)
        std::vector<std::vector<std::complex<double>>>restore_dm(const K_Vectors& kv,
            const std::vector<std::vector<std::complex<double>>>& dm_k_ibz,
            const Parallel_2D& pv)const;
        std::vector<std::vector<double>>restore_dm(const K_Vectors& kv,
            const std::vector<std::vector<double>>& dm_k_ibz,
            const Parallel_2D& pv)const;
        std::vector<std::complex<double>> rot_matrix_ao(const std::vector<std::complex<double>>& DMkibz,
            const int ik_ibz, const int kstar_size, const int isym, const Parallel_2D& pv, const bool TRS_conj = false) const;

        /// (nspin=4) build the 2*nao spin operator Sigma_y = I_nao (x) sigma_y in 2d-block layout.
        std::vector<std::complex<double>> set_sigma_y_2d(const Parallel_2D& pv) const;

        /// (nspin=4) time-reversal on the spin density matrix: D(k) = sigma_y D^*(-k) sigma_y,
        /// realized distribution-safely as scale * Sigma_y * conj(X) * Sigma_y (X is the already
        /// space-group-rotated D(-k) stored in the transposed 2d-block convention).
        std::vector<std::complex<double>> trs_spin_rotate(const std::vector<std::complex<double>>& X,
            const std::vector<std::complex<double>>& sigma_y, const Parallel_2D& pv, const double scale) const;

        /// calculate Wigner D matrix
        double wigner_d(const double beta, const int l, const int m1, const int m2) const;
        std::complex<double> wigner_D(const TCdouble& euler_angle, const int l, const int m1, const int m2, const bool inv) const;

        /// c^l_{m1, m2}=<Y_l^m1|S_l^m2>
        std::complex<double> ovlp_Ylm_Slm(const int l, const int m1, const int m2) const;

        /// calculate euler angle from rotation matrix
        TCdouble get_euler_angle(const ModuleBase::Matrix3& gmatc) const;

        /// T_mm' = [c^\dagger D c]_mm', the rotation matrix in the representation of real sphere harmonics
        /// @param nop  number of operations in gmatc; <0 means nsym_ (the unitary ones only).
        ///             Pass nsym_+nanti_ to also build the antiunitary operations' T_l.
        void cal_rotmat_Slm(const ModuleBase::Matrix3* gmatc, const int lmax, const int nop);

        /// set a block matrix onto a 2d-parallelized matrix(col-maj), at the position (starti, startj)
        /// if trans=true, the block matrix is transposed before setting
        void set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
            std::vector<std::complex<double>>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;
        void set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
            std::vector<double>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;

        /// 2d-block parallized rotation matrix in AO-representation, denoted as M.
        /// finally we will use D(k)=M(R, k)^\dagger*D(Rk)*M(R, k) to recover D(k) from D(Rk).
        std::vector<std::complex<double>> contruct_2d_rot_mat_ao(const Symmetry& symm, const Atom* atoms, const Statistics& cell_st,
            const TCdouble& kvec_d_ibz, int isym, const Parallel_2D& pv,
            const SpinRotation::Su2& spin_U /*= SpinRotation::Su2{ 1.0, 0.0, 0.0, 1.0 }*/) const;

        std::vector<std::vector<ModuleBase::ComplexMatrix>>& get_rotmat_Slm() { return this->rotmat_Slm_; }

        /// test-only: inject Ms_/little_groups_/nsym_ directly, bypassing cal_Ms(), so restore_dm()
        /// can be unit-tested against synthetic k-stars without a real UnitCell/K_Vectors setup.
        void set_density_rotations_for_testing(const std::vector<std::map<int, std::vector<std::complex<double>>>>& Ms,
            const std::vector<std::vector<int>>& little_groups, const int nsym, const int nspin)
        {
            this->Ms_ = Ms;
            this->little_groups_ = little_groups;
            this->nsym_ = nsym;
            this->nspin_ = nspin;
        }

        //--------------------------------------------------------------------------------
        /// list all cells in a Born-von-Karman supercell of the given period (no LibRI dependency,
        /// unlike RI_Util::get_Born_von_Karmen_cells which this mirrors for 3D periods).
        static std::vector<TC> get_bvk_cells(const TC& period);

    protected:
        /// set by cal_Ms() (or set_density_rotations_for_testing()); avoids a PARAM.inp.nspin
        /// read in restore_dm()/contruct_2d_rot_mat_ao(), which would otherwise pull a
        /// module_parameter link dependency into every target that links this LibRI-free class.
        int nspin_ = 1;

        int nsym_ = 1;
        /// (nspin=4, magnetic) number of ANTIUNITARY elements Theta*g of the Shubnikov group.
        /// Their orbital rotations / return lattices / Ms are appended after the nsym_ unitary
        /// ones, so the raw index isym in [nsym_, nsym_+nanti_) addresses gmatrix_anti[isym-nsym_].
        int nanti_ = 0;
        /// (nspin=4) true when the configuration carries a non-zero local moment. Then pure time
        /// reversal is NOT a symmetry (it reverses m) and the k-star must be restored with the
        /// Shubnikov elements Theta*gmatrix_anti[] instead of the generic -k shortcut.
        bool magnetic_nspin4_ = false;

        double eps_ = 1e-6;

        int abfs_Lmax_ = 0;

        /// the rotation matrix under the basis of S_l^m. size: [nsym][lmax][nm*nm]
        std::vector<std::vector<ModuleBase::ComplexMatrix>> rotmat_Slm_;

        /// bumped every time cal_rotmat_Slm() (re)fills rotmat_Slm_, so derived classes caching a
        /// converted copy of rotmat_Slm_ (e.g. EXX's RI::Tensor mirror) can detect staleness
        /// without recomparing the whole matrix.
        int rotmat_Slm_version_ = 0;

        /// The unitary matrix associate D(Rk) with D(k) for each ibz-kpoint Rk and each symmetry operation.
        /// size: [nks_ibz][nsym][nbasis*nbasis], only need to calculate once.
        std::vector<std::map<int, std::vector<std::complex<double>>>> Ms_;

        /// The little group of each ibz-kpoint: the subset of unitary space-group operations that
        /// fix kvec_d_ibz modulo a reciprocal lattice vector. D(k_ibz) is averaged over this group
        /// before star-expansion, since a finite-grid SCF density need not exactly respect it.
        /// size: [nks_ibz][<=nsym_], always non-empty (identity is always a member).
        std::vector<std::vector<int>> little_groups_;

        /// (nspin=4) the SU(2) spin-1/2 rotation U(isym) for each symmetry operation, size [nsym].
        /// The spinor AO rotation is T(isym) (x) U(isym); restore_HR_nspin4 (EXX) uses it to mix
        /// the 4 spin channels of the real-space H(R). Filled in cal_Ms (identity for nspin<4).
        std::vector<SpinRotation::Su2> spin_U_;

        /// irreducible sector
        Irreducible_Sector irs_;
    };
}

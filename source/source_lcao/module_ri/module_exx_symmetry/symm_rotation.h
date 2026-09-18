#pragma once
#include "source_cell/module_symmetry/symm_rotation_k.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <RI/global/Tensor.h>
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"

namespace ModuleSymmetry
{
    /// Real-space (RI::Tensor / HContainer) H(R) and RI-coefficient symmetry restoration for
    /// EXX/RPA, built on top of the LibRI-independent k-space restoration in
    /// ModuleSymmetry::Symmetry_rotation_k (source_cell/module_symmetry/symm_rotation_k.h),
    /// which provides cal_Ms/restore_dm/rot_matrix_ao and the shared rotation-matrix machinery
    /// (rotmat_Slm_, irs_, Ms_, spin_U_, ...). Only the parts that genuinely need LibRI (RI::Tensor
    /// atom-pair maps, HContainer real-space rotation) live here.
    class Symmetry_rotation : public Symmetry_rotation_k
    {
    public:
        Symmetry_rotation() {};
        ~Symmetry_rotation() {};

        //--------------------------------------------------------------------------------
        // setters
        void set_Cs_rotation(const std::vector<std::vector<int>>& abfs_l_nchi);
        //--------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------
        /// The main functions to rotate matrices
        /// Given H(R) in the irreduceble sector, calculate H(R) for all the atompairs and cells.
        template<typename Tdata>    // RI::Tensor type
        std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>> restore_HR(
            const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_irreduceble)const;
        template<typename TR>   // HContainer type
        void restore_HR(
            const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const hamilt::HContainer<TR>& HR_irreduceble, hamilt::HContainer<TR>& HR_rotated)const;
        /// (nspin=4) spinor overload: rotate all 4 spin channels of H(R) together. On top of the
        /// orbital rotation T1^dagger(.)T2 (mode 'H') / T1^T(.)T2^* (mode 'D') applied to every
        /// channel, the SU(2) spin part U(isym) mixes them:  H'^{ab}=sum_{cd} conj(U_{ca}) U_{db} [T1^dagger H^{cd} T2].
        /// The 4 channels are ordered is=a*2+b (a=row spin, b=col spin), matching RI_2D_Comm::split_is_block.
        /// (nspin=4 magnetic) The atom-pair reduction may also use the ANTIUNITARY elements of the
        /// Shubnikov group, flagged by isym >= nsym_. In real space time reversal acts as
        /// H(R) -> sigma_y H^*(R) sigma_y (R and the orbital indices untouched), which becomes a
        /// remap of the 4 channels applied after the SU(2) mixing; see symm_rotation_r.hpp.
        template<typename Tdata>    // RI::Tensor type
        std::array<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>, 4> restore_HR_nspin4(
            const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const std::array<std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>, 4>& HR_irreducible_soc)const;

        //--------------------------------------------------------------------------------
        /// test functions
        /// test H(R) rotation: giver a full H(R), pick out H(R) in the irreducible sector, rotate it, and compare with the original full H(R)
        template<typename Tdata>    // RI::Tensor type, using col-major implementation
        void test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_full);
        template<typename Tdata>    // test the rotation of RI coefficients
        void test_Cs_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& Cs_full)const;
        template<typename TR>   // HContainer type, using row-major implementation
        void test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const hamilt::HContainer<TR>& HR_full);
        template<typename Tdata>    // HContainer type
        void print_HR(const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR, const std::string name, const double& threshold = 0.0);
        //--------------------------------------------------------------------------------

    private:
        //--------------------------------------------------------------------------------
        std::vector<TC> get_Rs_from_BvK(const K_Vectors& kv)const;
        std::vector<TC> get_Rs_from_adjacent_list(const UnitCell& ucell,
                                                  const Grid_Driver& gd,
                                                  const Parallel_Orbitals& pv) const;
        //--------------------------------------------------------------------------------

        /// The sub functions to rotate matrices
        /// mode='H': H_12(R)=T^\dagger(V)H_1'2'(VR+O_1-O_2)T(V)
        /// mode='D': D_12(R)=T^T(V)D_1'2'(VR+O_1-O_2)T^*(V)
        template<typename Tdata>    // RI::Tensor type, blas
        RI::Tensor<Tdata> rotate_atompair_serial(const RI::Tensor<Tdata>& A, const int isym,
            const Atom& a1, const Atom& a2, const char mode, bool output = false)const;
        template<typename Tdata>    // pointer type, blas
        void rotate_atompair_serial(Tdata* TAT, const Tdata* A, const int& nw1, const int& nw2, const int isym,
            const Atom& a1, const Atom& a2, const char mode)const;
        template<typename TR>    // HContainer type, pblas
        void rotate_atompair_parallel(const TR* Alocal_in, const int isym, const Atom* atoms, const Statistics& st,
            const Tap& ap_in, const Tap& ap_out, const char mode, const Parallel_Orbitals& pv, TR* Alocal_out, const bool output = false)const;

        /// rotate a 3-dim C tensor in RI technique
        template<typename Tdata>
        RI::Tensor<Tdata> rotate_singleC_serial(const RI::Tensor<Tdata>& C, const int isym,
            const Atom& a1, const Atom& a2, const int& type1, bool output = false)const;

        template<typename Tdata>
        RI::Tensor<Tdata> set_rotation_matrix(const Atom& a, const int& isym)const;
        template<typename Tdata>
        RI::Tensor<Tdata> set_rotation_matrix_abf(const int& type, const int& isym)const;

        /// RI::Tensor mirror of rotmat_Slm_ (which is stored as ModuleBase::ComplexMatrix, shared
        /// with the LibRI-free k-space code), rebuilt lazily and cached across the many
        /// set_rotation_matrix/set_rotation_matrix_abf calls within one ion step (one per atom
        /// pair/cell), instead of reconverting the same small block every time.
        const RI::Tensor<std::complex<double>>& get_rotmat_Slm_tensor(const int isym, const int l)const;
        mutable std::vector<std::vector<RI::Tensor<std::complex<double>>>> rotmat_Slm_tensor_;
        mutable int rotmat_Slm_tensor_version_ = -1;
        //--------------------------------------------------------------------------------

        bool reduce_Cs_ = false;

        std::vector<std::vector<int>> abfs_l_nchi_;///< number of abfs for each angular momentum
    };

    template<typename T>  std::string vec3_fmt(const T& x, const T& y, const T& z)
    {
        return  "(" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ")";
    }
    template<typename T>  std::string vec3_fmt(const ModuleBase::Vector3<T>& v)
    {
        return vec3_fmt(v.x, v.y, v.z);
    }
    // output k stars and the rotation matrices of Bloch orbitals
    void print_symrot_info_k(const ModuleSymmetry::Symmetry_rotation& symrot,
        const K_Vectors& kv, const UnitCell& ucell);
    void print_symrot_info_R(const Symmetry_rotation& symrot, const Symmetry& symm,
        const int lmax_ao, const std::vector<TC>& Rs);
}

#include "symm_rotation_r.hpp"
#include "symm_rotation_r_hcontainer.hpp"

#pragma once
#include "irreducible_sector.h"
#include "module_basis/module_ao/parallel_orbitals.h"
#include <RI/global/Tensor.h>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"

namespace ModuleSymmetry
{
    using Tap = std::pair<int, int>;
    using TC = std::array<int, 3>;
    using TapR = std::pair<Tap, TC>;
    using TCdouble = Abfs::Vector3_Order<double>;

    class Symmetry_rotation
    {
    public:
        Symmetry_rotation() {};
        ~Symmetry_rotation() {};

        //--------------------------------------------------------------------------------
        // getters
        const std::map<Tap, std::set<TC>>& get_irreducible_sector()const { return this->irs_.get_irreducible_sector(); }
        TCdouble get_return_lattice(const Symmetry& symm,
            const ModuleBase::Matrix3& gmatd, const TCdouble gtransd,
            const TCdouble& posd_a1, const TCdouble& posd_a2)const
        {
            return this->irs_.get_return_lattice(symm, gmatd, gtransd, posd_a1, posd_a2);
        }
        //--------------------------------------------------------------------------------
        // setters
        void find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::vector<TC>& Rs, const TC& period, const Lattice& lat)
        {
            this->irs_.find_irreducible_sector(symm, atoms, st, Rs, period, lat);
        }
        void set_Cs_rotation(const std::vector<std::vector<int>>& abfs_l_nchi);
        //--------------------------------------------------------------------------------
        /// functions  to contruct rotation matrix in AO-representation

        /// The top-level calculation interface of this class. calculate the rotation matrix in AO representation: M
        /// only need once call in each ion step (decided by the configuration)
        /// @param kstars  equal k points to each ibz-kpont, corresponding to a certain symmetry operations. 
        void cal_Ms(const K_Vectors& kv,
            //const std::vector<std::map<int, TCdouble>>& kstars,
            const UnitCell& ucell, const Parallel_2D& pv);

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

        /// calculate Wigner D matrix
        double wigner_d(const double beta, const int l, const int m1, const int m2) const;
        std::complex<double> wigner_D(const TCdouble& euler_angle, const int l, const int m1, const int m2, const bool inv) const;

        /// c^l_{m1, m2}=<Y_l^m1|S_l^m2>
        std::complex<double> ovlp_Ylm_Slm(const int l, const int m1, const int m2) const;

        /// calculate euler angle from rotation matrix
        TCdouble get_euler_angle(const ModuleBase::Matrix3& gmatc) const;

        /// T_mm' = [c^\dagger D c]_mm', the rotation matrix in the representation of real sphere harmonics
        void cal_rotmat_Slm(const ModuleBase::Matrix3* gmatc, const int lmax);

        /// set a block matrix onto a 2d-parallelized matrix(col-maj), at the position (starti, startj) 
        /// if trans=true, the block matrix is transposed before setting
        void set_block_to_mat2d(const int starti, const int startj, const RI::Tensor<std::complex<double>>& block,
            std::vector<std::complex<double>>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;
        void set_block_to_mat2d(const int starti, const int startj, const RI::Tensor<std::complex<double>>& block,
            std::vector<double>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;

        /// 2d-block parallized rotation matrix in AO-representation, denoted as M.
        /// finally we will use D(k)=M(R, k)^\dagger*D(Rk)*M(R, k) to recover D(k) from D(Rk).
        std::vector<std::complex<double>> contruct_2d_rot_mat_ao(const Symmetry& symm, const Atom* atoms, const Statistics& cell_st,
            const TCdouble& kvec_d_ibz, int isym, const Parallel_2D& pv) const;

        std::vector<std::vector<RI::Tensor<std::complex<double>>>>& get_rotmat_Slm() { return this->rotmat_Slm_; }

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
        std::vector<TC> get_Rs_from_adjacent_list(const UnitCell& ucell, Grid_Driver& gd, const Parallel_Orbitals& pv)const;
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
        //--------------------------------------------------------------------------------

        int nsym_ = 1;

        double eps_ = 1e-6;

        bool TRS_first_ = true; //if R(k)=-k, firstly use TRS to restore D(k) from D(R(k)), i.e conjugate D(R(k)).

        bool reduce_Cs_ = false;
        int abfs_Lmax_ = 0;

        std::vector<std::vector<int>> abfs_l_nchi_;///< number of abfs for each angular momentum

        /// the rotation matrix under the basis of S_l^m. size: [nsym][lmax][nm*nm]
        std::vector<std::vector<RI::Tensor<std::complex<double>>>> rotmat_Slm_;
        // [natom][nsym], phase factor corresponding to a certain kvec_d_ibz
        // std::vector<std::vector<std::complex<double>>> phase_factor_;

        /// The unitary matrix associate D(Rk) with D(k) for each ibz-kpoint Rk and each symmetry operation. 
        /// size: [nks_ibz][nsym][nbasis*nbasis], only need to calculate once.
        std::vector<std::map<int, std::vector<std::complex<double>>>> Ms_;

        /// irreducible sector
        Irreducible_Sector irs_;

    };
}

#include "symmetry_rotation_R.hpp"
#include "symmetry_rotation_R_hcontainer.hpp"      
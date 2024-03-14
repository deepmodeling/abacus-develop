#pragma once
#include "module_base/abfs-vector3_order.h"
#include "module_base/complexmatrix.h"
#include "module_base/matrix3.h"
#include <vector>
#include <map>
#include <set>
#include "module_basis/module_ao/parallel_2d.h"
#include "module_cell/unitcell.h" 
#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/klist.h"
#include <RI/global/Tensor.h>
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
 
namespace ModuleSymmetry
{
    using Tap = std::pair<int, int>;
    using TC = std::array<int, 3>;
    using TapR = std::pair<Tap, TC>;
    using TCdouble = Abfs::Vector3_Order<double>;

    struct ap_less_func
    {
        bool operator()(const Tap& lhs, const Tap& rhs) const
        {
            if (lhs.first < rhs.first)return true;
            else if (lhs.first > rhs.first)return false;
            else return lhs.second < rhs.second;
        }
    };
    struct apR_less_func
    {
        bool operator()(const TapR& lhs, const TapR& rhs) const
        {
            if (lhs.first < rhs.first)return true;
            else if (lhs.first > rhs.first)return false;
            else return lhs.second < rhs.second;
        }
    };
    struct len_less_func
    {
        int norm2(const TC& R)const
        {
            return R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
        }
        bool operator()(const TC& lhs, const TC& rhs) const
        {
            if (norm2(lhs) < norm2(rhs))return true;
            else if (norm2(lhs) > norm2(rhs))return false;
            else return lhs < rhs;
        }
    };

    class Symmetry_rotation
    {
    public:
        Symmetry_rotation() {};
        ~Symmetry_rotation() {};
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

        /// Perfoming {R|t} to atom position r in the R=0 lattice, we get Rr+t, which may get out of R=0 lattice, 
        /// whose image in R=0 lattice is r'=Rr+t-O. This function is to get O for each atom and each symmetry operation.
        /// the range of direct position is [-0.5, 0.5).
        TCdouble get_return_lattice(const Symmetry& symm,
            const ModuleBase::Matrix3& gmatd, const TCdouble gtransd,
            const TCdouble& posd_a1, const TCdouble& posd_a2)const;
        void get_return_lattice_all(const Symmetry& symm, const Atom* atoms, const Statistics& st);

        /// set a block matrix onto a 2d-parallelized matrix(col-maj), at the position (starti, startj) 
        /// if trans=true, the block matrix is transposed before setting
        void set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
            std::vector<std::complex<double>>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;
        void set_block_to_mat2d(const int starti, const int startj, const ModuleBase::ComplexMatrix& block,
            std::vector<double>& obj_mat, const Parallel_2D& pv, const bool trans = false) const;

        /// 2d-block parallized rotation matrix in AO-representation, denoted as M.
        /// finally we will use D(k)=M(R, k)^\dagger*D(Rk)*M(R, k) to recover D(k) from D(Rk).
        std::vector<std::complex<double>> contruct_2d_rot_mat_ao(const Symmetry& symm, const Atom* atoms, const Statistics& cell_st,
            const TCdouble& kvec_d_ibz, int isym, const Parallel_2D& pv) const;

        ModuleBase::Matrix3 direct_to_cartesian(const ModuleBase::Matrix3& d, const ModuleBase::Matrix3& latvec)const;

        std::vector<std::vector<ModuleBase::ComplexMatrix>>& get_rotmat_Slm() { return this->rotmat_Slm_; }

        //--------------------------------------------------------------------------------
        /// The main function to find irreducible sector: {abR}
        void find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::vector<TC>& Rs, const TC& period, const Lattice& lat);
        std::vector<TC> get_Rs_from_BvK(const K_Vectors& kv)const;
        std::vector<TC> get_Rs_from_adjacent_list(const UnitCell& ucell, Grid_Driver& gd, const Parallel_Orbitals& pv)const;
        const std::map<Tap, std::set<TC>>& get_irreducible_sector()const { return this->irreducible_sector_; }
        // const std::map<int, std::set<std::pair<int, TC>>> convirt_irreducible_sector() {};
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

        //--------------------------------------------------------------------------------
        /// test functions
        /// test H(R) rotation: giver a full H(R), pick out H(R) in the irreducible sector, rotate it, and compare with the original full H(R)
        template<typename Tdata>    // RI::Tensor type, using col-major implementation
        void test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR_full);
        template<typename TR>   // HContainer type, using row-major implementation
        void test_HR_rotation(const Symmetry& symm, const Atom* atoms, const Statistics& st, const char mode,
            const hamilt::HContainer<TR>& HR_full);
        template<typename Tdata>    // HContainer type
        void print_HR(const std::map<int, std::map<std::pair<int, TC>, RI::Tensor<Tdata>>>& HR, const std::string name, const double& threshold = 0.0);
        //--------------------------------------------------------------------------------

    private:
        int round2int(const double x)const;

        //--------------------------------------------------------------------------------
        /// The sub functions to find irreducible sector: {abR}

        /// gauge='L' means H(R)=<R|H|0>; gauge='R' means H(R)=<0|H|R>
        /// gauge='L': R'=R+O_1-O_2; gauge='R': R'=R+O_2-O_1
        TC rotate_R_by_formula(const Symmetry& symm, const int isym, const int iat1, const int iat2, const TC& R, const char gauge = 'R')const;
        TapR rotate_R_by_formula(const Symmetry& symm, const int isym, const TapR& iapR, const char gauge = 'R')const;
        /// gauge='L': tau_a + R - tau_b; gauge='R': tau_a - tau_b - R (direct)
        TCdouble get_aRb_direct(const Atom* atoms, const Statistics& st, const int iat1, const int iat2, const TC& R, const char gauge = 'R')const;
        TCdouble get_aRb_direct(const Atom* atoms, const Statistics& st, const int iat1, const int iat2, const TCdouble& R, const char gauge = 'R')const;

        // /// find the irreducible atom pairs
        // /// algorithm 1: the way finding irreducible k-points
        // void find_irreducible_atom_pairs(const Symmetry& symm);
        // /// algorithm 2: taking out atom pairs from the initial set
        // void find_irreducible_atom_pairs_set(const Symmetry& symm);
        // /// double check between the two algorithms
        // void test_irreducible_atom_pairs(const Symmetry& symm);

        void output_full_map_to_irreducible_sector(const int nat);
        void output_sector_star();

        //--------------------------------------------------------------------------------

        /// The sub functions to rotate matrices
        /// mode='H': H_12(R)=T^\dagger(V)H_1'2'(VR+O_1-O_2)T(V)
        /// mode='D': D_12(R)=T^T(V)D_1'2'(VR+O_1-O_2)T^*(V)
        template<typename Tdata>    // RI::Tensor type, blas
        RI::Tensor<Tdata> rotate_atompair_serial(const RI::Tensor<Tdata>& t, const int isym,
            const Atom& a1, const Atom& a2, const char mode, bool output = false)const;
        template<typename TR>    // HContainer type, pblas
        void rotate_atompair_parallel(const TR* Alocal_in, const int isym, const Atom* atoms, const Statistics& st,
            const Tap& ap_in, const Tap& ap_out, const char mode, const Parallel_Orbitals& pv, TR* Alocal_out, const bool output = false)const;
        //--------------------------------------------------------------------------------

        /// The sub functions judge special symmetry
        void gen_symmetry_BvK(const Symmetry& symm, const Atom* atoms, const Lattice& lat, const Statistics& st, const TC bvk_period);
        /// whether in 2D plain or not for each symmetry operation
        std::vector<bool> in_plain(const ModuleSymmetry::Symmetry& symm, const ModuleBase::Matrix3& latvec)const;
        //--------------------------------------------------------------------------------

        int nsym_ = 1;

        double eps_ = 1e-6;

        bool TRS_first_ = true; //if R(k)=-k, firstly use TRS to restore D(k) from D(R(k)), i.e conjugate D(R(k)).

        /// the rotation matrix under the basis of S_l^m. size: [nsym][lmax][nm*nm]
        std::vector<std::vector<ModuleBase::ComplexMatrix>> rotmat_Slm_;
        // [natom][nsym], phase factor corresponding to a certain kvec_d_ibz
        // std::vector<std::vector<std::complex<double>>> phase_factor_;

        /// The unitary matrix associate D(Rk) with D(k) for each ibz-kpoint Rk and each symmetry operation. 
        /// size: [nks_ibz][nsym][nbasis*nbasis], only need to calculate once.
        std::vector<std::map<int, std::vector<std::complex<double>>>> Ms_;

        /// irreducible sector
        /// irreducible atom pairs: [n_iap][(isym, ap=(iat1, iat2))]
        // std::vector<std::map<int, Tap>> atompair_stars_;

        ///The index range of the orbital matrix to be calculated: irreducible R in irreducible atom pairs
        // (including R in other atom pairs that cannot rotate into R_stars_[irreducebule_ap])
        std::map<Tap, std::set<TC>> irreducible_sector_;

        // //[natoms*natoms](R, (isym, irreducible_R))
        // std::vector<std::map<TC, std::pair<int, TC>>> full_map_to_irreducible_sector_;
        // (abR) -> (isym, abR)
        std::map<TapR, std::pair<int, TapR>, apR_less_func> full_map_to_irreducible_sector_;

        // all the {abR}s , where the isym=0 one in each star forms the irreducible sector.
        // [irreducible sector size][isym, ((ab),R)]
        std::vector<std::map<int, TapR>> sector_stars_;

        /// the direct lattice vector of {R|t}\tau-\tau' for each atoms and each symmetry operation. [natom][nsym]
        std::vector<std::vector<TCdouble>> return_lattice_;

        std::vector<int> invmap_;
        std::vector<int> isymbvk_to_isym_;
        std::vector<ModuleBase::Matrix3> bvk_gmatrix_;
        std::vector<ModuleBase::Vector3<double>> bvk_gtrans_;
        int bvk_nsym_;
    };
}

#include "symmetry_rotation_R.hpp"
#include "symmetry_rotation_R_hcontainer.hpp"
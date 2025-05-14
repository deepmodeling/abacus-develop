#pragma once
#include <vector>
#include <map>
#include <set>
#include "module_base/abfs-vector3_order.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell.h" 
#include "module_cell/module_symmetry/symmetry.h"
#include "module_cell/klist.h"

namespace ModuleSymmetry
{
    using Tap = std::pair<int, int>;
    using TC = std::array<int, 3>;
    using TapR = std::pair<Tap, TC>;
    using TCdouble = Abfs::Vector3_Order<double>;

    class Irreducible_Sector
    {
    public:
        struct ap_less_func
        {
            bool operator()(const Tap& lhs, const Tap& rhs) const
            {
                if (lhs.first < rhs.first) {return true;
                } else if (lhs.first > rhs.first) {return false;
                } else { return lhs.second < rhs.second;
}
            }
        };
        struct apR_less_func
        {
            bool operator()(const TapR& lhs, const TapR& rhs) const
            {
                if (lhs.first < rhs.first) {return true;
                } else if (lhs.first > rhs.first) {return false;
                } else { return lhs.second < rhs.second;
}
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
                if (norm2(lhs) < norm2(rhs)) {return true;
                } else if (norm2(lhs) > norm2(rhs)) {return false;
                } else { return lhs < rhs;
}
            }
        };

        //--------------------------------------------------------------------------------
        /// The main function to find irreducible sector: {abR}
        void find_irreducible_sector(const Symmetry& symm, const Atom* atoms, const Statistics& st,
            const std::vector<TC>& Rs, const TC& period, const Lattice& lat);
        const std::map<Tap, std::set<TC>>& get_irreducible_sector()const { return this->irreducible_sector_; }
        // const std::map<int, std::set<std::pair<int, TC>>> convirt_irreducible_sector() {};
        //--------------------------------------------------------------------------------

        /// Perfoming {R|t} to atom position r in the R=0 lattice, we get Rr+t, which may get out of R=0 lattice, 
        /// whose image in R=0 lattice is r'=Rr+t-O. This function is to get O for each atom and each symmetry operation.
        /// the range of direct position is [-0.5, 0.5).
        TCdouble get_return_lattice(const Symmetry& symm,
            const ModuleBase::Matrix3& gmatd, const TCdouble gtransd,
            const TCdouble& posd_a1, const TCdouble& posd_a2)const;

    protected:
        void cal_return_lattice_all(const Symmetry& symm, const Atom* atoms, const Statistics& st);

        //--------------------------------------------------------------------------------
        /// The sub functions to find irreducible sector: {abR}

        /// gauge='L' means H(R)=<R|H|0>; gauge='R' means H(R)=<0|H|R>
        /// gauge='L': R'=R+O_1-O_2; gauge='R': R'=R+O_2-O_1
        TC rotate_R(const Symmetry& symm, const int isym, const int iat1, const int iat2, const TC& R, const char gauge = 'R')const;
        TapR rotate_apR_by_formula(const Symmetry& symm, const int isym, const TapR& iapR, const char gauge = 'R')const;
        /// gauge='L': tau_a + R - tau_b; gauge='R': tau_a - tau_b - R (direct)
        TCdouble get_aRb_direct(const Atom* atoms, const Statistics& st, const int iat1, const int iat2, const TC& R, const char gauge = 'R')const;
        TCdouble get_aRb_direct(const Atom* atoms, const Statistics& st, const int iat1, const int iat2, const TCdouble& R, const char gauge = 'R')const;

        ModuleBase::Matrix3 direct_to_cartesian(const ModuleBase::Matrix3& d, const ModuleBase::Matrix3& latvec)const;

        // /// find the irreducible atom pairs
        // /// algorithm 1: the way finding irreducible k-points
        // void find_irreducible_atom_pairs(const Symmetry& symm);
        // /// algorithm 2: taking out atom pairs from the initial set
        // void find_irreducible_atom_pairs_set(const Symmetry& symm);
        // /// double check between the two algorithms
        // void test_irreducible_atom_pairs(const Symmetry& symm);

        void output_full_map_to_irreducible_sector(const int nat);
        void output_sector_star();
        void write_irreducible_sector();

        //--------------------------------------------------------------------------------
        /// The sub functions judge special symmetry
        void gen_symmetry_BvK(const Symmetry& symm, const Atom* atoms, const Lattice& lat, const Statistics& st, const TC bvk_period);
        /// whether in 2D plain or not for each symmetry operation
        // std::vector<bool> in_plain(const ModuleSymmetry::Symmetry& symm, const ModuleBase::Matrix3& latvec)const;
        //--------------------------------------------------------------------------------

        //--------------------------------------------------------------------------------
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
        std::map<TapR, std::map<int, TapR>> sector_stars_;

        /// the direct lattice vector of {R|t}\tau-\tau' for each atoms and each symmetry operation. [natom][nsym]
        std::vector<std::vector<TCdouble>> return_lattice_;

        std::vector<int> invmap_;

        /// symmetry info for BvK supercell
        std::vector<int> isymbvk_to_isym_;
        int bvk_nsym_;

        friend class Symmetry_rotation;
    };
}
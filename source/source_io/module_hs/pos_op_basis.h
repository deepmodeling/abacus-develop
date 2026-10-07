#ifndef POS_OP_BASIS_H
#define POS_OP_BASIS_H

#include "source_base/sph_bessel_recursive.h"
#include "source_base/vector3.h"
#include "source_base/ylm.h"
#include "source_basis/module_ao/orb_atomic_lm.h"
#include "source_basis/module_ao/orb_gaunt_table.h"
#include "source_basis/module_ao/orb_read.h"
#include "source_cell/unitcell.h"
#include "source_lcao/center2orb_orb11.h"
#include "source_lcao/center2orb_orb21.h"
#include "source_lcao/center2orb.h"

#include <map>
#include <vector>

/**
 * @brief Build and hold the orbital basis and two-center integral tables
 *        needed to evaluate position-operator matrix elements
 *        <phi_mu | r_hat | phi_nu>.
 *
 * The class prepares the numerical atomic orbitals, the auxiliary r-orbital,
 * the non-local projectors (when needed), the spherical-Bessel/Gaunt tables,
 * and the Center2_Orb integral tables. It also builds the iw2* index maps
 * that translate a global orbital index into (type, atom, L, N, m).
 */
class PosOpBasis
{
  public:
    PosOpBasis();
    ~PosOpBasis();

    void build(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force, int nlocal);
    void build_nonlocal(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force, int nlocal);

    // const accessors for PosOpCalc / PosOpWriter
    const Numerical_Orbital_Lm& get_orb_r() const { return orb_r; }

    const std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>>&
    get_center2_orb11() const
    {
        return center2_orb11;
    }

    const std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>>&
    get_center2_orb21_r() const
    {
        return center2_orb21_r;
    }

    const std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>&
    get_center2_orb11_nonlocal() const
    {
        return center2_orb11_nonlocal;
    }

    const std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>&
    get_center2_orb21_r_nonlocal() const
    {
        return center2_orb21_r_nonlocal;
    }

    int get_iw2it(int iw) const { return iw2it[iw]; }
    int get_iw2ia(int iw) const { return iw2ia[iw]; }
    int get_iw2iL(int iw) const { return iw2iL[iw]; }
    int get_iw2iN(int iw) const { return iw2iN[iw]; }
    int get_iw2im(int iw) const { return iw2im[iw]; }

  private:
    void setup_tables(const LCAO_Orbitals& orb);
    void build_orbs(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force);
    void build_nonlocal_orbs(const UnitCell& ucell, const LCAO_Orbitals& orb, bool cal_force);
    void build_iw_map(const UnitCell& ucell, int nlocal);

    std::vector<int> iw2ia;
    std::vector<int> iw2iL;
    std::vector<int> iw2im;
    std::vector<int> iw2iN;
    std::vector<int> iw2it;

    ModuleBase::Sph_Bessel_Recursive::D2* psb_ = nullptr;
    ORB_gaunt_table MGT;

    Numerical_Orbital_Lm orb_r;
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> orbs;
    std::vector<std::vector<Numerical_Orbital_Lm>> orbs_nonlocal;

    std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>>
        center2_orb11;

    std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>>
        center2_orb21_r;

    std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb11>>>>>
        center2_orb11_nonlocal;

    std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, std::map<size_t, Center2_Orb::Orb21>>>>>
        center2_orb21_r_nonlocal;
};

#endif // POS_OP_BASIS_H

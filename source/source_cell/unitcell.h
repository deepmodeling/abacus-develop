#ifndef UNITCELL_H
#define UNITCELL_H

#include <memory>
#include "source_base/global_function.h"
#include "source_cell/sep_cell.h"
#include "source_cell/magnetism.h"
#include "module_symmetry/symmetry.h"
#include "source_cell/module_neighlist/atom_provider.h"
#include "source_cell/base_cell.h"
#include "source_cell/nonlocal_info_base.h"

/**
 * @brief Provide the basic information about unitcell.
 */
class UnitCell : public AtomProvider, public BaseCell {
  public:
    UnitCell();
    ~UnitCell();

    /// @name BaseCell / AtomProvider interface overrides
    /// @{
    double get_lat0() const override {
        return lat0;
    }

    double get_omega() const override {
        return omega;
    }

    const ModuleBase::Matrix3& get_latvec() const override {
        return latvec;
    }

    int get_natom() const override {
        return nat;
    }

    int get_na(int i) const override {
        return atoms[i].na;
    }

    int get_ntype() const override {
        return ntype;
    }

    ModuleBase::Vector3<double> get_tau(int i, int j) const override {
        return atoms[i].tau[j];
    }
    /// @}

    /// @brief Initialize basic cell parameters (latname, ntype, lmaxmax, init_vel)
    ///        from INPUT and parse fixed_axes into lat_axis_free flags.
    void setup_from_input(const std::string& latname_in,
               const int& ntype_in,
               const int& lmaxmax_in,
               const bool& init_vel_in,
               const std::string& fixed_axes_in);

    void setup_cell(const std::string& fn, std::ofstream& log, const double symmetry_prec, 
		    const int dfthalf_type, const std::string& pseudo_dir, const int nspin,
                    const std::string& basis_type, const std::string& orbital_dir, const std::string& init_wfc,
                    const double onsite_radius, const bool deepks_setorb, const bool rpa,
                    const bool fixed_atoms, const bool noncolin, const std::string& calculation, 
		    const std::string& esolver_type, const int symmetry);

    void set_iat2itia();

    void set_iat2iwt(const int& npol_in);

    /// iat2iwt is the atom index iat to the first global index for orbital of
    /// this atom the size of iat2iwt is nat, the value should be
    /// sum_{i=0}^{iat-1} atoms[it].nw * npol where the npol is the number of
    /// polarizations, 1 for non-magnetic(NSPIN=1 or 2), 2 for magnetic(only
    /// NSPIN=4) this part only used for Atomic Orbital based calculation
    /// @brief Indexing tool for find orbital global index from it,ia,iw
    template <typename Tiait>
    inline Tiait
        itiaiw2iwt(const Tiait& it, const Tiait& ia, const Tiait& iw) const {
        return Tiait(this->iat2iwt[this->itia2iat(it, ia)] + iw);
    }
    /// @brief Get iat2iwt
    inline const int* get_iat2iwt() const { return iat2iwt.data(); }
    /// @brief Get npol
    inline const int& get_npol() const { return npol; }

    /// @brief Indexing tools for ia and it
    /// @return true if the last out is reset
    template <typename Tiat, typename Tiait>
    inline bool iat2iait(const Tiat iat, Tiait* ia, Tiait* it) const {
        if (iat >= nat) {
            *ia = 0;
            *it = ntype;
            return false;
        }
        *ia = (Tiait)iat2ia[iat];
        *it = (Tiait)iat2it[iat];
        return true;
    }

    template <typename Tiat, typename Tiait>
    inline bool ijat2iaitjajt(const Tiat ijat,
                              Tiait* ia,
                              Tiait* it,
                              Tiait* ja,
                              Tiait* jt) const {
        Tiat iat = ijat / nat;
        Tiat jat = ijat % nat;
        iat2iait(iat, ia, it);
        iat2iait(jat, ja, jt);
        return true;
    }

    /// @brief Get tau for atom iat
    inline const ModuleBase::Vector3<double>& get_tau(const int& iat) const {
        return atoms[iat2it[iat]].tau[iat2ia[iat]];
    }

    /// @brief Calculate vector between two atoms with R cell
    inline const ModuleBase::Vector3<double>
        cal_dtau(const int& iat1,
                 const int& iat2,
                 const ModuleBase::Vector3<int>& R) const {
        return get_tau(iat2) + double(R.x) * a1 + double(R.y) * a2
               + double(R.z) * a3 - get_tau(iat1);
    }

    /// @brief get atomCounts, which is a map from element type to atom number
    std::map<int, int> get_atom_Counts() const;
    /// @brief get orbitalCounts, which is a map from element type to orbital
    /// number
    std::map<int, int> get_orbital_Counts() const;
    /// @brief get lnchiCounts, which is a map from element type to the l:nchi
    /// map
    std::map<int, std::map<int, int>> get_lnchi_Counts() const;
    /// these are newly added functions, the three above functions are
    /// deprecated and will be removed in the future

  public:

    /// @name Lattice
    /// @{
    Lattice lat;
    std::string& Coordinate = lat.Coordinate;
    std::string& latName = lat.latName;
    double& lat0 = lat.lat0;
    double& lat0_angstrom = lat.lat0_angstrom;
    double& tpiba = lat.tpiba;
    double& tpiba2 = lat.tpiba2;
    double& omega = lat.omega;
    std::vector<int>& lat_axis_free = lat.lat_axis_free;

    ModuleBase::Matrix3& latvec = lat.latvec;
    ModuleBase::Vector3<double>&a1 = lat.a1, &a2 = lat.a2, &a3 = lat.a3;
    ModuleBase::Vector3<double>& latcenter = lat.latcenter;
    ModuleBase::Matrix3& latvec_supercell = lat.latvec_supercell;
    ModuleBase::Matrix3& G = lat.G;
    ModuleBase::Matrix3& GT = lat.GT;
    ModuleBase::Matrix3& GGT = lat.GGT;
    ModuleBase::Matrix3& invGGT = lat.invGGT;
    /// @}

    /// @name Statistics
    /// @{
    Statistics st;
    int& ntype = st.ntype;
    int& nat = st.nat;
    int*& iat2it = st.iat2it;
    int*& iat2ia = st.iat2ia;
    int*& iwt2iat = st.iwt2iat;
    int*& iwt2iw = st.iwt2iw;
    ModuleBase::IntArray& itia2iat = st.itia2iat;
    int& namax = st.namax;
    int& nwmax = st.nwmax;
    /// @}

    Atom* atoms = nullptr;
    Sep_Cell sep_cell;

    /// @name Magnetism
    /// @{
    Magnetism magnet;                               ///< magnetism Yu Liu 2021-07-03
    std::vector<std::vector<double>> atom_mulliken; ///< [nat][nspin]
    int n_mag_at = 0;
    /// @}

    /// @name Symmetry
    /// @{
    ModuleSymmetry::Symmetry symm;
    /// @}

    /// @name Pseudopotential / Orbital parameters
    /// @brief meshx : max number of mesh point in pseudopotential file
    /// @brief natomwfc : number of starting wavefunctions
    /// @brief lmax : Max L used for localized orbital
    /// @brief nmax : Max N used for localized orbital
    /// @brief lmax_ppwf : Max L of pseudo wave functions
    /// @brief lmaxmax : revert from INPUT
    /// @{
    int meshx = 0;
    int natomwfc = 0;
    int lmax = 0;
    int nmax = 0;
    int nmax_total = 0; ///< mohan add 2009-09-10
    int lmax_ppwf = 0;
    int lmaxmax = 0;   ///< liuyu 2021-07-04
    bool init_vel = false; ///< liuyu 2021-07-15
                       // double nelec;
    /// @}

    /// @name File lists
    /// @{
    std::vector<std::string> pseudo_fn;
    std::vector<std::string> pseudo_type;

    std::vector<std::string> orbital_fn;  ///< filenames of orbitals, liuyu add 2022-10-19
    std::string  descriptor_file; ///< filenames of descriptor_file, liuyu add 2023-04-06
    std::vector<std::string> abfs_orbital_files; ///< ABFS orbital filenames read from STRU "ABFS_ORBITAL" (used by LCAO EXX)
    std::vector<std::string> jle_orbital_files;  ///< JLE orbital filenames read from STRU "ABFS_JLES_ORBITAL" (used by LCAO EXX)
    /// @}

    /// @name State flags
    /// @todo Encapsulate ionic_position_updated and cell_parameter_updated with
    /// setters that enforce state invariants; currently exposed as mutable
    /// flags that can be toggled from anywhere.
    /// @{
    bool set_atom_flag = false;                     ///< added on 2009-3-8 by mohan
    bool ionic_position_updated
        = false; ///< whether the ionic position has been updated
    bool cell_parameter_updated
        = false; ///< whether the cell parameters are updated
    /// @}

    /// @name Nonlocal info
    /// @{
    /**
     * @brief Pointer to non-local pseudopotential information.
     *
     * This pointer is set during LCAO initialization and provides access
     * to non-local projector data. It is null for non-LCAO calculations.
     */
    std::unique_ptr<NonlocalInfoBase> infoNL;
    /// @}

  private:
    // --------------------- Private Data ---------------------

    std::vector<int> iat2iwt; ///< iat ==> iwt, the first global index for orbital of this atom
    int npol = 1; ///< number of spin polarizations, initialized in set_iat2iwt
                  /// ----------------- END of iat2iwt part -----------------

    ModuleBase::Matrix3 stress; ///< calculate stress on the cell

    /// @name BaseCell private overrides
    /// @{
    Kind get_kind() const override
    {
        return Kind::unit_cell;
    }

    int get_nat() const override
    {
        return nat;
    }

    const ModuleBase::Matrix3& get_GT() const override
    {
        return GT;
    }
    /// @}
};

#endif // unitcell class

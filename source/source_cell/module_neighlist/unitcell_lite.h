#ifndef UNITCELL_LITE_H
#define UNITCELL_LITE_H

#include "source_cell/module_neighlist/atom_provider.h"
#include <vector>

/**
 * @brief A lightweight unit cell class for molecular dynamics simulations.
 *
 * This class provides a minimal set of unit cell information needed for
 * large-scale molecular dynamics simulations (e.g., billion-atom simulations).
 * It implements the AtomProvider interface and stores only essential data:
 * lattice parameters and atomic coordinates.
 *
 * Compared to the full UnitCell class, UnitCellLite has significantly lower
 * memory overhead by omitting electronic structure-related data such as
 * pseudopotentials, orbitals, magnetism, and symmetry information.
 *
 * @see AtomProvider
 * @see UnitCell
 */
class UnitCellLite : public AtomProvider
{
public:
    /**
     * @brief Default constructor.
     *
     * Initializes all data members to zero/empty state.
     */
    UnitCellLite() = default;

    /**
     * @brief Default destructor.
     */
    ~UnitCellLite() = default;

    // ========== AtomProvider interface implementation ==========

    /**
     * @brief Get the lattice constant in Bohr.
     * @return Lattice constant lat0.
     */
    double get_lat0() const override;

    /**
     * @brief Get the unit cell volume.
     * @return Cell volume omega in Bohr^3.
     */
    double get_omega() const override;

    /**
     * @brief Get the lattice vectors.
     * @return Reference to the 3x3 matrix of lattice vectors.
     */
    const ModuleBase::Matrix3& get_latvec() const override;

    /**
     * @brief Get the total number of atoms.
     * @return Total atom count nat.
     */
    int get_natom() const override;

    /**
     * @brief Get the number of atoms for a given type.
     * @param i Atom type index (0-based).
     * @return Number of atoms of type i.
     * @note Asserts that i is in valid range [0, ntype_).
     */
    int get_na(int i) const override;

    /**
     * @brief Get the number of atom types.
     * @return Number of atom types ntype.
     */
    int get_ntype() const override;

    /**
     * @brief Get the coordinate of atom (type i, index j).
     * @param i Atom type index (0-based).
     * @param j Atom index within type i (0-based).
     * @return Cartesian coordinate of the atom in Bohr.
     * @note Asserts that i and j are in valid ranges.
     */
    ModuleBase::Vector3<double> get_tau(int i, int j) const override;

    // ========== Setter methods ==========

    /**
     * @brief Set the lattice constant.
     * @param lat0 Lattice constant in Bohr.
     */
    void set_lat0(double lat0);

    /**
     * @brief Set the unit cell volume.
     * @param omega Cell volume in Bohr^3.
     */
    void set_omega(double omega);

    /**
     * @brief Set the lattice vectors.
     * @param latvec 3x3 matrix of lattice vectors.
     */
    void set_latvec(const ModuleBase::Matrix3& latvec);

    /**
     * @brief Set all lattice parameters together.
     * @param lat0 Lattice constant in Bohr.
     * @param omega Cell volume in Bohr^3.
     * @param latvec 3x3 matrix of lattice vectors.
     */
    void set_lattice(double lat0, double omega, const ModuleBase::Matrix3& latvec);

    /**
     * @brief Set atom information for all types.
     *
     * This method sets the number of atom types, the count of atoms per type,
     * and all atomic coordinates. It automatically computes the total atom
     * count (nat_) and the cumulative atom counts (naa_).
     *
     * @param ntype Number of atom types.
     * @param na Vector of atom counts for each type [ntype].
     * @param tau Vector of all atomic coordinates [nat].
     *
     * @note Asserts that na.size() == ntype and tau.size() == sum(na).
     */
    void set_atoms(int ntype,
                   const std::vector<int>& na,
                   const std::vector<ModuleBase::Vector3<double>>& tau);

private:
    // ========== Data members ==========

    /// Lattice constant in Bohr
    double lat0_ = 0.0;

    /// Unit cell volume in Bohr^3
    double omega_ = 0.0;

    /// Total number of atoms
    int nat_ = 0;

    /// Number of atom types
    int ntype_ = 0;

    /// Lattice vectors (3x3 matrix)
    ModuleBase::Matrix3 latvec_;

    /// Number of atoms for each type [ntype]
    std::vector<int> na_;

    /// Cumulative sum of na: naa_[i] = na_[0] + na_[1] + ... + na_[i]
    std::vector<int> naa_;

    /// Atomic coordinates in Cartesian (Bohr) [nat]
    std::vector<ModuleBase::Vector3<double>> tau_;

    // ========== Internal methods ==========

    /**
     * @brief Compute cumulative atom counts from na_.
     *
     * Updates naa_ such that naa_[i] = sum of na_[0] to na_[i].
     * Called internally by set_atoms().
     */
    void compute_naa_();
};

#endif // UNITCELL_LITE_H
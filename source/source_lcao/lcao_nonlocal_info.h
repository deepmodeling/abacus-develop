#ifndef LCAO_NONLOCAL_INFO_H
#define LCAO_NONLOCAL_INFO_H

#include "../source_cell/nonlocal_info_base.h"
#include "setup_nonlocal.h"

/**
 * @brief LCAO-specific implementation of non-local pseudopotential information.
 *
 * Adapts the InfoNonlocal class to the NonlocalInfoBase interface, enabling
 * the UnitCell module to access LCAO non-local projector data without direct
 * dependency on LCAO-specific implementations.
 */
class LCAONonlocalInfo : public NonlocalInfoBase {
    InfoNonlocal nonlocal;

public:
    /**
     * @brief Default constructor.
     */
    LCAONonlocalInfo() = default;

    /**
     * @brief Destructor.
     */
    ~LCAONonlocalInfo() = default;

    /**
     * @brief Get the maximum cutoff radius among all non-local projectors.
     * @return const reference to rcutmax_Beta.
     */
    const double& get_rcutmax_Beta() const override {
        return nonlocal.get_rcutmax_Beta();
    }

    /**
     * @brief Get the number of projectors for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return Number of projectors.
     */
    int get_nproj(const int& type_in) const override {
        return nonlocal.nproj[type_in];
    }

    /**
     * @brief Get the maximum number of projectors across all atom types.
     * @return Maximum nproj value.
     */
    int get_nprojmax() const override {
        return nonlocal.nprojmax;
    }

    /**
     * @brief Get the cutoff radius for a specific atom type's projectors.
     * @param[in] type_in Atom type index.
     * @return Cutoff radius.
     */
    double get_rcut_max(const int& type_in) const override {
        return nonlocal.Beta[type_in].get_rcut_max();
    }

    /**
     * @brief Get the element label for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return const reference to label string.
     */
    const std::string& get_label(const int& type_in) const override {
        return nonlocal.Beta[type_in].getLabel();
    }

    /**
     * @brief Get the type index for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return Type index.
     */
    int get_type(const int& type_in) const override {
        return nonlocal.Beta[type_in].getType();
    }

    /**
     * @brief Get the angular momentum L for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Angular momentum L.
     */
    int get_proj_L(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getL();
    }

    /**
     * @brief Get the number of radial mesh points for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Number of radial mesh points.
     */
    int get_proj_Nr(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getNr();
    }

    /**
     * @brief Get the radial mesh array for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return const pointer to radial mesh array.
     */
    const double* get_proj_radial(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getRadial();
    }

    /**
     * @brief Get the beta radial function array for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return const pointer to beta_r array.
     */
    const double* get_proj_beta_r(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getBeta_r();
    }

    /**
     * @brief Get the number of k-space mesh points for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Number of k-space mesh points.
     */
    int get_proj_Nk(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getNk();
    }

    /**
     * @brief Get the k-space spacing for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Delta k value.
     */
    double get_proj_dk(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getDk();
    }

    /**
     * @brief Get the uniform real-space spacing for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Delta r uniform value.
     */
    double get_proj_dr_uniform(const int& type_in, const int& ip_in) const override {
        return nonlocal.Beta[type_in].Proj[ip_in].getDruniform();
    }

    /**
     * @brief Setup non-local projectors for LCAO basis.
     * @param[in] ntype_in Number of atom types.
     * @param[in] atoms_in Pointer to atoms array.
     * @param[in] log Reference to log file stream.
     * @param[in] orb Reference to LCAO orbitals object.
     */
    void setupNonlocal(
        const int& ntype_in,
        Atom* atoms_in,
        std::ofstream& log,
        LCAO_Orbitals& orb,
        const std::string& basis_type,
        const bool& out_element_info,
        const bool& lspinorb,
        const int& nspin) {
        nonlocal.setupNonlocal(ntype_in, atoms_in, log, orb, basis_type, out_element_info, lspinorb, nspin);
    }

    /**
     * @brief Get non-const reference to internal InfoNonlocal object.
     * Used for special operations that require direct access.
     * @return Reference to internal InfoNonlocal.
     */
    InfoNonlocal& get_nonlocal() {
        return nonlocal;
    }

    /**
     * @brief Get const reference to internal InfoNonlocal object.
     * Used for read-only access in special cases.
     * @return Const reference to internal InfoNonlocal.
     */
    const InfoNonlocal& get_nonlocal() const {
        return nonlocal;
    }
};

#endif
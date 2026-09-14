#ifndef NONLOCAL_INFO_BASE_H
#define NONLOCAL_INFO_BASE_H

#include <string>

/**
 * @brief Abstract base class for non-local pseudopotential information.
 *
 * Provides a common interface for different basis set implementations
 * to store and access non-local projector data. This class enables
 * the UnitCell module to be independent of LCAO-specific implementations
 * by using polymorphism.
 */
class NonlocalInfoBase {
public:
    /**
     * @brief Virtual destructor for proper cleanup of derived classes.
     */
    virtual ~NonlocalInfoBase() = default;

    /**
     * @brief Get the maximum cutoff radius among all non-local projectors.
     * @return const reference to rcutmax_Beta.
     */
    virtual const double& get_rcutmax_Beta() const = 0;

    /**
     * @brief Get the number of projectors for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return Number of projectors.
     */
    virtual int get_nproj(const int& type_in) const = 0;

    /**
     * @brief Get the maximum number of projectors across all atom types.
     * @return Maximum nproj value.
     */
    virtual int get_nprojmax() const = 0;

    /**
     * @brief Get the cutoff radius for a specific atom type's projectors.
     * @param[in] type_in Atom type index.
     * @return Cutoff radius.
     */
    virtual double get_rcut_max(const int& type_in) const = 0;

    /**
     * @brief Get the element label for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return const reference to label string.
     */
    virtual const std::string& get_label(const int& type_in) const = 0;

    /**
     * @brief Get the type index for a specific atom type.
     * @param[in] type_in Atom type index.
     * @return Type index.
     */
    virtual int get_type(const int& type_in) const = 0;

    /**
     * @brief Get the angular momentum L for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Angular momentum L.
     */
    virtual int get_proj_L(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the number of radial mesh points for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Number of radial mesh points.
     */
    virtual int get_proj_Nr(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the radial mesh array for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return const pointer to radial mesh array.
     */
    virtual const double* get_proj_radial(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the beta radial function array for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return const pointer to beta_r array.
     */
    virtual const double* get_proj_beta_r(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the number of k-space mesh points for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Number of k-space mesh points.
     */
    virtual int get_proj_Nk(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the k-space spacing for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Delta k value.
     */
    virtual double get_proj_dk(const int& type_in, const int& ip_in) const = 0;

    /**
     * @brief Get the uniform real-space spacing for a specific projector.
     * @param[in] type_in Atom type index.
     * @param[in] ip_in Projector index.
     * @return Delta r uniform value.
     */
    virtual double get_proj_dr_uniform(const int& type_in, const int& ip_in) const = 0;
};

#endif
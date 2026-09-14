/**
 * @file atom_pseudo.h
 * @brief Atom_pseudo class for atom pseudopotential data.
 */
#ifndef ATOM_PSEUDO_H
#define ATOM_PSEUDO_H

#include "source_base/vector3.h"
#include "source_base/complexarray.h"
#include "source_base/complexmatrix.h"
#include "pseudo.h"

/**
 * @brief Atom_pseudo class for atom pseudopotential data.
 */
class Atom_pseudo : public pseudo
{
public:

    Atom_pseudo();
    ~Atom_pseudo();

    /// @brief spin-orbit coupling data (mohan add 2021-05-07)
    ModuleBase::ComplexArray d_so; ///< (:,:,:), spin-orbit case
    ModuleBase::matrix d_real;     ///< (:,:), non-spin-orbit case
    int nproj;                     ///< number of projectors
    int nproj_soc;                 ///< dimension of D_ij^so

    std::vector<int> non_zero_count_soc = {0, 0, 0, 0}; ///< non-zero count for SOC
    std::vector<std::vector<int>> index1_soc = {{}, {}, {}, {}}; ///< index1 for SOC
    std::vector<std::vector<int>> index2_soc = {{}, {}, {}, {}}; ///< index2 for SOC

    /**
     * @brief Set spin-orbit coupling matrix.
     *
     * @param d_so_in input complex matrix for SOC
     * @param nproj_in number of projectors
     * @param nproj_in_so number of projectors for SOC
     * @param has_so whether SOC is present
     * @param lspinorb whether spin-orbit is enabled
     * @param nspin number of spin states
     */
    void set_d_so(
        ModuleBase::ComplexMatrix &d_so_in,
        const int &nproj_in,
        const int &nproj_in_so,
        const bool has_so,
        const bool lspinorb,
        const int nspin);


    /**
     * @brief Get spin-orbit coupling matrix.
     *
     * @param is spin index
     * @param p1 first projector index
     * @param p2 second projector index
     * @param tmp_d pointer to the matrix element (output)
     */
    inline void get_d(const int& is, const int& p1, const int& p2, const std::complex<double>*& tmp_d)
    {
        tmp_d = &this->d_so(is, p1, p2);
        return;
    }

    /**
     * @brief Get real coupling matrix.
     *
     * @param is spin index
     * @param p1 first projector index
     * @param p2 second projector index
     * @param tmp_d pointer to the matrix element (output)
     */
    inline void get_d(const int& is, const int& p1, const int& p2, const double*& tmp_d)
    {
        tmp_d = &this->d_real(p1, p2);
        return;
    }
    

#ifdef __MPI
    /**
     * @brief Broadcast atom pseudopotential data.
     *
     * For UPF201 format.
     */
    void bcast_atom_pseudo(void);
#endif

};

#endif

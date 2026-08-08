#ifndef MAGNETISM_H
#define MAGNETISM_H

#include "source_base/global_function.h"
#include "source_base/vector3.h"
#include <vector>

/**
 * @brief Class for magnetism calculations.
 */
class Magnetism
{
public:
    /// @brief Constructor
    Magnetism();
    /// @brief Destructor
    ~Magnetism();

    /// @brief notice : bcast (MPI operation) is done in unitcell
    std::vector<double> start_mag;

    /// @brief tot_mag : majority spin - minority spin (nelup - neldw)
    double tot_mag;

    /// @brief non-collinear total magnetic moment
    double tot_mag_nc[3]={0.0};

    /// @brief absolute magnetic moment
    double abs_mag;

    /**
     * @brief Compute the magnetic moment.
     *
     * @param omega unit cell volume
     * @param nrxx number of grid points in real space
     * @param nxyz total number of grid points
     * @param rho charge density
     * @param nspin number of spin states
     * @param two_fermi whether to use two Fermi levels
     * @param nelec total number of electrons
     * @param nelec_spin number of electrons per spin channel
     */
    void compute_mag(const double& omega,
            const int& nrxx,
            const int& nxyz,
            const double* const * rho,
            const int& nspin,
            const bool& two_fermi,
            const double& nelec,
            double* nelec_spin);

    /// @brief ux_
    double ux_[3]={0.0};

    /// @brief lsign_
    bool lsign_=false;

private:
    /**
     * @brief Judge if two vectors are parallel.
     *
     * @param a first vector
     * @param b second vector
     * @return true if vectors are parallel
     */
    bool judge_parallel(const double a[3], const ModuleBase::Vector3<double> &b);

};

/**
 * @brief A comment about variables nelup, neldw, multiplicity and tot_mag.
 *
 * All these variables contain the same information and must be kept harmonized.
 * Variables nelup and neldw will be removed in future versions of the code.
 * Variables multiplicity and tot_mag, though redundent will probably
 * coexist since multiplicity is the more natural way (?)for defining the spin
 * configuratio in the quantum-chemistry community while tot_mag is
 * more natural (?) when dealing with extended systems.
 */

#endif

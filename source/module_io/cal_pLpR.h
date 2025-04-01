/**
 * calculate the <phi_i|Lx/Ly/Lz|phi_j> matrix elements, in which the Lx/Ly/Lz
 * are the angular momentum operators, |phi_i> and |phi_j> are the numerical
 * atomic orbitals (NAOs).
 * 
 * Formulation
 * -----------
 * 
 * Calculate the <phi_i|Lx/Ly/Lz|phi_j> with ladder operator L+ and L-.
 * 
 * The relation between Lx, Ly and L+, L- are:
 * 
 *                                  Lx = (L+ + L-) / 2
 *                                  Ly = (L+ - L-) / 2i
 * 
 * With L+, the spherical harmonic function Ylm (denoted as |l, m> in the following)
 * can be raised:
 * 
 *                          L+|l, m> = sqrt((l-m)(l+m+1))|l, m+1>
 * 
 * Likely, with L-, the spherical harmonic function Ylm can be lowered:
 * 
 *                          L-|l, m> = sqrt((l+m)(l-m+1))|l, m-1>
 * 
 * Therefore the Lx matrix element can be calculated as:
 * 
 *                <l, m|Lx|l, m'> =   sqrt((l-m)(l+m+1)) * delta(m, m'+1) / 2
 *                                  + sqrt((l+m)(l-m+1)) * delta(m, m'-1) / 2
 * 
 * The Ly matrix element can be calculated as:
 * 
 *                <l, m|Ly|l, m'> =   sqrt((l-m)(l+m+1)) * delta(m, m'+1) / 2i
 *                                  - sqrt((l+m)(l-m+1)) * delta(m, m'-1) / 2i
 * 
 * The Lz matrix element can be calculated as:
 * 
 *                          <l, m|Lz|l, m'> = m * delta(m, m')
 * 
 * However, things will change when there are more than one centers.
 * 
 * Technical Details
 * -----------------
 * 
 * 0. The calculation of matrix elements involves the two-center-integral calculation,
 *    this is supported by ABACUS built-in class TwoCenterIntegrator.
 *    see: source/module_basis/module_nao/two_center_integrator.h.
 *  
 * 1. The interface of it is RadialCollection, which is a collection of radial functions. 
 *    see: source/module_basis/module_nao/radial_collection.h
 * 
 * 2. The radial functions are stored in AtomicRadials class,
 *    see: source/module_basis/module_nao/atomic_radials.h
 * 
 * 3. The construction of AtomicRadials involves the filename of orbital, it is stored
 *    in the UnitCell instance
 * 
 * 4. The calculation will run over all supercells in which the two-center-integral
 *    is not zero. This is done by the class SltkGridDriver, which is a driver for
 *    searching the neighboring cells.
 */
#include <cmath>
#include <vector>
#include <map>
#include <tuple>
#include <complex>
#include <memory>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"

namespace ModuleIO
{
    std::complex<double> cal_LzijR(
        const std::unique_ptr<TwoCenterIntegrator>& calculator,
        const int it, const int ia, const int il, const int iz, const int mi,
        const int jt, const int ja, const int jl, const int jz, const int mj,
        const ModuleBase::Vector3<double>& vR);

    std::complex<double> cal_LyijR(
        const std::unique_ptr<TwoCenterIntegrator>& calculator,
        const int it, const int ia, const int il, const int iz, const int im,
        const int jt, const int ja, const int jl, const int jz, const int jm,
        const ModuleBase::Vector3<double>& vR);

    std::complex<double> cal_LxijR(
        const std::unique_ptr<TwoCenterIntegrator>& calculator,
        const int it, const int ia, const int il, const int iz, const int im,
        const int jt, const int ja, const int jl, const int jz, const int jm,
        const ModuleBase::Vector3<double>& vR);
    
    // the calculation of <phi_i|Lx/Ly/Lz|phi_j> matrix elements will be outputted
    // in the way that indexed by:
    // it, ia, il, iz, im, iRx, iRy, iRz, jt, ja, jl, jz, jm, in which the
    // iRx, iRy, iRz are the indices of the supercell in which the two-center-integral
    // it and jt are indexes of atomtypes,
    // ia and ja are indexes of atoms within the atomtypes,
    // il and jl are indexes of the angular momentum,
    // iz and jz are indexes of the zeta functions
    // im and jm are indexes of the magnetic quantum numbers.
    // The output is a complex number, which is the value of the matrix element.
    // Always the matrix is quite large, so direct print to file.
    class AngularMomentumExpectationCalculator
    {
        public:
            AngularMomentumExpectationCalculator() = delete; // the default constructor is meaningless
            AngularMomentumExpectationCalculator(
                const std::string& orbital_dir,
                const UnitCell& ucell,
                const double& search_radius,
                const int tdestructor,
                const int tgrid,
                const int tatom,
                const bool searchpbc,
                std::ofstream* ptr_log = nullptr);
            ~AngularMomentumExpectationCalculator() = default;

            void calculate(std::ofstream* ofs, 
                           const UnitCell& ucell, 
                           const char dir = 'x',
                           const int precision = 10);

        private:
            // ofsrunning
            std::ofstream* ofs_;
            // the two-center-integrator
            std::unique_ptr<TwoCenterIntegrator> calculator_;
            // the spherical bessel transformer
            ModuleBase::SphericalBesselTransformer sbt_;
            // the radial collection
            std::unique_ptr<RadialCollection> orb_;

            // neighboring searcher
            std::unique_ptr<Grid_Driver> neighbor_searcher_;
    };
} // namespace ModuleIO
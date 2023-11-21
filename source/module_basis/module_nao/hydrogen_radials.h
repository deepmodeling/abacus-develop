#ifndef HYDROGEN_RADIALS_H_
#define HYDROGEN_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"
#include "module_base/assoc_laguerre.h"
// include pair container
#include <utility>
// include map container
#include <map>

/*
    Hydrogen-like radial functions generator

    This class is a derived class of RadialSet, the key function to call is build().
    With build(), the hydrogen-like radial functions will be generated and stored in data member NumericalRadials.
    The workflow is like:

    1. build() calls hydrogen() to start the hydrogen-like radial functions generation task.
    2. hydrogen() calls generate_orb() to generate the radial functions for each n, l pair.
    3. generate_orb() calls generate_hydrogen_radial_toconv() to generate the radial function for set of n, l pairs.
    4. with more details, unzip_strategy() is called to parse the strategy string to get the n, l pairs.
    5. for each n, l pair, generate_hydrogen_radial_toconv() is called to generate the radial function for a given n, l pair.
    6. note that the "toconv" means each radial function should have a norm converged to 1 at a certain radius, this is
       also controlled by a convergence threshold (conv_thr), now set to 1e-6. Therefore each radial function will have different
       rmax.
    7. if the norm of radial function is not converged to 1 at the given radius, the radius will be increased by dr, and the
       radial function will be recalculated until the norm is converged to 1.
    8. the radial function is calculated by generate_hydrogen_radial_segment(), which is a wrapper of assoc_laguerre_.
    9. assoc_laguerre_ is a class that can calculate the radial function of hydrogen-like atom, with the help of associated
       Laguerre polynomials.
    10. the radial function is calculated from 0.0 to a radius where the norm of radial function is converged, and the radial
        function is stored in a pair of vectors, first vector is the radial grid, second vector is the radial function.
    11. after the first generation step, all radial functions are zero-padded to the same length, and stored in a map, with
        the key being the n, l pair, and the value being the pair of vectors.
    12. to store the radial functions, they must indiced by the l, zeta pairs, so mapping_nl_lzeta() is called to map the n, l
        pairs to the l, zeta pairs.
    13. finally, the radial functions are stored in NumericalRadials, with the key being the l, zeta pairs, and the value being
        the pair of vectors.
    
    User capable settings:
    1. charge of the nucleus charge according to pseudopotential, it is also a rescale of radius. Higher charge result in smaller
       radius.
    2. generation strategy. minimal will only generate 1 orbital per angular momentum, double will generate 2 orbitals per angular
       momentum, full will generate all orbitals per angular momentum up to nmax. Note: nmax is read from atom_database.
    3. conv_thr, the convergence threshold of the norm of radial function, if not reached, will continue to increase the radius.
       user use this to control the accuracy of radial function. Too large conv_thr will result in inaccurate spherical Bessel
       transformation results (in FFT-two_center_integrator)
*/
class HydrogenRadials : public RadialSet
{
    public:
        /// @brief default constructor
        HydrogenRadials() {}
        /// @brief overloaded assignment operator
        /// @param rhs HydrogenRadials object to be assigned
        /// @return a reference to the assigned HydrogenRadials object
        HydrogenRadials& operator=(const HydrogenRadials& rhs);
        /// @brief copy constructor
        /// @return a copy of the HydrogenRadials object
        HydrogenRadials* clone() const { return new HydrogenRadials(*this); } // covariant return type
        /// @brief destructor
        ~HydrogenRadials() {}

        /// @brief build the hydrogen-like radial functions and push into NumericalRadials
        /// @param itype index of the atom type
        /// @param charge charge of the nucleus
        /// @param nmax maxmium principal quantum number
        /// @param rcut cutoff radius of the radial function (not used anymore)
        /// @param dr step size of the radial grid
        /// @param rank MPI rank
        /// @param strategy strategy string
        /// @param ptr_log pointer to the log ofstream
        void build(const int itype = 0,
                   const double charge = 1.0,
                   const int nmax = 0,
                   const double rcut = 10.0,
                   const double dr = 0.01,
                   const double conv_thr = 1e-6,
                   const int rank = 0,
                   const std::string strategy = "minimal",
                   std::ofstream* ptr_log = nullptr        
        );

        /// @brief parse the strategy string to get the n, l pairs
        /// @param nmax maxmium principal quantum number
        /// @param strategy strategy string
        /// @return a vector of n, l pairs
        std::vector<std::pair<int, int>> unzip_strategy(const int nmax = 0,
                                                        const std::string strategy = "minimal");
        /// @brief generate hydrogen-like radial functions for a given n, l, from 0.0 to a radius where the norm of radial function is converged
        /// @param charge charge of the nucleus
        /// @param n principal quantum number
        /// @param l angular momentum quantum number
        /// @param converge_threshold the threshold of norm of radial function, if not reached, will continue to increase the radius
        /// @param rank MPI rank
        /// @param rgrid returned radial grid
        /// @param rvalue returned radial function
        /// @param ptr_log pointer to the log ofstream
        /// @return the rmax of present radial function
        double generate_hydrogen_radial_toconv(const double charge,
                                               const int n,
                                               const int l,
                                               const double conv_thr,
                                               const int rank,
                                               std::vector<double>& rgrid,
                                               std::vector<double>& rvalue,
                                               std::ofstream* ptr_log = nullptr);
        /// @brief returns the norm of the radial function
        /// @param rgrid radial grid
        /// @param rvalue radial function
        /// @return norm of the radial function
        double radial_norm(const std::vector<double> rgrid,
                           const std::vector<double> rvalue);

        /// @brief generate set of hydrogen-like radial functions for a given charge, nmax, dr, rank, strategy
        /// @param charge charge of the nucleus
        /// @param nmax maxmium principal quantum number
        /// @param dr step size of the radial grid
        /// @param rank MPI rank
        /// @param strategy strategy string
        /// @param ptr_log pointer to the log ofstream
        std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>>
        generate_orb(const double charge = 1.0,
                     const int nmax = 0,
                     const double dr = 0.01,
                     const double conv_thr = 1e-6,
                     const int rank = 0,
                     const std::string strategy = "minimal",
                     std::ofstream* ptr_log = nullptr);
        /// @brief mapping the n, l pairs to the l, zeta pairs
        /// @param nmax maxmium principal quantum number
        /// @param strategy strategy string
        /// @return a map of n, l pairs to l, zeta pairs
        std::map<std::pair<int, int>, std::pair<int, int>>
        mapping_nl_lzeta(const int nmax = 0,
                         const std::string strategy = "minimal");
        /// @brief kernel function of hydrogen-like radial functions
        /// @param charge charge of the nucleus
        /// @param nmax maxmium principal quantum number
        /// @param dr step size of the radial grid
        /// @param conv_thr convergence threshold of the norm of radial function
        /// @param rank MPI rank
        /// @param strategy strategy string
        /// @param ptr_log pointer to the log ofstream
        void hydrogen(const double charge = 1.0,
                      const int nmax = 0,
                      const double dr = 0.01,
                      const double conv_thr = 1e-6,
                      const int rank = 0,
                      const std::string strategy = "minimal",
                      std::ofstream* ptr_log = nullptr);

    private:
        /// @brief [DEPRECATED IN ACTUAL USE!!!] generate hydrogen-like radial functions for a given charge, nmax, dr, rank, strategy
        /// @param charge charge of the nucleus
        /// @param nmax maxmium principal quantum number
        /// @param rcut the maximal radius
        /// @param dr step size of the radial grid
        /// @param rank MPI rank
        /// @param ptr_log pointer to the log ofstream
        void generate_hydrogen_radials(const double charge = 1.0,
                                       const int nmax = 0,
                                       const double rcut = 10.0,
                                       const double dr = 0.01,
                                       const int rank = 0,
                                       std::ofstream* ptr_log = nullptr);
        /// @brief generate hydrogen-like radial functions for a given n, l, in a given range [rmin, rmax]
        /// @param charge charge of the nucleus
        /// @param n principal quantum number
        /// @param l angular momentum quantum number
        /// @param rmin the minimal radius
        /// @param rmax the maximal radius
        /// @param dr step size of the radial grid
        /// @param rank MPI rank
        /// @param ptr_log pointer to the log ofstream
        /// @return the radial function stored in std::vector<double>
        std::vector<double> generate_hydrogen_radial_segment(const double charge = 1.0,
                                                             const int n = 0,
                                                             const int l = 0,
                                                             const double rmin = 0.0,
                                                             const double rmax = 10.0,
                                                             const double dr = 0.01,
                                                             const int rank = 0,
                                                             std::ofstream* ptr_log = nullptr);
        /// @brief Associated Laguerre polynomials generator
        Assoc_Laguerre assoc_laguerre_;
};
#endif // HYDROGEN_RADIALS_H_
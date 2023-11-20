#ifndef HYDROGEN_RADIALS_H_
#define HYDROGEN_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"
#include "module_base/assoc_laguerre.h"
// include pair container
#include <utility>

class HydrogenRadials : public RadialSet
{
    public:
        HydrogenRadials() {}
        
        HydrogenRadials& operator=(const HydrogenRadials& rhs);
        HydrogenRadials* clone() const { return new HydrogenRadials(*this); } // covariant return type
        
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
        void generate_hydrogen_radials(const double charge = 1.0,
                                       const int nmax = 0,
                                       const double rcut = 10.0,
                                       const double dr = 0.01,
                                       const int rank = 0,
                                       std::ofstream* ptr_log = nullptr);
        /// @brief generate hydrogen-like radial functions for a given n, l, in a given range [rmin, rmax]
        /// @param charge 
        /// @param n 
        /// @param l 
        /// @param rmin 
        /// @param rmax 
        /// @param dr 
        /// @param rank 
        /// @param ptr_log 
        /// @return 
        std::vector<double> generate_hydrogen_radial_segment(const double charge = 1.0,
                                                             const int n = 0,
                                                             const int l = 0,
                                                             const double rmin = 0.0,
                                                             const double rmax = 10.0,
                                                             const double dr = 0.01,
                                                             const int rank = 0,
                                                             std::ofstream* ptr_log = nullptr);

        Assoc_Laguerre assoc_laguerre_;
};
#endif // HYDROGEN_RADIALS_H_
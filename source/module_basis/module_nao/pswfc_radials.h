#ifndef PSWFC_RADIALS_H_
#define PSWFC_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"
#include <map>
#include <vector>

class PswfcRadials : public RadialSet {
    public:
        PswfcRadials() {};
        PswfcRadials& operator=(const PswfcRadials& rhs);
        PswfcRadials* clone() const { return new PswfcRadials(*this); }
        ~PswfcRadials() {};
        /// @brief central function to build RadialCollection from ONCVPSP program generated pseudopotential file
        /// @param file file name of pseudopotential file
        /// @param itype atomic type, indiced in UnitCell class
        /// @param screening_coeff screening coefficient of pseudowavefunction
        /// @param conv_thr convergence threshold of norm of pseudowavefunction, see function cut_to_convergence for details
        /// @param ptr_log output file stream for logging
        /// @param rank MPI rank
        void build(const std::string& file = "", 
                   const int itype = 0,
                   const double screening_coeff = 0.1,
                   const double conv_thr = 1e-10,
                   std::ofstream* ptr_log = nullptr, 
                   const int rank = 0);

        /// @brief read ONCVPSP program generated pseudopotential file, and store the radial functions in RadialCollection
        /// @param ifs input file stream from orbital file
        /// @param screening_coeff screening coefficient of pseudowavefunction
        /// @param conv_thr convergence threshold of norm of pseudowavefunction, see function cut_to_convergence for details
        /// @param ptr_log output file stream for logging
        /// @param rank MPI rank
        void read_upf_pswfc(std::ifstream& ifs,               //!< input file stream from orbital file
                            const double screening_coeff,     //!< screening coefficient
                            const double conv_thr,            //!< convergence threshold
                            std::ofstream* ptr_log = nullptr, //!< output file stream for logging
                            const int rank = 0                //!< MPI rank
        );

        /// @brief returns the norm of the radial function
        /// @param rgrid radial grid
        /// @param rvalue radial function
        /// @return norm of the radial function
        double radial_norm(const std::vector<double> rgrid,
                           const std::vector<double> rvalue);
        /// @brief python-like startswith function
        /// @param word as it is
        /// @param pattern pattern to be matched
        /// @return true if word starts with pattern
        bool startswith(std::string word, std::string pattern);
        /// @brief read value from attributes in HTML-like format
        /// @param ifs input file stream
        /// @param word as it is
        /// @return value of the attribute
        std::string read_keyword_value(std::ifstream& ifs, std::string word);
        /// @brief steal string from quotes
        /// @param word as it is
        /// @return string between quotes
        std::string steal_from_quotes(std::string word);
        /// @brief steal string from quotes
        /// @param ifs input file stream
        /// @param word as it is
        /// @return string between quotes
        std::string steal_from_quotes(std::ifstream& ifs, std::string word);

        /// @brief cut radial function to convergence
        /// @param rgrid radial grid
        /// @param rvalue radial function
        /// @param conv_thr convergence of norm of radial function
        /// @return cutoff radius
        double cut_to_convergence(const std::vector<double>& rgrid, 
                                  std::vector<double>& rvalue, 
                                  const double& conv_thr);
        /// @brief smooth the radial function to avoid high frequency noise in FFT-spherical bessel transform
        /// @param rgrid radial grid
        /// @param rvalue radial function
        /// @param sigma sigma of the Gaussian kernel
        void smooth(std::vector<double>& rgrid,
                    std::vector<double>& rvalue,
                    const double sigma = 0.1);
        /// @brief call cut_to_convergence for each (l,zeta) corresponding orbital in std::map, then zero-padding to the maximal r, generate a grid
        /// @param pswfc_map a map of (l,zeta) corresponding orbital
        /// @return a vector of radial grid
        std::vector<double> pswfc_prepossess(std::map<std::pair<int, int>, std::vector<double>>& pswfc_map,
                                             const double conv_thr = 1e-6);
    private:

};
#endif // PSWFC_RADIALS_H_
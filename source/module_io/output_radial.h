#ifndef RADIAL_IO_H_
#define RADIAL_IO_H_

#include <fstream>
#include <vector>

namespace ModuleIO
{
    class OutputRadial
    {
        public:
            OutputRadial() {};
            ~OutputRadial() {};
            /// @brief will create a ofstream object and open the file
            /// @param filename name of the file to be opened
            void initialize(const std::string& filename);
            /// @brief configure basic information of the file
            /// @param symbol element symbol
            /// @param ecut energy cutoff of numerical atomic orbital
            /// @param rcut radius cutoff of numerical atomic orbital
            /// @param lmax maximum angular momentum
            /// @param l_nchi number of radial functions for each angular momentum
            /// @param ngrid number of grid points
            /// @param dr grid spacing
            void configure(const std::string& symbol = "H",
                        const double& ecut = 100.0,
                        const double& rcut = 10.0,
                        const int& lmax = 0, 
                        const int* l_nchi = nullptr,
                        const int& ngrid = 0,
                        const double& dr = 0.01
                        );
            /// @brief import radialset data in
            /// @param ngrid number of grid points
            /// @param rgrid radial grid
            /// @param chi zeta function
            void push(int ngrid = 0, 
                    const double* rgrid = nullptr, 
                    const double* chi = nullptr);
            /// @brief will close the file and delete the ofstream object
            void finalize();

            // getters
            std::string symbol() const {return symbol_;};
            double ecut() const {return ecut_;};
            double rcut() const {return rcut_;};
            int lmax() const {return lmax_;};
            std::vector<int> l_nchi() const {return l_nchi_;};
            int ngrid() const {return ngrid_;};
            double dr() const {return dr_;};
            int current_l() const {return current_l_;};
            int current_chi() const {return current_chi_;};

        private:
            std::ofstream file_to_;
            std::string symbol_ = "H";
            double ecut_ = 100.0;
            double rcut_ = 10.0;
            int lmax_ = 0;
            std::vector<int> l_nchi_;
            int ngrid_ = 0;
            double dr_ = 0.01;

            int current_l_ = 0;
            int current_chi_ = 0;

            void write_header();
    };
} // namespace ModuleIO
#endif // RADIAL_IO_H_
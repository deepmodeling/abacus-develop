#ifndef PSI_INITIALIZER_NAO_H
#define PSI_INITIALIZER_NAO_H
#include "psi_initializer.h"
#include "module_base/realarray.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital method
*/
template <typename T, typename Device>
class psi_initializer_nao : public psi_initializer<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        psi_initializer_nao() {this->set_method("nao");};
        ~psi_initializer_nao() {};

        virtual void proj_ao_onkG(int ik) override;

        #ifdef __MPI // MPI additional implementation
        /// @brief initialize the psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                Parallel_Kpoints*,                      //< parallel kpoints
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr,          //< nonlocal pseudopotential
                                const int& = 0) override;               //< MPI rank
        #else
        /// @brief serial version of initialize function, link psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr) override;//< nonlocal pseudopotential
        #endif

        virtual void read_external_orbs(std::string* orbital_files,
                                        const int& rank) override;
        void allocate_table() override;
        void tabulate() override;
        
        std::vector<std::string> external_orbs() const { return orbital_files_; }
        std::vector<std::vector<int>> n_rgrid() const { return n_rgrid_; }
        std::vector<int> n_rgrid(const int& itype) const { return n_rgrid_[itype]; }
        int n_rgrid(const int& itype, const int& ichi) const { return n_rgrid_[itype][ichi]; }
        std::vector<std::vector<std::vector<double>>> rvalue() const { return rvalue_; }
        std::vector<std::vector<double>> rvalue(const int& itype) const { return rvalue_[itype]; }
        std::vector<double> rvalue(const int& itype, const int& ichi) const { return rvalue_[itype][ichi]; }
        double rvalue(const int& itype, const int& ichi, const int& ir) const { return rvalue_[itype][ichi][ir]; }
        std::vector<std::vector<std::vector<double>>> rgrid() const { return rgrid_; }
        std::vector<std::vector<double>> rgrid(const int& itype) const { return rgrid_[itype]; }
        std::vector<double> rgrid(const int& itype, const int& ichi) const { return rgrid_[itype][ichi]; }
        double rgrid(const int& itype, const int& ichi, const int& ir) const { return rgrid_[itype][ichi][ir]; }

    private:
        std::vector<std::string> orbital_files_;
        ModuleBase::realArray ovlp_flzjlq_;
        /// @brief number of realspace grids per type per chi, [itype][ichi]
        std::vector<std::vector<int>> n_rgrid_;
        /// @brief data of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
        std::vector<std::vector<std::vector<double>>> rvalue_;
        /// @brief r of numerical atomic orbital per type per chi per position, [itype][ichi][ir]
        std::vector<std::vector<std::vector<double>>> rgrid_;
};
#endif
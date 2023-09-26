#pragma once
#include "Exx_LRI.h"
#include <memory>

class Local_Orbital_Charge;
class LCAO_Matrix;
class Charge_Mixing;
namespace elecstate
{
class ElecState;
}

template<typename Tdata>
class Exx_LRI_Interface
{
public:
    /// @brief  Constructor for Exx_LRI_Interface
    /// @param exx_lri
    Exx_LRI_Interface(std::shared_ptr<Exx_LRI<Tdata>> exx_lri) : exx_lri(exx_lri) {}
    Exx_LRI_Interface() = delete;

    void write_Hexxs(const std::string &file_name) const;
    void read_Hexxs(const std::string& file_name);
    
    using TAC = std::pair<int, std::array<int, 3>>;
    std::vector< std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& get_Hexxs() const { return this->exx_lri->Hexxs; }
    
    Tdata& get_Eexx() const { return this->exx_lri->Eexx; }

    // Processes in ESolver_KS_LCAO
    /// @brief in beforescf: set xc type, opt_orb, do DM mixing
    void exx_beforescf(const K_Vectors& kv, const Charge_Mixing& chgmix);

    /// @brief in eachiterinit:  do DM mixing and calculate Hexx when entering 2nd SCF
    void exx_eachiterinit(const Local_Orbital_Charge& loc, const Charge_Mixing& chgmix, const int& iter);

    /// @brief in hamilt2density: calculate Hexx and Eexx
    void exx_hamilt2density(elecstate::ElecState& elec, const Parallel_Orbitals& pv);

    /// @brief: in do_after_converge: add exx operators; do DM mixing if seperate loop
    bool exx_after_converge(
        hamilt::Hamilt<std::complex<double>>& hamilt,
        LCAO_Matrix& lm,
        const Local_Orbital_Charge& loc,
        const K_Vectors& kv,
        int& iter);
    
private:
    std::shared_ptr<Exx_LRI<Tdata>> exx_lri;
};
#include "Exx_LRI_interface.hpp"
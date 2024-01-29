#ifndef EXX_LRI_INTERFACE_H
#define EXX_LRI_INTERFACE_H

#include "Exx_LRI.h"
#include "module_ri/Mix_DMk_2D.h"
#include <memory>

class Local_Orbital_Charge;
class LCAO_Matrix;
class Charge_Mixing;
namespace elecstate
{
    class ElecState;
    template <typename TK, typename TR>
    class DensityMatrix;
}

template<typename T, typename Tdata>
class Exx_LRI_Interface
{
public:
    /// @brief  Constructor for Exx_LRI_Interface
    /// @param exx_ptr
    Exx_LRI_Interface(std::shared_ptr<Exx_LRI<Tdata>> exx_ptr) : exx_ptr(exx_ptr) {}
    Exx_LRI_Interface() = delete;

    void write_Hexxs(const std::string &file_name) const;
    void read_Hexxs(const std::string& file_name);
    
    using TAC = std::pair<int, std::array<int, 3>>;
    std::vector< std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& get_Hexxs() const { return this->exx_ptr->Hexxs; }
    
    Tdata& get_Eexx() const { return this->exx_ptr->Eexx; }

    // Processes in ESolver_KS_LCAO
    /// @brief in beforescf: set xc type, opt_orb, do DM mixing
    void exx_beforescf(const K_Vectors& kv, const Charge_Mixing& chgmix);

    /// @brief in eachiterinit:  do DM mixing and calculate Hexx when entering 2nd SCF
    void exx_eachiterinit(const elecstate::DensityMatrix<T, double>& dm/**< double should be Tdata if complex-PBE-DM is supported*/,
        const int& iter);

    /// @brief in hamilt2density: calculate Hexx and Eexx
    void exx_hamilt2density(elecstate::ElecState& elec, const Parallel_Orbitals& pv);

    /// @brief: in do_after_converge: add exx operators; do DM mixing if seperate loop
    bool exx_after_converge(
        hamilt::Hamilt<T>& hamilt,
        LCAO_Matrix& lm,
        const elecstate::DensityMatrix<T, double>& dm/**< double should be Tdata if complex-PBE-DM is supported*/,
        const K_Vectors& kv,
        int& iter);

private:
    std::shared_ptr<Exx_LRI<Tdata>> exx_ptr;
    Mix_DMk_2D mix_DMk_2D;
    int two_level_step = 0;
};

#include "Exx_LRI_interface.hpp"

#endif
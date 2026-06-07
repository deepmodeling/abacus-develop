#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H

#ifdef __EXX

#include "operator_lcao.h"
#include "source_cell/klist.h"
#include "source_hamilt/module_xc/exx_info.h"

#include <RI/global/Tensor.h>
#include <RI/ri/Cell_Nearest.h>

// Forward declaration to avoid circular include (Exx_LRI_interface.hpp includes op_exx_lcao.h)
template <typename T, typename Tdata>
class Exx_LRI_Interface;

namespace hamilt
{

#ifndef __OPEXXTEMPLATE
#define __OPEXXTEMPLATE

template <class T>
class OperatorEXX : public T
{
};

#endif
enum Add_Hexx_Type
{
    R,
    k
};
template <typename TK, typename TR>
class OperatorEXX<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
    using TAC = std::pair<int, std::array<int, 3>>;

  public:
    /// @brief Full-workflow Constructor  that takes Exx_LRI_Interface objects directly.
    /// Used in the main project (HamiltLCAO) for both scf and nscf.
    OperatorEXX<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                      hamilt::HContainer<TR>* hR_in,
                                      const UnitCell& ucell,
                                      const K_Vectors& kv_in,
                                      Exx_LRI_Interface<TK, double>* exd_in,
                                      Exx_LRI_Interface<TK, std::complex<double>>* exc_in,
                                      Add_Hexx_Type add_hexx_type_in = Add_Hexx_Type::R,
                                      const int istep_in = 0,
                                      const bool restart_in = false);

    /// @brief One-shot operator constructor, only for adding Hexxs, without exd/exc workflow
    /// Used in write_Vxc
    OperatorEXX<OperatorLCAO<TK, TR>>(
        HS_Matrix_K<TK>* hsk_in,
        hamilt::HContainer<TR>* hR_in,
        const UnitCell& ucell,
        const K_Vectors& kv_in,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in = nullptr,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in = nullptr,
        Add_Hexx_Type add_hexx_type_in = Add_Hexx_Type::R);

    virtual void contributeHk(int ik) override;
    virtual void contributeHR() override;

    template <typename Tdata>
    void cal_dH(const int ispin,
        std::array<std::vector<hamilt::HContainer<double>*>, 3>& dhR,
        const std::array<std::vector<std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>>, 3>& dHexxs);

  private:
    Add_Hexx_Type add_hexx_type = Add_Hexx_Type::R;
    int current_spin = 0;
    bool HR_fixed_done = false;
    bool initial_gga_done = false; // Taoni Bao add 2026-05-18, to fix RT-TDDFT EXX missing problem in the evolution

    /// @brief Non-owning pointers to the EXX interface objects.
    /// When set (via the interface-based constructor), Hexxd/Hexxc/two_level_step
    /// are sourced from these objects rather than stored as separate members.
    Exx_LRI_Interface<TK, double>* exd = nullptr;
    Exx_LRI_Interface<TK, std::complex<double>>* exc = nullptr;

    std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr;

    /// @brief if restart, read and save Hexx, and directly use it during the first outer loop.
    bool restart = false;

    const int istep = 0; // the ion step

    void add_loaded_Hexx(const int ik);

    const UnitCell& ucell;

    const K_Vectors& kv;

    // if k points has no shift, use cell_nearest to reduce the memory cost
    RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
    bool use_cell_nearest = true;

    /// @brief Hexxk for all k-points, only for the 1st scf loop ofrestart load
    std::vector<std::vector<double>> Hexxd_k_load;
    std::vector<std::vector<std::complex<double>>> Hexxc_k_load;
};

using TAC = std::pair<int, std::array<int, 3>>;

RI::Cell_Nearest<int, int, 3, double, 3> init_cell_nearest(const UnitCell& ucell, const std::array<int, 3>& Rs_period);

// allocate according to the read-in HexxR, used in nscf
template <typename Tdata, typename TR>
void reallocate_hcontainer(const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
                           HContainer<TR>* hR,
                           const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest = nullptr);

/// allocate according to BvK cells, used in scf
template <typename TR>
void reallocate_hcontainer(const int nat,
                           HContainer<TR>* hR,
                           const std::array<int, 3>& Rs_period,
                           const RI::Cell_Nearest<int, int, 3, double, 3>* const cell_nearest = nullptr);

} // namespace hamilt
#endif // __EXX
#include "op_exx_lcao.hpp"
#endif // OPEXXLCAO_H
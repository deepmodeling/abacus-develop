#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H

#ifdef __EXX

#include "operator_lcao.h"
#include "source_cell/klist.h"
#include "source_hamilt/module_xc/exx_info.h"

#include <RI/global/Tensor.h>
#include <RI/ri/Cell_Nearest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

// Forward declaration to avoid circular include (exx_lri_interface.hpp includes op_exx_lcao.h)
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
                                      const Exx_Info& exx_info,
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
        std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in,
        const Exx_Info* exx_info,
        Add_Hexx_Type add_hexx_type_in);

    // Retained for the one-shot Add_Hexx_Type::k path used by write_Vxc/RDMFT;
    // the main LCAO path uses contributeHR().
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

    /// @brief EXX info, passed from ESolver
    const Exx_Info* exx_info_ptr = nullptr;

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

template <typename Tdata>
inline std::pair<bool, std::array<int, 3>> infer_complete_Rs_period_from_Hexxs(
    const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
{
    std::array<int, 3> Rs_period = {0, 0, 0};
    if (Hexxs.empty())
    {
        return {false, Rs_period};
    }

    std::set<std::array<int, 3>> unique_Rs;
    std::array<int, 3> min_R = {
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max(),
        std::numeric_limits<int>::max()};
    std::array<int, 3> max_R = {
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min(),
        std::numeric_limits<int>::min()};

    for (const auto& Htmp1 : Hexxs[0])
    {
        for (const auto& Htmp2 : Htmp1.second)
        {
            const std::array<int, 3>& cell = Htmp2.first.second;
            if (unique_Rs.insert(cell).second)
            {
                for (int idim = 0; idim < 3; ++idim)
                {
                    min_R[idim] = std::min(min_R[idim], cell[idim]);
                    max_R[idim] = std::max(max_R[idim], cell[idim]);
                }
            }
        }
    }

    if (unique_Rs.empty())
    {
        return {false, Rs_period};
    }

    std::size_t expected_nR = 1;
    for (int idim = 0; idim < 3; ++idim)
    {
        Rs_period[idim] = max_R[idim] - min_R[idim] + 1;
        if (Rs_period[idim] <= 0)
        {
            return {false, {0, 0, 0}};
        }
        expected_nR *= static_cast<std::size_t>(Rs_period[idim]);
    }

    if (expected_nR != unique_Rs.size())
    {
        return {false, {0, 0, 0}};
    }

    return {true, Rs_period};
}

inline bool can_remap_wigner_seitz_for_nscf(const Add_Hexx_Type add_hexx_type,
                                             const bool period_from_kmesh,
                                             const bool zero_koffset)
{
    return add_hexx_type == Add_Hexx_Type::R && (!period_from_kmesh || zero_koffset);
}

struct WignerSeitzRemapStats
{
    std::size_t input_blocks = 0;
    std::size_t output_blocks = 0;
    std::size_t remapped_blocks = 0;
    std::size_t split_blocks = 0;
    std::size_t max_images = 0;
};

inline std::vector<std::array<int, 3>> find_nearest_bvk_cells(const UnitCell& ucell,
                                                               const int iat0,
                                                               const int iat1,
                                                               const std::array<int, 3>& R,
                                                               const std::array<int, 3>& period)
{
    // Boundary-equivalent BvK images must all carry the same matrix weight.
    double min_distance_sq = std::numeric_limits<double>::max();
    std::vector<std::array<int, 3>> nearest_cells;

    for (int ix = -1; ix <= 1; ++ix)
    {
        for (int iy = -1; iy <= 1; ++iy)
        {
            for (int iz = -1; iz <= 1; ++iz)
            {
                const std::array<int, 3> candidate = {
                    R[0] + ix * period[0], R[1] + iy * period[1], R[2] + iz * period[2]};
                const ModuleBase::Vector3<int> candidate_vector(candidate[0], candidate[1], candidate[2]);
                const double distance_sq = ucell.cal_dtau(iat0, iat1, candidate_vector).norm2();
                const double tolerance = 1.0e-10
                                         * std::max(1.0,
                                                    std::max(std::abs(distance_sq), std::abs(min_distance_sq)));

                if (distance_sq + tolerance < min_distance_sq)
                {
                    min_distance_sq = distance_sq;
                    nearest_cells.clear();
                    nearest_cells.push_back(candidate);
                }
                else if (std::abs(distance_sq - min_distance_sq) <= tolerance)
                {
                    nearest_cells.push_back(candidate);
                }
            }
        }
    }
    return nearest_cells;
}

template <typename Tdata>
inline WignerSeitzRemapStats remap_Hexxs_wigner_seitz(
    const UnitCell& ucell,
    const std::array<int, 3>& period,
    std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs)
{
    WignerSeitzRemapStats stats;
    std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>> remapped(Hexxs.size());

    for (std::size_t is = 0; is < Hexxs.size(); ++is)
    {
        for (const auto& Htmp1 : Hexxs[is])
        {
            const int iat0 = Htmp1.first;
            for (const auto& Htmp2 : Htmp1.second)
            {
                ++stats.input_blocks;
                const int iat1 = Htmp2.first.first;
                const std::array<int, 3>& R = Htmp2.first.second;
                const auto nearest_cells = find_nearest_bvk_cells(ucell, iat0, iat1, R, period);
                if (nearest_cells.empty())
                {
                    throw std::runtime_error("No nearest BvK image found for HexxR block");
                }

                stats.max_images = std::max(stats.max_images, nearest_cells.size());
                if (nearest_cells.size() != 1 || nearest_cells.front() != R)
                {
                    ++stats.remapped_blocks;
                }
                if (nearest_cells.size() > 1)
                {
                    ++stats.split_blocks;
                }

                const Tdata weight = static_cast<Tdata>(1.0 / static_cast<double>(nearest_cells.size()));
                auto& target = remapped[is][iat0];
                for (const auto& R_bvk : nearest_cells)
                {
                    const TAC key = {iat1, R_bvk};
                    auto it = target.find(key);
                    if (it == target.end())
                    {
                        target.emplace(key, weight * Htmp2.second);
                    }
                    else
                    {
                        it->second += weight * Htmp2.second;
                    }
                }
            }
        }
    }

    Hexxs.swap(remapped);
    for (const auto& Hspin : Hexxs)
    {
        for (const auto& Htmp1 : Hspin)
        {
            stats.output_blocks += Htmp1.second.size();
        }
    }
    return stats;
}

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
#endif // OPEXXLCAO_H

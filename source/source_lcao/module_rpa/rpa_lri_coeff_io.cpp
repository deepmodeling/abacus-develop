#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_coeff_v1.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_basis/module_ao/elem_basis_idx_orb.h"
#include "source_estate/elecstate_lcao.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_lcao/module_ri/module_exx_symmetry/symm_rotation.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs(const UnitCell& ucell,
                               std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in,
                               std::string filename)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs");
    ModuleBase::timer::start("RPA_LRI", "out_Cs");

    std::stringstream ss;
    ss << filename << (this->runtime.rank + 1) << ".txt";
    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + ss.str().c_str(), std::ios::out);
    ofs << ucell.nat << "    " << 0 << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Cs_in)
    {
        size_t I = Ip.first;
        size_t i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            size_t J = JPp.first.first;
            auto R = JPp.first.second;
            auto& tmp_Cs = JPp.second;
            size_t j_num = ucell.atoms[ucell.iat2it[J]].nw;

            ofs << I + 1 << "   " << J + 1 << "   " << R[0] << "   " << R[1] << "   " << R[2] << "   " << i_num
                << std::endl;
            ofs << j_num << "   " << tmp_Cs.shape[0] << std::endl;
            for (int i = 0; i != i_num; i++)
            {
                for (int j = 0; j != j_num; j++)
                {
                    for (int mu = 0; mu != tmp_Cs.shape[0]; mu++)
                    {
                        ofs << std::setw(30) << tmp_Cs(mu, i, j) << std::endl;
                    }
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_Cs");
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs_v1(const UnitCell& ucell,
                                  std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in,
                                  std::string filename)
try
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs_v1");
    ModuleBase::timer::start("RPA_LRI", "out_Cs_v1");

    std::vector<RpaLriDetail::LriCoefV1Block> records;
    records.reserve(Cs_in.size());
    for (auto& Ip: Cs_in)
    {
        const int I = static_cast<int>(Ip.first);
        const int i_num = ucell.atoms[ucell.iat2it[I]].nw;
        for (auto& JPp: Ip.second)
        {
            const int J = static_cast<int>(JPp.first.first);
            auto& tmp_Cs = JPp.second;
            if (tmp_Cs.shape.size() != 3 || tmp_Cs.shape[0] <= 0 || tmp_Cs.shape[1] <= 0 || tmp_Cs.shape[2] <= 0)
            {
                continue;
            }
            const int j_num = ucell.atoms[ucell.iat2it[J]].nw;
            if (static_cast<int>(tmp_Cs.shape[1]) != i_num || static_cast<int>(tmp_Cs.shape[2]) != j_num)
            {
                throw std::runtime_error("LibRPA v1 Cs output encountered an inconsistent tensor shape.");
            }
            RpaLriDetail::LriCoefV1Block record;
            record.ia1 = I + 1;
            record.ia2 = J + 1;
            for (int axis = 0; axis < 3; ++axis)
            {
                record.r[axis] = JPp.first.second[axis];
            }
            record.tensor = &tmp_Cs;
            record.nw1 = i_num;
            record.nw2 = j_num;
            record.naux = static_cast<int>(tmp_Cs.shape[0]);
            for (int i = 0; i != record.nw1; ++i)
            {
                for (int j = 0; j != record.nw2; ++j)
                {
                    for (int mu = 0; mu != record.naux; ++mu)
                    {
                        record.max_abs
                            = std::max(record.max_abs, std::abs(RpaLriDetail::real_as_double(tmp_Cs(mu, i, j))));
                    }
                }
            }
            records.push_back(record);
        }
    }

    std::stringstream ss;
    ss << filename << this->runtime.rank << ".dat";
    const std::string out_name = this->runtime.input.rpa_outdir + ss.str();
    RpaLriDetail::write_lri_coef_v1(out_name, static_cast<std::int32_t>(ucell.nat), records);
    ModuleBase::timer::end("RPA_LRI", "out_Cs_v1");
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_Cs_v1");
}

template void RPA_LRI<double, double>::out_Cs(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Cs_in,
    std::string filename);
template void RPA_LRI<std::complex<double>, double>::out_Cs(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Cs_in,
    std::string filename);
template void RPA_LRI<double, double>::out_Cs_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Cs_in,
    std::string filename);
template void RPA_LRI<std::complex<double>, double>::out_Cs_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& Cs_in,
    std::string filename);

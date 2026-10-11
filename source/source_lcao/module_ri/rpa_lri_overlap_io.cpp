#include "exx_lri.h"
#include "rpa_lri.h"
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
void RPA_LRI<T, Tdata>::out_abfs_overlap(const UnitCell& ucell,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                         std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                                         std::string filename,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                                         const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs)
{
    ModuleBase::TITLE("RPA_LRI", "out_abfs_overlap");
    ModuleBase::timer::start("RPA_LRI", "out_abfs_overlap");
    const double threshold = 1e-15;
    const auto format = std::scientific;
    int prec = 15;

    int all_mu_s = 0;
    int all_mu = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    std::vector<int> mu_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        mu_shift[I] = all_mu;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
        all_mu += index_abfs[ucell.iat2it[I]].count_size;
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (this->runtime.rank + 1) << ".txt";

    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;

    // Fourier of ss(R->k), s(R->k)
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_ss;
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_s;
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                RI::Tensor<std::complex<double>> tmp_olp_ss
                    = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                RI::Tensor<std::complex<double>> tmp_olp_s
                    = RI::Global_Func::convert<std::complex<double>>(overlap_abfs_abf[I][{J, R}]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_ss.shape[0], tmp_olp_ss.shape[1]});
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_s.shape[0], tmp_olp_s.shape[1]});
                }
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));

                olp_q_ss[I][{J, q}] = olp_q_ss[I][{J, q}] + tmp_olp_ss * kphase;
                olp_q_s[I][{J, q}] = olp_q_s[I][{J, q}] + tmp_olp_s * kphase;
            }
        }
    }
    // for multi-mpi
    for (int I = 0; I != ucell.nat; I++)
    {
        for (int J = 0; J != ucell.nat; J++)
        {
            for (int ik = 0; ik != nks_tot; ik++)
            {
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs_s[ucell.iat2it[J]].count_size;
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                if (olp_q_s[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs[ucell.iat2it[J]].count_size;
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                for (int ir = 0; ir < olp_q_ss[I][{J, q}].shape[0]; ir++)
                {
                    for (int ic = 0; ic < olp_q_ss[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_ss[I][{J, q}](ir, ic));
                    }
                    for (int ic = 0; ic < olp_q_s[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_s[I][{J, q}](ir, ic));
                    }
                }
            }
        }
    }

    // out_ri_tensor("olp_ss.txt", olp_q_ss, 0.);
    // Inverse of overlap(q)
    inverse_olp(ucell, olp_q_ss, index_abfs_s);
    // out_ri_tensor("olp_ss_inv.txt", olp_q_ss, 0.);
    // out_ri_tensor("olp_s.txt", olp_q_s, 0.);
    for (auto& Ip: overlap_abfs_abf)
    {
        auto I = Ip.first;
        size_t mu_num_s = index_abfs_s[ucell.iat2it[I]].count_size;
        size_t mu_num = index_abfs[ucell.iat2it[I]].count_size;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> sinvS;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                if (sinvS[J].empty())
                {
                    sinvS[J] = RI::Tensor<std::complex<double>>(
                        {overlap_abfs_abfs[I][{J, R}].shape[0], overlap_abfs_abf[I][{J, R}].shape[1]});
                }
            }
            for (const auto& pair: sinvS)
            {
                auto J = pair.first;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                for (int K = 0; K != ucell.nat; K++)
                {
                    sinvS[J] += olp_q_ss.at(I).at({K, q}) * olp_q_s.at(K).at({J, q});
                }
            }
            for (auto& iJU: sinvS)
            {
                auto iJ = iJU.first;
                auto& vq_J = iJU.second;
                size_t nu_num = index_abfs[ucell.iat2it[iJ]].count_size;
                ofs << all_mu_s << "   " << all_mu << "   " << mu_s_shift[I] + 1 << "   " << mu_s_shift[I] + mu_num_s
                    << "  " << mu_shift[iJ] + 1 << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * this->runtime.input.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    // ofs << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].real()
                    //     << std::setw(25) << std::fixed << std::setprecision(15) << (*vq_J.data)[i].imag() <<
                    //     std::endl;
                    // if (fabs((*vq_J.data)[i].real()) > threshold || fabs((*vq_J.data)[i].imag()) > threshold)
                    ofs << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].real() << " "
                        << std::showpoint << format << std::setprecision(prec) << (*vq_J.data)[i].imag() << "\n";
                    // else
                    //     ofs << std::showpoint << format << std::setprecision(prec) << 0.0 << " " << std::showpoint
                    //         << format << std::setprecision(prec) << 0.0 << "\n";
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_abfs_overlap");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap_v1(const UnitCell& ucell,
                                            std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
                                            std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abf,
                                            std::string filename,
                                            const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
                                            const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs)
try
{
    ModuleBase::TITLE("RPA_LRI", "out_abfs_overlap_v1");
    ModuleBase::timer::start("RPA_LRI", "out_abfs_overlap_v1");

    struct SinvSRecord
    {
        std::int32_t iq = 0;
        std::int32_t nrow_total = 0;
        std::int32_t ncol_total = 0;
        std::int32_t begin_row = 0;
        std::int32_t end_row = 0;
        std::int32_t begin_col = 0;
        std::int32_t end_col = 0;
        double q_weight = 0.0;
        std::int64_t offset = 0;
        std::vector<std::complex<double>> payload;
    };

    int all_mu_s = 0;
    int all_mu = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    std::vector<int> mu_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        mu_shift[I] = all_mu;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
        all_mu += index_abfs[ucell.iat2it[I]].count_size;
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();

    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_ss;
    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>> olp_q_s;
    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                RI::Tensor<std::complex<double>> tmp_olp_ss
                    = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                RI::Tensor<std::complex<double>> tmp_olp_s
                    = RI::Global_Func::convert<std::complex<double>>(overlap_abfs_abf[I][{J, R}]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_ss.shape[0], tmp_olp_ss.shape[1]});
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({tmp_olp_s.shape[0], tmp_olp_s.shape[1]});
                }
                const double arg
                    = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec)) * ModuleBase::TWO_PI;
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                olp_q_ss[I][{J, q}] = olp_q_ss[I][{J, q}] + tmp_olp_ss * kphase;
                olp_q_s[I][{J, q}] = olp_q_s[I][{J, q}] + tmp_olp_s * kphase;
            }
        }
    }

    for (int I = 0; I != ucell.nat; I++)
    {
        for (int J = 0; J != ucell.nat; J++)
        {
            for (int ik = 0; ik != nks_tot; ik++)
            {
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                if (olp_q_ss[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs_s[ucell.iat2it[J]].count_size;
                    olp_q_ss[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                if (olp_q_s[I][{J, q}].empty())
                {
                    auto mu = index_abfs_s[ucell.iat2it[I]].count_size;
                    auto nu = index_abfs[ucell.iat2it[J]].count_size;
                    olp_q_s[I][{J, q}] = RI::Tensor<std::complex<double>>({mu, nu});
                }
                for (int ir = 0; ir < olp_q_ss[I][{J, q}].shape[0]; ir++)
                {
                    for (int ic = 0; ic < olp_q_ss[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_ss[I][{J, q}](ir, ic));
                    }
                    for (int ic = 0; ic < olp_q_s[I][{J, q}].shape[1]; ic++)
                    {
                        Parallel_Reduce::reduce_all<std::complex<double>>(olp_q_s[I][{J, q}](ir, ic));
                    }
                }
            }
        }
    }

    inverse_olp(ucell, olp_q_ss, index_abfs_s);

    std::vector<SinvSRecord> records;
    for (auto& Ip: overlap_abfs_abf)
    {
        auto I = Ip.first;
        size_t mu_num_s = index_abfs_s[ucell.iat2it[I]].count_size;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> sinvS;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto R = JPp.first.second;
                if (sinvS[J].empty())
                {
                    sinvS[J] = RI::Tensor<std::complex<double>>(
                        {overlap_abfs_abfs[I][{J, R}].shape[0], overlap_abfs_abf[I][{J, R}].shape[1]});
                }
            }
            for (const auto& pair: sinvS)
            {
                auto J = pair.first;
                auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
                for (int K = 0; K != ucell.nat; K++)
                {
                    sinvS[J] += olp_q_ss.at(I).at({K, q}) * olp_q_s.at(K).at({J, q});
                }
            }
            for (auto& iJU: sinvS)
            {
                auto iJ = iJU.first;
                auto& vq_J = iJU.second;
                size_t nu_num = index_abfs[ucell.iat2it[iJ]].count_size;
                SinvSRecord record;
                record.iq = static_cast<std::int32_t>(ik + 1);
                record.nrow_total = static_cast<std::int32_t>(all_mu_s);
                record.ncol_total = static_cast<std::int32_t>(all_mu);
                record.begin_row = static_cast<std::int32_t>(mu_s_shift[I] + 1);
                record.end_row = static_cast<std::int32_t>(mu_s_shift[I] + mu_num_s);
                record.begin_col = static_cast<std::int32_t>(mu_shift[iJ] + 1);
                record.end_col = static_cast<std::int32_t>(mu_shift[iJ] + nu_num);
                record.q_weight = p_kv->wk[ik] / 2.0 * this->runtime.input.nspin;
                record.payload.reserve(static_cast<std::size_t>(vq_J.shape[0])
                                       * static_cast<std::size_t>(vq_J.shape[1]));
                for (int i = 0; i != vq_J.shape[0]; ++i)
                {
                    for (int j = 0; j != vq_J.shape[1]; ++j)
                    {
                        if (i >= static_cast<int>(mu_num_s) || j >= static_cast<int>(nu_num))
                        {
                            throw std::runtime_error("LibRPA v1 shrink_sinvS encountered an inconsistent block shape.");
                        }
                        record.payload.push_back(vq_J(i, j));
                    }
                }
                records.push_back(std::move(record));
            }
        }
    }

    const std::int64_t record_bytes = 7 * static_cast<std::int64_t>(sizeof(std::int32_t))
                                      + static_cast<std::int64_t>(sizeof(double))
                                      + static_cast<std::int64_t>(sizeof(std::int64_t));
    std::int64_t offset = 2 * static_cast<std::int64_t>(sizeof(std::int32_t))
                          + static_cast<std::int64_t>(records.size()) * record_bytes;
    for (auto& record: records)
    {
        record.offset = offset;
        offset += static_cast<std::int64_t>(record.payload.size() * sizeof(std::complex<double>));
    }

    std::stringstream ss;
    ss << filename << this->runtime.rank << ".dat";
    const std::string out_name = this->runtime.input.rpa_outdir + ss.str();
    std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + out_name);
    }

    const std::int32_t marker = RpaLriDetail::LIBRPA_SHRINK_SINVS_V1_MARKER;
    const std::int32_t nrecords = static_cast<std::int32_t>(records.size());
    RpaLriDetail::write_scalar(ofs, marker, out_name);
    RpaLriDetail::write_scalar(ofs, nrecords, out_name);
    for (const auto& record: records)
    {
        RpaLriDetail::write_scalar(ofs, record.iq, out_name);
        RpaLriDetail::write_scalar(ofs, record.nrow_total, out_name);
        RpaLriDetail::write_scalar(ofs, record.ncol_total, out_name);
        RpaLriDetail::write_scalar(ofs, record.begin_row, out_name);
        RpaLriDetail::write_scalar(ofs, record.end_row, out_name);
        RpaLriDetail::write_scalar(ofs, record.begin_col, out_name);
        RpaLriDetail::write_scalar(ofs, record.end_col, out_name);
        RpaLriDetail::write_scalar(ofs, record.q_weight, out_name);
        RpaLriDetail::write_scalar(ofs, record.offset, out_name);
    }
    for (const auto& record: records)
    {
        RpaLriDetail::checked_write(ofs,
                                    record.payload.data(),
                                    record.payload.size() * sizeof(std::complex<double>),
                                    out_name);
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_abfs_overlap_v1");
}
catch (const std::exception& error)
{
    RpaLriDetail::abort_output(this->mpi_comm, error.what());
}
catch (...)
{
    RpaLriDetail::abort_output(this->mpi_comm, "Unknown exception in out_abfs_overlap_v1");
}

template void RPA_LRI<double, double>::out_abfs_overlap(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abf,
    std::string filename,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
template void RPA_LRI<std::complex<double>, double>::out_abfs_overlap(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abf,
    std::string filename,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
template void RPA_LRI<double, double>::out_abfs_overlap_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abf,
    std::string filename,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);
template void RPA_LRI<std::complex<double>, double>::out_abfs_overlap_v1(
    const UnitCell& ucell,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abf,
    std::string filename,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs);

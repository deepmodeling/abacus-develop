#include "exx_lri.h"
#include "rpa_lri.h"
#include "rpa_lri_detail.h"
#include "source_base/global_function.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"
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
                // Fourier transform the real-space overlap block with
                // exp(i 2*pi k.R), using the lattice-vector convention.
                const double arg = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI;
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
                    // S^{-1}S(I,J;q) = sum_K S_ss(I,K;q) S_s(K,J;q).
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
                // Fourier transform the real-space overlap block with
                // exp(i 2*pi k.R), using the lattice-vector convention.
                const double arg
                    = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec)) * ModuleBase::TWO_PI;
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

    if (this->runtime.input.out_abf_overlap)
    {
        out_abfs_overlap_raw_v1(ucell, overlap_abfs_abfs, index_abfs_s);
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
                    // S^{-1}S(I,J;q) = sum_K S_ss(I,K;q) S_s(K,J;q).
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


template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_abfs_overlap_raw_v1(
    const UnitCell& ucell,
    const std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    if (!this->runtime.input.rpa || this->runtime.input.out_librpa_ver != 1
        || this->info.shrink_abfs_pca_thr < 0.0)
    {
        throw std::runtime_error("raw active-ABF overlap writer requires rpa=true, "
                                 "out_librpa_ver=1, and shrink ABFs.");
    }

    const int natom = ucell.nat;
    if (natom <= 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer found no atoms.");
    }
    std::vector<int> atom_naux(static_cast<std::size_t>(natom), 0);
    std::vector<int> atom_shift(static_cast<std::size_t>(natom), 0);
    int naux = 0;
    for (int I = 0; I < natom; ++I)
    {
        const int count = index_abfs_s[ucell.iat2it[I]].count_size;
        if (count <= 0 || naux > std::numeric_limits<int>::max() - count)
        {
            throw std::runtime_error("raw active-ABF overlap writer found an invalid active basis layout.");
        }
        atom_shift[static_cast<std::size_t>(I)] = naux;
        atom_naux[static_cast<std::size_t>(I)] = count;
        naux += count;
    }

    const int nks_tot = this->runtime.input.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    if (nks_tot <= 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer found no q points.");
    }
    const double hermitian_tol = 1e-10;
    const double rank_tol = 1e-12;
    const int natom_metadata = natom;
    int natom_min = natom_metadata;
    int natom_max = natom_metadata;
    Parallel_Reduce::reduce_min(mpi_comm, natom_min);
    Parallel_Reduce::reduce_max(mpi_comm, natom_max);
    if (natom_min != natom_metadata || natom_max != natom_metadata)
    {
        throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent atom metadata.");
    }
    for (const int count: atom_naux)
    {
        int min_count = count;
        int max_count = count;
        Parallel_Reduce::reduce_min(mpi_comm, min_count);
        Parallel_Reduce::reduce_max(mpi_comm, max_count);
        if (min_count != count || max_count != count)
        {
            throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent active basis metadata.");
        }
    }

    // comm_map2_first may concentrate all R blocks for one (I,J) on one rank.
    // Check only actual post-communication keys; no complete R coverage is assumed.
    int local_invalid_metadata = 0;
    std::string local_metadata_error;
    std::vector<std::array<int, 5>> local_keys;
    std::vector<int> local_pair_seen(static_cast<std::size_t>(natom) * static_cast<std::size_t>(natom), 0);
    for (const auto& Ip: overlap_abfs_abfs)
    {
        const int I = Ip.first;
        if (I < 0 || I >= natom)
        {
            local_invalid_metadata = 1;
            if (local_metadata_error.empty())
            {
                local_metadata_error = "raw active-ABF overlap writer found an invalid row atom.";
            }
            continue;
        }
        for (const auto& JPp: Ip.second)
        {
            const int J = JPp.first.first;
            const auto R = JPp.first.second;
            if (J < 0 || J >= natom)
            {
                local_invalid_metadata = 1;
                if (local_metadata_error.empty())
                {
                    local_metadata_error = "raw active-ABF overlap writer found invalid block (I,J,R)=("
                                            + std::to_string(I) + "," + std::to_string(J) + ","
                                            + std::to_string(R[0]) + "," + std::to_string(R[1]) + ","
                                            + std::to_string(R[2]) + ").";
                }
                continue;
            }
            const auto& tensor = JPp.second;
            if (tensor.shape.size() != 2
                || tensor.shape[0] != atom_naux[static_cast<std::size_t>(I)]
                || tensor.shape[1] != atom_naux[static_cast<std::size_t>(J)])
            {
                local_invalid_metadata = 1;
                if (local_metadata_error.empty())
                {
                    local_metadata_error = "raw active-ABF overlap writer found inconsistent metadata for block (I,J,R)=("
                                            + std::to_string(I) + "," + std::to_string(J) + ","
                                            + std::to_string(R[0]) + "," + std::to_string(R[1]) + ","
                                            + std::to_string(R[2]) + ").";
                }
                continue;
            }
            local_keys.push_back({I, J, R[0], R[1], R[2]});
            local_pair_seen[static_cast<std::size_t>(I) * static_cast<std::size_t>(natom)
                            + static_cast<std::size_t>(J)] = 1;
        }
    }
    int global_invalid_metadata = local_invalid_metadata;
    Parallel_Reduce::reduce_max(mpi_comm, global_invalid_metadata);
    if (global_invalid_metadata != 0)
    {
        throw std::runtime_error(local_metadata_error.empty()
                                     ? "raw active-ABF overlap writer found invalid block metadata on another rank."
                                     : local_metadata_error);
    }
    for (int& seen: local_pair_seen)
    {
        Parallel_Reduce::reduce_max(mpi_comm, seen);
    }
    for (std::size_t pair_index = 0; pair_index < local_pair_seen.size(); ++pair_index)
    {
        if (local_pair_seen[pair_index] == 0)
        {
            throw std::runtime_error("raw active-ABF overlap writer found a missing atom pair.");
        }
    }

    // Sorting is also a defensive local-map duplicate check; std::map normally
    // makes such a duplicate impossible. The gathered check below rejects both
    // local and cross-rank repeats and reports the same concrete key on all ranks.
    std::sort(local_keys.begin(), local_keys.end());
    const bool local_duplicate = std::adjacent_find(local_keys.begin(), local_keys.end()) != local_keys.end();
    const std::size_t key_width = 5;
    const bool local_count_overflow
        = local_keys.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / key_width;
    int global_duplicate_flag = local_count_overflow ? 1 : 0;
    const int mpi_size = Parallel_Common::nproc(mpi_comm);
    Parallel_Reduce::reduce_max(mpi_comm, global_duplicate_flag);
    if (global_duplicate_flag != 0)
    {
        throw std::runtime_error("raw active-ABF overlap writer cannot represent MPI key count.");
    }
    int local_int_count = static_cast<int>(local_keys.size() * key_width);
    std::vector<int> recv_counts(static_cast<std::size_t>(mpi_size), 0);
    Parallel_Reduce::gather_int_all(mpi_comm, local_int_count, recv_counts.data());
    std::vector<int> displacements(static_cast<std::size_t>(mpi_size), 0);
    int total_int_count = 0;
    for (int rank = 0; rank < mpi_size; ++rank)
    {
        if (recv_counts[static_cast<std::size_t>(rank)] < 0
            || recv_counts[static_cast<std::size_t>(rank)] > std::numeric_limits<int>::max() - total_int_count)
        {
            throw std::runtime_error("raw active-ABF overlap writer cannot represent gathered MPI key counts.");
        }
        displacements[static_cast<std::size_t>(rank)] = total_int_count;
        total_int_count += recv_counts[static_cast<std::size_t>(rank)];
    }
    std::vector<int> local_packed;
    local_packed.reserve(static_cast<std::size_t>(local_int_count));
    for (const auto& key: local_keys)
    {
        local_packed.insert(local_packed.end(), key.begin(), key.end());
    }
    std::vector<int> gathered(static_cast<std::size_t>(total_int_count));
    Parallel_Reduce::gatherv_int_all(mpi_comm, local_packed.data(), local_int_count, gathered.data(),
                                     recv_counts.data(), displacements.data());
    std::vector<std::array<int, 5>> gathered_keys(static_cast<std::size_t>(total_int_count) / key_width);
    for (std::size_t index = 0; index < gathered_keys.size(); ++index)
    {
        std::copy_n(gathered.begin() + index * key_width, key_width, gathered_keys[index].begin());
    }
    std::sort(gathered_keys.begin(), gathered_keys.end());
    const auto duplicate = std::adjacent_find(gathered_keys.begin(), gathered_keys.end());
    if (local_duplicate || duplicate != gathered_keys.end())
    {
        throw std::runtime_error("raw active-ABF overlap writer found post-communication duplicate contributor for block (I,J,R)=("
                                 + std::to_string((*duplicate)[0]) + "," + std::to_string((*duplicate)[1]) + ","
                                 + std::to_string((*duplicate)[2]) + "," + std::to_string((*duplicate)[3]) + ","
                                 + std::to_string((*duplicate)[4]) + ").");
    }

    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const auto q = RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]);
        const double q_weight = p_kv->wk[ik] / 2.0 * this->runtime.input.nspin;
        if (!std::isfinite(q_weight)
            || !std::isfinite(q[0]) || !std::isfinite(q[1]) || !std::isfinite(q[2]))
        {
            throw std::runtime_error("raw active-ABF overlap writer found non-finite q metadata.");
        }
        for (const double coordinate: q)
        {
            double min_coordinate = coordinate;
            double max_coordinate = coordinate;
            Parallel_Reduce::reduce_min(mpi_comm, min_coordinate);
            Parallel_Reduce::reduce_max(mpi_comm, max_coordinate);
            if (max_coordinate - min_coordinate > rank_tol * std::max(1.0, std::abs(coordinate)))
            {
                throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent q coordinates.");
            }
        }
        std::vector<std::complex<double>> overlap(
            static_cast<std::size_t>(naux) * static_cast<std::size_t>(naux), std::complex<double>(0.0, 0.0));
        for (const auto& Ip: overlap_abfs_abfs)
        {
            const int I = Ip.first;
            for (const auto& JPp: Ip.second)
            {
                const int J = JPp.first.first;
                const auto& tensor = JPp.second;
                const auto R = JPp.first.second;
                const double arg = (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                    * ModuleBase::TWO_PI;
                const std::complex<double> phase(std::cos(arg), std::sin(arg));
                for (int ir = 0; ir < atom_naux[static_cast<std::size_t>(I)]; ++ir)
                {
                    for (int ic = 0; ic < atom_naux[static_cast<std::size_t>(J)]; ++ic)
                    {
                        const std::size_t row = static_cast<std::size_t>(atom_shift[static_cast<std::size_t>(I)] + ir);
                        const std::size_t col = static_cast<std::size_t>(atom_shift[static_cast<std::size_t>(J)] + ic);
                        overlap[row * static_cast<std::size_t>(naux) + col]
                            += static_cast<std::complex<double>>(tensor(ir, ic)) * phase;
                    }
                }
            }
        }
        for (std::complex<double>& value: overlap)
        {
            Parallel_Reduce::reduce_all<std::complex<double>>(value);
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            {
                throw std::runtime_error("raw active-ABF overlap writer found a non-finite overlap value.");
            }
        }
        for (int i = 0; i < naux; ++i)
        {
            for (int j = 0; j < naux; ++j)
            {
                const std::complex<double> difference
                    = overlap[static_cast<std::size_t>(i) * naux + j]
                    - std::conj(overlap[static_cast<std::size_t>(j) * naux + i]);
                const double scale = std::max(
                    1.0,
                    std::max(std::abs(overlap[static_cast<std::size_t>(i) * naux + j]),
                             std::abs(overlap[static_cast<std::size_t>(j) * naux + i])));
                if (std::abs(difference) > hermitian_tol * scale)
                {
                    throw std::runtime_error("raw active-ABF overlap writer found a non-Hermitian S(q).");
                }
            }
        }

        double min_weight = q_weight;
        double max_weight = q_weight;
        Parallel_Reduce::reduce_min(mpi_comm, min_weight);
        Parallel_Reduce::reduce_max(mpi_comm, max_weight);
        if (max_weight - min_weight > rank_tol * std::max(1.0, std::abs(q_weight)))
        {
            throw std::runtime_error("raw active-ABF overlap writer found rank-inconsistent q weight.");
        }

        if (this->runtime.rank != 0)
        {
            continue;
        }
        const std::string out_name = this->runtime.input.rpa_outdir + "v1_abf_overlap_active_iq_" + std::to_string(ik + 1) + ".dat";
        const std::string tmp_name = out_name + ".tmp";
        std::ofstream ofs(tmp_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
        if (!ofs.good())
        {
            throw std::runtime_error("Failed to open " + tmp_name);
        }
        const std::int32_t marker = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_MARKER;
        const std::int32_t version = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_VERSION;
        const std::int32_t iq = ik + 1;
        const std::int32_t kind = RpaLriDetail::LIBRPA_ABF_OVERLAP_V1_KIND_ACTIVE;
        const std::int32_t naux_i32 = naux;
        const std::int32_t natom_i32 = natom;
        RpaLriDetail::write_scalar(ofs, marker, tmp_name);
        RpaLriDetail::write_scalar(ofs, version, tmp_name);
        RpaLriDetail::write_scalar(ofs, iq, tmp_name);
        RpaLriDetail::write_scalar(ofs, kind, tmp_name);
        RpaLriDetail::write_scalar(ofs, naux_i32, tmp_name);
        RpaLriDetail::write_scalar(ofs, natom_i32, tmp_name);
        RpaLriDetail::write_scalar(ofs, q_weight, tmp_name);
        for (const double coordinate: q)
        {
            RpaLriDetail::write_scalar(ofs, coordinate, tmp_name);
        }
        for (const int count: atom_naux)
        {
            const std::int32_t count_i32 = count;
            RpaLriDetail::write_scalar(ofs, count_i32, tmp_name);
        }
        RpaLriDetail::checked_write(ofs, overlap.data(), overlap.size() * sizeof(std::complex<double>), tmp_name);
        ofs.close();
        if (!ofs.good() || std::rename(tmp_name.c_str(), out_name.c_str()) != 0)
        {
            throw std::runtime_error("Failed to finalize " + out_name);
        }
    }
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
template void RPA_LRI<double, double>::out_abfs_overlap_raw_v1(
    const UnitCell& ucell,
    const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);
template void RPA_LRI<std::complex<double>, double>::out_abfs_overlap_raw_v1(
    const UnitCell& ucell,
    const std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>& overlap_abfs_abfs,
    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s);

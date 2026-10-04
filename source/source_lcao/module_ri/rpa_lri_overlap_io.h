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
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI;
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
                record.payload.reserve(static_cast<std::size_t>(vq_J.shape[0]) *
                                       static_cast<std::size_t>(vq_J.shape[1]));
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
        + static_cast<std::int64_t>(sizeof(double)) + static_cast<std::int64_t>(sizeof(std::int64_t));
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
        RpaLriDetail::checked_write(ofs, record.payload.data(),
                                    record.payload.size() * sizeof(std::complex<double>),
                                    out_name);
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_abfs_overlap_v1");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::inverse_olp(const UnitCell& ucell,
                                    std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& overlap_abfs_abfs,
                                    const ModuleBase::Element_Basis_Index::IndexLNM& index_abfs_s)
{
    ModuleBase::TITLE("RPA_LRI", "inverse_olp");
    ModuleBase::timer::start("RPA_LRI", "inverse_olp");
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    size_t all_mu_s = 0;
    std::vector<int> mu_s_shift(ucell.nat);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_s_shift[I] = all_mu_s;
        all_mu_s += index_abfs_s[ucell.iat2it[I]].count_size;
    }
    RI::Tensor<std::complex<double>> olp_all = RI::Tensor<std::complex<double>>({all_mu_s, all_mu_s});
    for (int ik = 0; ik < nks_tot; ik++)
    {
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                // std::cout << "IJ: " << I << "," << J << std::endl;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;
                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                    {
                        olp_all(mu_s_shift[I] + ir, mu_s_shift[J] + ic) = JPp.second(ir, ic);
                    }
                }
            }
        }
        // for multi-mpi
        // for (int ir = 0; ir < all_mu_s; ir++)
        // {
        //     for (int ic = 0; ic < all_mu_s; ic++)
        //     {
        //         Parallel_Reduce::reduce_all<std::complex<double>>(olp_all(ir, ic));
        //     }
        // }

        // check Hermitian
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                auto delta = std::abs(olp_all(ir, ic) - std::conj(olp_all(ic, ir)));
                if (delta > 1e-10)
                {
                    std::cout << "Warning: olp_all is not Hermitian!" << std::endl;
                    std::cout << "ik,ir,ic: " << ik << "," << ir << "," << ic << std::endl;
                    std::cout << "delta(ir, ic): " << delta << std::endl;
                }
            }
        }
        // out_pure_ri_tensor("olp_all.txt", olp_all, 0.);
        auto olp_inv = LRI_CV_Tools::cal_I(olp_all,
                                           Inverse_Matrix<std::complex<double>>::Method::syev,
                                           this->info.shrink_LU_inv_thr);
        for (int ir = 0; ir < all_mu_s; ir++)
        {
            for (int ic = ir; ic < all_mu_s; ic++)
            {
                olp_inv(ic, ir) = std::conj(olp_inv(ir, ic));
            }
        }
        // out_pure_ri_tensor("olp_inv.txt", olp_inv, 0.);
        for (auto& Ip: overlap_abfs_abfs)
        {
            auto I = Ip.first;
            size_t mu_s_I = index_abfs_s[ucell.iat2it[I]].count_size;
            for (auto& JPp: Ip.second)
            {
                auto q = JPp.first.second;
                if (q != RI_Util::Vector3_to_array3(p_kv->kvec_c[ik]))
                    continue;
                auto J = JPp.first.first;
                auto mu_s_J = index_abfs_s[ucell.iat2it[J]].count_size;

                for (int ir = 0; ir < mu_s_I; ir++)
                {
                    for (int ic = 0; ic < mu_s_J; ic++)
                        JPp.second(ir, ic) = olp_inv(mu_s_shift[I] + ir, mu_s_shift[J] + ic);
                }
            }
        }
    }
    ModuleBase::timer::end("RPA_LRI", "inverse_olp");
}

// debug function
// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn,
//                                            RI::Tensor<std::complex<double>>& olp,
//                                            const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real()
//                    << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_pure_ri_tensor(const std::string fn, RI::Tensor<double>& olp, const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     int nr = olp.shape[0];
//     int nc = olp.shape[1];
//     size_t nnz = nr * nc;
//     fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//     fs << "%" << std::endl;

//     fs << nr << " " << nc << " " << nnz << std::endl;

//     for (int j = 0; j < nc; j++)
//     {
//         for (int i = 0; i < nr; i++)
//         {
//             auto v = olp(i, j);
//             if (fabs(v) > threshold)
//                 fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v << "\n";
//         }
//     }

//     fs.close();
// }

// template <typename T, typename Tdata>
// void RPA_LRI<T, Tdata>::out_ri_tensor(const std::string fn,
//                                       std::map<TA, std::map<TAq, RI::Tensor<std::complex<double>>>>& olp,
//                                       const double threshold)
// {
//     std::ofstream fs;
//     auto format = std::scientific;
//     int prec = 15;
//     fs.open(fn);
//     for (auto& IJq: olp)
//     {
//         int I = IJq.first;
//         for (auto& Jq: IJq.second)
//         {
//             int J = Jq.first.first;
//             auto q = Jq.first.second;
//             auto mat = Jq.second;
//             int nr = mat.shape[0];
//             int nc = mat.shape[1];
//             size_t nnz = nr * nc;
//             fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
//             fs << I << " " << J << " " << q.at(0) << " " << q.at(1) << " " << q.at(2) << std::endl;
//             fs << "%" << std::endl;

//             fs << nr << " " << nc << " " << nnz << std::endl;

//             for (int j = 0; j < nc; j++)
//             {
//                 for (int i = 0; i < nr; i++)
//                 {
//                     auto v = mat(i, j);
//                     if (fabs(v.real()) > threshold || fabs(v.imag()) > threshold)
//                         fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec)
//                            << v.real() << " " << std::showpoint << format << std::setprecision(prec) << v.imag()
//                            << "\n";
//                 }
//             }
//         }
//     }

//     fs.close();
// }

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_eigen_vector(const Parallel_Orbitals& parav, const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = this->runtime.input.nspin == 2 ? p_kv->get_nks() / 2 : p_kv->get_nks();
    const int npsin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    const int nbands = parav.get_wfc_global_nbands();
    const int nbasis = parav.get_wfc_global_nbasis();
    const std::size_t values_per_iw = static_cast<std::size_t>(nbands) * npsin_tmp;

    if (this->runtime.input.out_librpa_ver == 1)
    {
#ifdef __MPI
        ModuleBase::timer::start("RPA_LRI", "out_eigen_vector_v1_mpi_io");
        try
        {
            const MPI_Comm io_comm = parav.comm();
            if (io_comm == MPI_COMM_NULL)
            {
                throw std::runtime_error("KS eigenvector MPI-IO writer has no wavefunction communicator.");
            }
            int communicator_relation = MPI_UNEQUAL;
            RpaLriDetail::collective_mpi_check(
                io_comm,
                MPI_Comm_compare(io_comm, this->mpi_comm, &communicator_relation),
                "Failed to compare KS eigenvector MPI communicators");
            RpaLriDetail::collective_require(
                io_comm,
                communicator_relation == MPI_IDENT || communicator_relation == MPI_CONGRUENT,
                "KS eigenvector wavefunction and RPA communicators are inconsistent");
#endif

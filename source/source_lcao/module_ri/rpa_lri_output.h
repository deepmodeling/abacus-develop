    {
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            const ModuleBase::Vector3<double> position = ucell.atoms[it].tau[ia] * ucell.lat0; // in unit of Bohr
            ofs << std::setw(15) << position.x << std::setw(15) << position.y
                << std::setw(15) << position.z << std::setw(15) << it + 1 << std::endl;
        }
    }

    // Reader-v1 stores k-point sampling in bz_sample.txt. Keep the legacy
    // k-point block in stru_out.txt for reader-v0 consumers.
    if (this->runtime.input.out_librpa_reader_version != 1)
    {
        ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2]
            << std::setw(6) << std::endl;

        for (int ik = 0; ik != nks_tot; ik++)
        {
            const ModuleBase::Vector3<double> kpoint =
                p_kv->kvec_c[ik] * (ModuleBase::TWO_PI / ucell.lat0); // in unit of 1/Bohr
            ofs << std::setw(15) << kpoint.x << std::setw(15) << kpoint.y
                << std::setw(15) << kpoint.z << std::endl;
        }
        if (this->runtime.input.symmetry == "-1")
        {
            for (int ik = 0; ik != nks_tot; ++ik)
            {
                ofs << (ik + 1) << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bz_sampling(const UnitCell& ucell)
{
    if (this->runtime.rank != 0)
    {
        return;
    }

    ModuleBase::TITLE("DFT_RPA_interface", "out_bz_sampling");
    const int nks_tot = this->runtime.input.nspin == 2 ? static_cast<int>(p_kv->get_nks()) / 2 : p_kv->get_nks();
    const int n_coulomb_irreducible = nks_tot;

    std::ofstream ofs(this->runtime.input.rpa_outdir + "bz_sample.txt", std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open bz_sample.txt");
    }
    ofs << p_kv->nmp[0] << std::setw(6) << p_kv->nmp[1] << std::setw(6) << p_kv->nmp[2] << std::endl;
    ofs << nks_tot << std::setw(8) << n_coulomb_irreducible << std::endl;
    double weight_sum = 0.0;
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        weight_sum += p_kv->wk[ik];
    }
    if (weight_sum <= 0.0)
    {
        throw std::runtime_error("Cannot write bz_sample.txt with non-positive total k-point weight.");
    }
    for (int ik = 0; ik < nks_tot; ++ik)
    {
        const ModuleBase::Vector3<double> kvec_cartesian
            = p_kv->kvec_c[ik] * (ModuleBase::TWO_PI / ucell.lat0);
        ofs << std::setw(8) << ik + 1
            << std::setw(24) << std::scientific << std::setprecision(15)
            << (p_kv->wk[ik] / weight_sum)
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].x
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].y
            << std::setw(24) << std::scientific << std::setprecision(15) << p_kv->kvec_d[ik].z
            << std::setw(24) << std::scientific << std::setprecision(15) << kvec_cartesian.x
            << std::setw(24) << std::scientific << std::setprecision(15) << kvec_cartesian.y
            << std::setw(24) << std::scientific << std::setprecision(15) << kvec_cartesian.z
            << std::setw(8) << ik + 1
            << std::setw(8) << ik + 1
            << std::endl;
    }
    ofs.close();
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_bands(const elecstate::ElecState* pelec)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (this->runtime.rank != 0)
    {
        return;
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    const int nspin_tmp = this->runtime.input.nspin == 2 ? 2 : 1;
    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + "band_out.txt", std::ios::out);
    ofs << std::fixed << std::setprecision(15);
    ofs << nks_tot << std::endl;
    ofs << nspin_tmp << std::endl;
    ofs << this->runtime.input.nbands << std::endl;
    ofs << this->runtime.nlocal << std::endl;
    ofs << (pelec->eferm.ef / 2.0) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
            for (int ib = 0; ib != this->runtime.input.nbands; ib++)
            {
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << pelec->wg(ik + is * nks_tot, ib) * nks_tot
                    << std::setw(25) << pelec->ekb(ik + is * nks_tot, ib) / 2.0
                    << std::setw(25)
                    << pelec->ekb(ik + is * nks_tot, ib) * ModuleBase::Ry_to_eV << std::endl;
            }
        }
    }
    ofs.close();
    return;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_Cs(const UnitCell& ucell, std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Cs_in, std::string filename)
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
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_Cs_v1");
    ModuleBase::timer::start("RPA_LRI", "out_Cs_v1");

    struct CsRecord
    {
        int ia1 = 0;
        int ia2 = 0;
        int R[3] = {0, 0, 0};
        double max_abs = 0.0;
        std::int64_t offset = 0;
        RI::Tensor<Tdata>* tensor = nullptr;
        int nw1 = 0;
        int nw2 = 0;
        int naux = 0;
    };

    std::vector<CsRecord> records;
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
            CsRecord record;
            record.ia1 = I + 1;
            record.ia2 = J + 1;
            record.R[0] = JPp.first.second[0];
            record.R[1] = JPp.first.second[1];
            record.R[2] = JPp.first.second[2];
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
                        record.max_abs = std::max(record.max_abs, std::abs(RpaLriDetail::real_as_double(tmp_Cs(mu, i, j))));
                    }
                }
            }
            records.push_back(record);
        }
    }

    const std::int64_t nblocks = static_cast<std::int64_t>(records.size());
    const std::int64_t record_bytes = static_cast<std::int64_t>(5 * sizeof(std::int32_t)
        + sizeof(double) + sizeof(std::int64_t));
    std::int64_t offset = static_cast<std::int64_t>(3 * sizeof(std::int32_t) + 2 * sizeof(std::int64_t))
        + nblocks * record_bytes;
    for (auto& record: records)
    {
        record.offset = offset;
        const unsigned long long nw_product = RpaLriDetail::checked_mul_u64(
            static_cast<unsigned long long>(record.nw1),
            static_cast<unsigned long long>(record.nw2),
            "LibRPA v1 Cs block size");
        const unsigned long long values = RpaLriDetail::checked_mul_u64(
            nw_product,
            static_cast<unsigned long long>(record.naux),
            "LibRPA v1 Cs block size");
        const unsigned long long bytes = RpaLriDetail::checked_mul_u64(
            values,
            static_cast<unsigned long long>(sizeof(double)),
            "LibRPA v1 Cs block size");
        offset += RpaLriDetail::checked_i64_from_u64(bytes, "LibRPA v1 Cs block size");
    }

    std::stringstream ss;
    ss << filename << this->runtime.rank << ".dat";
    const std::string out_name = this->runtime.input.rpa_outdir + ss.str();
    std::ofstream ofs(out_name.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + out_name);
    }

    const std::int32_t marker = RpaLriDetail::LIBRPA_LRICOEF_V1_MARKER;
    const std::int32_t natom = static_cast<std::int32_t>(ucell.nat);
    const std::int32_t ncell = 0;
    RpaLriDetail::write_scalar(ofs, marker, out_name);
    RpaLriDetail::write_scalar(ofs, natom, out_name);
    RpaLriDetail::write_scalar(ofs, ncell, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    RpaLriDetail::write_scalar(ofs, nblocks, out_name);
    for (const auto& record: records)
    {
        const std::int32_t ia1 = record.ia1;
        const std::int32_t ia2 = record.ia2;
        const std::int32_t r0 = record.R[0];
        const std::int32_t r1 = record.R[1];
        const std::int32_t r2 = record.R[2];
        RpaLriDetail::write_scalar(ofs, ia1, out_name);
        RpaLriDetail::write_scalar(ofs, ia2, out_name);
        RpaLriDetail::write_scalar(ofs, r0, out_name);
        RpaLriDetail::write_scalar(ofs, r1, out_name);
        RpaLriDetail::write_scalar(ofs, r2, out_name);
        RpaLriDetail::write_scalar(ofs, record.max_abs, out_name);
        RpaLriDetail::write_scalar(ofs, record.offset, out_name);
    }
    for (const auto& record: records)
    {
        for (int i = 0; i != record.nw1; ++i)
        {
            for (int j = 0; j != record.nw2; ++j)
            {
                for (int mu = 0; mu != record.naux; ++mu)
                {
                    const double value = RpaLriDetail::real_as_double((*record.tensor)(mu, i, j));
                    RpaLriDetail::write_scalar(ofs, value, out_name);
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_Cs_v1");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_coulomb_k(const UnitCell& ucell,
                                      std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>& Vs,
                                      std::string filename,
                                      Exx_LRI<double>* exx_lri)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_coulomb_k");
    ModuleBase::timer::start("RPA_LRI", "out_coulomb_k");

    int all_mu = 0;
    std::vector<int> mu_shift(ucell.nat);
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    for (int I = 0; I != ucell.nat; I++)
    {
        mu_shift[I] = all_mu;
        all_mu += exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
    }
    const int nks_tot = this->runtime.input.nspin == 2 ? (int)p_kv->get_nks() / 2 : p_kv->get_nks();
    std::stringstream ss;
    ss << filename << (this->runtime.rank + 1) << ".txt";

    std::ofstream ofs;
    ofs.open(this->runtime.input.rpa_outdir + ss.str().c_str(), std::ios::out);

    ofs << nks_tot << std::endl;
    ofs << std::fixed << std::setprecision(15);
    for (auto& Ip: Vs)
    {
        auto I = Ip.first;
        size_t mu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);

        for (int ik = 0; ik != nks_tot; ik++)
        {
            std::map<size_t, RI::Tensor<std::complex<double>>> Vq_k_IJ;
            for (auto& JPp: Ip.second)
            {
                auto J = JPp.first.first;

                auto R = JPp.first.second;
                if (J < I)
                {
                    continue;
                }
                RI::Tensor<std::complex<double>> tmp_VR = RI::Global_Func::convert<std::complex<double>>(JPp.second);
                const double arg = 1 * (p_kv->kvec_c[ik] * (RI_Util::array3_to_Vector3(R) * ucell.latvec))
                                   * ModuleBase::TWO_PI; // latvec
                const std::complex<double> kphase = std::complex<double>(cos(arg), sin(arg));
                if (Vq_k_IJ[J].empty())
                {
                    Vq_k_IJ[J] = RI::Tensor<std::complex<double>>({tmp_VR.shape[0], tmp_VR.shape[1]});
                }
                Vq_k_IJ[J] = Vq_k_IJ[J] + tmp_VR * kphase;
            }
            for (auto& vq_Jp: Vq_k_IJ)
            {
                auto iJ = vq_Jp.first;
                auto& vq_J = vq_Jp.second;
                size_t nu_num = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[iJ]);
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[iJ] + 1
                    << "   " << mu_shift[iJ] + nu_num << std::endl;
                ofs << ik + 1 << "  " << p_kv->wk[ik] / 2.0 * this->runtime.input.nspin << std::endl;
                for (int i = 0; i != vq_J.data->size(); i++)
                {
                    ofs << std::setw(25) << (*vq_J.data)[i].real()
                        << std::setw(25) << (*vq_J.data)[i].imag() << std::endl;
                }
            }
        }
    }
    ofs.close();
    ModuleBase::timer::end("RPA_LRI", "out_coulomb_k");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::out_librpa_basis_v1(const UnitCell& ucell,
                                            Exx_LRI<double>* exx_lri,
                                            const std::string& aux_filename,
                                            const std::string& legacy_filename)
{
    if (this->runtime.rank != 0)
    {
        return;
    }
    ModuleBase::TITLE("DFT_RPA_interface", "out_librpa_basis_v1");

    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> type_naux(static_cast<std::size_t>(ucell.ntype), 0);
    std::vector<int> type_nw(static_cast<std::size_t>(ucell.ntype), 0);
    int total_wfc = 0;
    int total_aux = 0;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        type_naux[static_cast<std::size_t>(it)] = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(it);
        if (type_naux[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type auxiliary size.");
        }
        type_nw[static_cast<std::size_t>(it)] = ucell.atoms[it].nw;
        if (type_nw[static_cast<std::size_t>(it)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 basis output found a non-positive per-type wave-function size.");
        }
        total_wfc += type_nw[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
        total_aux += type_naux[static_cast<std::size_t>(it)] * ucell.atoms[it].na;
    }

    const auto wfc_l_nchi = RpaLriDetail::collect_wfc_l_nchi(ucell);
    const auto aux_l_nchi = RpaLriDetail::collect_abfs_l_nchi(exx_lri->abfs);
    const std::string wfc_filename = this->runtime.input.rpa_outdir + "wfc_basis.txt";
    const std::string aux_output_filename = this->runtime.input.rpa_outdir + aux_filename;
    const std::string legacy_output_filename = this->runtime.input.rpa_outdir + legacy_filename;
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_nw, wfc_l_nchi, wfc_filename);
    RpaLriDetail::write_librpa_split_basis_file(ucell, type_naux, aux_l_nchi, aux_output_filename);

    std::ofstream ofs(legacy_output_filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + legacy_output_filename);
    }
    ofs << std::setw(10) << ucell.ntype
        << std::setw(10) << total_wfc
        << std::setw(10) << total_aux
        << "    fallback" << std::endl;
    for (int it = 0; it != ucell.ntype; ++it)
    {
        ofs << std::setw(10) << it + 1
            << std::setw(10) << type_nw[static_cast<std::size_t>(it)]
            << std::setw(10) << type_naux[static_cast<std::size_t>(it)]
            << std::endl;
    }
    const auto write_legacy_shells = [&ofs](const std::vector<std::vector<int>>& shell_counts_by_type)
    {
        for (std::size_t itype = 0; itype < shell_counts_by_type.size(); ++itype)
        {
            const auto& shell_counts = shell_counts_by_type[itype];
            int nshell = 0;
            for (const int count : shell_counts)
            {
                nshell += count;
            }
            ofs << std::setw(10) << itype + 1
                << std::setw(10) << nshell
                << std::endl;
            for (std::size_t l = 0; l < shell_counts.size(); ++l)
            {
                for (int iradial = 0; iradial < shell_counts[l]; ++iradial)
                {
                    ofs << std::setw(10) << l << std::endl;
                }
            }

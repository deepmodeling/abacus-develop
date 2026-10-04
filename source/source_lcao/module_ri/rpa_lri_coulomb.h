        }
        total_basis += increment;
    }

    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs.good())
    {
        throw std::runtime_error("Failed to open " + filename);
    }
    ofs << std::setw(10) << ucell.ntype
        << std::setw(10) << total_basis
        << "    abacus" << std::endl;
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        ofs << std::setw(10) << itype + 1
            << std::setw(10) << type_sizes[static_cast<std::size_t>(itype)]
            << std::endl;
    }
    for (int itype = 0; itype < ucell.ntype; ++itype)
    {
        const auto& shell_counts = type_shell_counts[static_cast<std::size_t>(itype)];
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
    }
}

template<typename Tdata>
inline bool has_valid_matrix_shape(const RI::Tensor<Tdata>& tensor)
{
    return tensor.shape.size() == 2 && tensor.shape[0] > 0 && tensor.shape[1] > 0;
}

template<typename Tdata>
inline std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>
collect_local_irreducible_abf_blocks(
    const std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>>& period_blocks,
    const std::map<ModuleSymmetry::Tap, std::set<ModuleSymmetry::TC>>& irreducible_sector,
    std::size_t& n_skipped_irreducible_blocks)
{
    std::map<RI_2D_Comm::TA, std::map<RI_2D_Comm::TAC, RI::Tensor<Tdata>>> irreducible_blocks;
    n_skipped_irreducible_blocks = 0;
    for (const auto& irap_Rs: irreducible_sector)
    {
        const auto period_iter = period_blocks.find(irap_Rs.first.first);
        for (const auto& irR: irap_Rs.second)
        {
            const RI_2D_Comm::TAC ir_key = {irap_Rs.first.second, irR};
            if (period_iter == period_blocks.end())
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            const auto block_iter = period_iter->second.find(ir_key);
            if (block_iter == period_iter->second.end()
                || !has_valid_matrix_shape(block_iter->second))
            {
                ++n_skipped_irreducible_blocks;
                continue;
            }
            irreducible_blocks[irap_Rs.first.first][ir_key] = block_iter->second;
        }
    }
    return irreducible_blocks;
}

inline std::size_t sum_skipped_irreducible_blocks(const MPI_Comm& mpi_comm,
                                                  const std::size_t local_count)
{
    unsigned long long global_count = static_cast<unsigned long long>(local_count);
    unsigned long long reduced_count = global_count;
    MPI_Allreduce(&global_count, &reduced_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm);
    return static_cast<std::size_t>(reduced_count);
}

}

template <typename T, typename Tdata>
RPA_LRI<T, Tdata>::~RPA_LRI() = default;

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::postSCF(const UnitCell& ucell,
                                const MPI_Comm& mpi_comm_in,
                                const module_dm::DensityMatrix<T, Tdata>& dm,
                                const elecstate::ElecState* pelec,
                                const K_Vectors& kv,
                                const LCAO_Orbitals& orb,
                                const Parallel_Orbitals& parav,
                                const psi::Psi<T>& psi)
{
    ModuleBase::TITLE("RPA_LRI", "postSCF");
    ModuleBase::timer::start("RPA_LRI", "postSCF");
    ModuleBase::GlobalFunc::MAKE_DIR(this->runtime.input.rpa_outdir);
    this->ccp_rmesh_times_cut = this->runtime.input.rpa_ccp_rmesh_times;
    this->ccp_rmesh_times_ewald = this->info.ccp_rmesh_times; // should be `exx_ccp_rmesh_times`

    this->cal_postSCF_exx(dm, mpi_comm_in, ucell, kv, orb, parav);
    this->init(mpi_comm_in, kv, orb.cutoffs());
    this->out_bands(pelec);
    this->out_eigen_vector(parav, psi);
    this->out_struc(ucell);
    this->out_bz_sampling(ucell);

    std::cout << "rpa_pca_threshold: " << this->info.pca_threshold << std::endl;
    std::cout << "rpa_ccp_rmesh_times_cut: " << this->ccp_rmesh_times_cut << std::endl;
    std::cout << "rpa_ccp_rmesh_times_ewald: " << this->ccp_rmesh_times_ewald << std::endl;
    std::cout << "rpa_lcao_exx(Ha): " << std::fixed << std::setprecision(15) << exx_cut_coulomb->Eexx / 2.0 << std::endl;

    std::cout << "etxc(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etxc / 2.0 << std::endl;
    std::cout << "etot(Ha): " << std::fixed << std::setprecision(15) << pelec->f_en.etot / 2.0 << std::endl;
    std::cout << "Etot_without_rpa(Ha): " << std::fixed << std::setprecision(15)
              << (pelec->f_en.etot - pelec->f_en.etxc + exx_cut_coulomb->Eexx) / 2.0 << std::endl;
    exx_cut_coulomb.reset();
    RpaLriDetail::trim_malloc_cache();

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        cal_large_Cs(ucell, orb, kv);
        cal_abfs_overlap(ucell, orb, kv);
        RpaLriDetail::trim_malloc_cache();
    }
    this->output_ewald_coulomb(ucell, kv, orb);

    ModuleBase::timer::end("RPA_LRI", "postSCF");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::init(const MPI_Comm& mpi_comm_in, const K_Vectors& kv_in, const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("RPA_LRI", "init");
    ModuleBase::timer::start("RPA_LRI", "init");
    this->mpi_comm = mpi_comm_in;
    this->orb_cutoff_ = orb_cutoff;
    this->lcaos = exx_cut_coulomb->lcaos;
    this->p_kv = &kv_in;
    this->MGT = exx_cut_coulomb->MGT;

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        this->abfs_shrink = exx_cut_coulomb->abfs;
    }
    else
    {
        this->abfs = exx_cut_coulomb->abfs;
    }
    //	this->cv = std::move(exx_lri_rpa.cv);
    //    exx_lri_rpa.cv = exx_lri_rpa.cv;
    ModuleBase::timer::end("RPA_LRI", "init");
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::cal_postSCF_exx(const module_dm::DensityMatrix<T, Tdata>& dm,
                                        const MPI_Comm& mpi_comm_in,
                                        const UnitCell& ucell,
                                        const K_Vectors& kv,
                                        const LCAO_Orbitals& orb,
                                        const Parallel_Orbitals& parav)
{
    ModuleBase::TITLE("RPA_LRI", "cal_postSCF_exx");
    ModuleBase::timer::start("RPA_LRI", "cal_postSCF_exx");

    this->mpi_comm = mpi_comm_in;
    this->p_kv = &kv;
    this->orb_cutoff_ = orb.cutoffs();

    Mix_DMk_2D<T> mix_DMk_2D;
    bool exx_spacegroup_symmetry = (this->runtime.input.nspin < 4 && ModuleSymmetry::Symmetry::symm_flag == 1);
    if (exx_spacegroup_symmetry)
        {mix_DMk_2D.set_nks(kv.get_nkstot_nospin() * (this->runtime.input.nspin == 2 ? 2 : 1));}
    else
        {mix_DMk_2D.set_nks(kv.get_nks());}

    mix_DMk_2D.set_mixing_plain(1.0);

    // --------------------------------------------------------------------------------------
    // NOTE: ABFs are constructed here BEFORE symmetry processing to calculate the correct
    // abfs_Lmax for symrot.set_abfs_Lmax(). Previously, abfs_Lmax was obtained either from
    // exx_cut_coulomb->abfs_Lmax() (which was nullptr at that point) or from this->info.abfs_Lmax
    // (which defaults to 0). This caused abfs_Lmax to degrade to 0 in RPA + symmetry path,
    // leading to incorrect rotation matrices for ABFs with higher angular momentum.
    // --------------------------------------------------------------------------------------
    std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_for_lmax;
    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell, orb, this->lcaos, this->info.kmesh_times, this->info.shrink_abfs_pca_thr);
        if (!this->info.files_shrink_abfs.empty())
        {
            abfs_for_lmax = Exx_Abfs::IO::construct_abfs(abfs_for_lmax, orb, this->info.files_shrink_abfs, this->info.kmesh_times);
        }
    }
    else
    {
        this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs(orb, this->info.kmesh_times);
        abfs_for_lmax = Exx_Abfs::Construct_Orbs::abfs_same_atom(ucell, orb, this->lcaos, this->info.kmesh_times, this->info.pca_threshold);
        if (!this->info.files_abfs.empty())
        {
            abfs_for_lmax = Exx_Abfs::IO::construct_abfs(abfs_for_lmax, orb, this->info.files_abfs, this->info.kmesh_times);
        }
    }

    ModuleSymmetry::Symmetry_rotation symrot;
    if (exx_spacegroup_symmetry)
    {
        const std::array<Tcell, Ndim> period = RI_Util::get_Born_vonKarmen_period(kv);
        const auto& Rs = RI_Util::get_Born_von_Karmen_cells(period);
        symrot.find_irreducible_sector(ucell.symm, ucell.atoms, ucell.st, Rs, period, ucell.lat, this->runtime.global_out_dir);
        // set Lmax of the rotation matrices to max(l_ao, l_abf), to support rotation under ABF
        // NOTE: Using Exx_Abfs::Construct_Orbs::get_Lmax() to compute Lmax from the actual ABFs
        // instead of relying on exx_cut_coulomb->abfs_Lmax() (not yet initialized) or
        // this->info.abfs_Lmax (defaults to 0). This ensures correct Lmax for symmetry rotation.
        symrot.set_abfs_Lmax(Exx_Abfs::Construct_Orbs::get_Lmax(abfs_for_lmax));
        symrot.cal_Ms(kv, ucell, parav, this->runtime.input.nspin);
        // output Ts (symrot_R.txt) and Ms (symrot_k.txt)
        ModuleSymmetry::print_symrot_info_R(symrot, ucell.symm, ucell.lmax, Rs);
        ModuleSymmetry::print_symrot_info_k(symrot, kv, ucell);
        mix_DMk_2D.mix(symrot.restore_dm(kv, dm.get_dmk_vec(), parav), true);
    }
    else { mix_DMk_2D.mix(dm.get_dmk_vec(), true); }

    const std::vector<std::map<TA, std::map<TAC, RI::Tensor<Tdata>>>>
        Ds = RI_2D_Comm::split_m2D_ktoR<Tdata>(
            ucell,
            kv,
            mix_DMk_2D.get_DMk_out(),
            parav,
            this->runtime.input.nspin,
            exx_spacegroup_symmetry);

    if (!exx_cut_coulomb)
    {
        Exx_Info_RI local_info = this->info;
        local_info.ccp_rmesh_times = this->ccp_rmesh_times_cut;
        Exx_LRI<double>* new_exx = new Exx_LRI<double>(local_info);
        exx_cut_coulomb.reset(new_exx);
    }

    if (this->info.shrink_abfs_pca_thr >= 0.0)
    {
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        // This ensures consistency between the ABFs used for Lmax calculation and the ABFs
        // used for actual EXX computation.
        this->abfs_shrink = abfs_for_lmax;
        Exx_Abfs::Construct_Orbs::print_orbs_size(ucell, abfs_shrink, this->runtime.log);
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_shrink);
    }
    else
        // NOTE: Reuse abfs_for_lmax constructed earlier to avoid redundant ABFs construction.
        exx_cut_coulomb->init_spencer(mpi_comm_in, ucell, kv, orb, abfs_for_lmax);
    // cal C and V for exx
    Exx_LRI<double>* cut_coulomb = exx_cut_coulomb.get();
    this->output_cut_coulomb_cs(ucell, cut_coulomb);
    // cal CVCD
    if (exx_spacegroup_symmetry && this->runtime.input.exx_symmetry_realspace)
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav, &symrot);
    }
    else
    {
        exx_cut_coulomb->cal_exx_elec(Ds, ucell, parav);
    }
    // cout<<"postSCF_Eexx: "<<exx_lri_rpa.Eexx<<endl;
    ModuleBase::timer::end("RPA_LRI", "cal_postSCF_exx");
}

// if use shrink, output Coulomb and Cs_data in small abfs
// otherwise, output in normal abfs
template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_cut_coulomb_cs(const UnitCell& ucell, Exx_LRI<double>* exx_lri_rpa)
{
    ModuleBase::TITLE("RPA_LRI", "output_cut_coulomb_cs");
    ModuleBase::timer::start("RPA_LRI", "output_cut_coulomb_cs");

    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    std::cout << "Use rpa_ccp_rmesh_times=" << this->ccp_rmesh_times_cut << " to calculate cut Coulomb" << std::endl;
    // Shrink_ABFS_ORBITAL cannot exceed this angular momentum of MGT
    exx_lri_rpa->cal_cut_coulomb_cs(Vs_cut_IJR, Cs, ucell, this->runtime.input.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_cut, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(this->mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_cut_IJ
        = RI_2D_Comm::comm_map2_first(this->mpi_comm, Vs_cut_IJR, atoms00, atoms01);
    Vs_cut_IJR.clear();
    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_cut_IJ, period);
    if (this->runtime.input.out_librpa_reader_version == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_lri_rpa,
                                  use_shrink ? "aux_basis_s.txt" : "aux_basis.txt",
                                  use_shrink ? "wfc_basis_s.txt" : "basis_map.txt");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "V_cut_", exx_lri_rpa);
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_cut_", exx_lri_rpa);
    }
    Vs_period.clear();
    Vs_period.swap(tmp);

    this->Cs_period = RI::RI_Tools::cal_period(Cs, period);
    this->Cs_period = exx_lri_rpa->exx_lri.post_2D.set_tensors_map2(this->Cs_period);

    if (this->runtime.input.out_librpa_reader_version == 1)
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
        {
            this->out_Cs_v1(ucell, this->Cs_period, "Cs_shrink_");
        }
        else
        {
            this->out_Cs_v1(ucell, this->Cs_period, "Cs_");
        }
    }
    else
    {
        if (this->info.shrink_abfs_pca_thr >= 0.0)
            this->out_Cs(ucell, this->Cs_period, "Cs_shrinked_data_");
        else
            this->out_Cs(ucell, this->Cs_period, "Cs_data_");
    }
    Cs_period.clear();
    Cs_period.swap(tmp);

    ModuleBase::timer::end("RPA_LRI", "output_cut_coulomb_cs");
}

template <typename T, typename Tdata>
Conv_Coulomb_Pot_K::Coulomb_Method RPA_LRI<T, Tdata>::select_coulomb_basis_method_(Exx_LRI<double>* exx_lri) const
{
    if (exx_lri == nullptr || exx_lri->exx_objs.empty())
    {
        throw std::invalid_argument("Cannot select Coulomb basis method from an empty Exx_LRI object.");
    }
    return exx_lri->exx_objs.count(Conv_Coulomb_Pot_K::Coulomb_Method::Center2)
        ? Conv_Coulomb_Pot_K::Coulomb_Method::Center2
        : exx_lri->exx_objs.begin()->first;
}

template <typename T, typename Tdata>
std::vector<int> RPA_LRI<T, Tdata>::collect_atom_naux_(const UnitCell& ucell, Exx_LRI<double>* exx_lri) const
{
    const auto basis_method = this->select_coulomb_basis_method_(exx_lri);
    std::vector<int> atom_naux(static_cast<std::size_t>(ucell.nat), 0);
    for (int I = 0; I != ucell.nat; ++I)
    {
        atom_naux[static_cast<std::size_t>(I)]
            = exx_lri->exx_objs.at(basis_method).cv.get_index_abfs_size(ucell.iat2it[I]);
        if (atom_naux[static_cast<std::size_t>(I)] <= 0)
        {
            throw std::runtime_error("LibRPA v1 output found a non-positive per-atom auxiliary size.");
        }
    }
    return atom_naux;
}

template <typename T, typename Tdata>
void RPA_LRI<T, Tdata>::output_ewald_coulomb(const UnitCell& ucell, const K_Vectors& kv, const LCAO_Orbitals& orb)
{
    ModuleBase::TITLE("RPA_LRI", "output_ewald_coulomb");
    ModuleBase::timer::start("RPA_LRI", "output_ewald_coulomb");

    Exx_Info_RI local_info = this->info;
    local_info.ccp_rmesh_times = this->ccp_rmesh_times_ewald;
    if (!exx_full_coulomb)
    {
        Exx_LRI<double>* new_exx = new Exx_LRI<double>(local_info);
        exx_full_coulomb.reset(new_exx);
    }

    if (this->info.shrink_abfs_pca_thr >= 0.0)
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs_shrink);
    else
        exx_full_coulomb->init(mpi_comm, ucell, kv, orb, this->abfs);
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJR;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Cs;
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> tmp;
    exx_full_coulomb->cal_ewald_coulomb(Vs_full_IJR, Cs, ucell, this->runtime.input.out_ri_cv);
    // MPI: {ia0, {ia1, R}} to {ia0, ia1}
    std::vector<TA> atoms(ucell.nat);
    for (int iat = 0; iat < ucell.nat; ++iat)
        atoms[iat] = iat;
    const std::array<Tcell, Ndim> period_Vs
        = LRI_CV_Tools::cal_latvec_range<Tcell>(1 + this->ccp_rmesh_times_ewald, ucell, this->orb_cutoff_);
    const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA, TC>>>> list_As_Vs_atoms
        = RI::Distribute_Equally::distribute_atoms(mpi_comm, atoms, period_Vs, 2, false);
    const auto list_A0_pair_R = list_As_Vs_atoms.first;
    const auto list_A1_pair_R = list_As_Vs_atoms.second[0];
    std::set<TA> atoms00;
    std::set<TA> atoms01;
    for (const auto& I: list_A0_pair_R)
    {
        atoms00.insert(I);
    }
    for (const auto& JR: list_A1_pair_R)
    {
        atoms01.insert(JR.first);
    }
    std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> Vs_full_IJ
        = RI_2D_Comm::comm_map2_first(mpi_comm, Vs_full_IJR, atoms00, atoms01);
    Vs_full_IJR.clear();

    const std::array<Tcell, Ndim> period = {p_kv->nmp[0], p_kv->nmp[1], p_kv->nmp[2]};
    this->Vs_period = RI::RI_Tools::cal_period(Vs_full_IJ, period);
    if (this->runtime.input.out_librpa_reader_version == 1)
    {
        const bool use_shrink = this->info.shrink_abfs_pca_thr >= 0.0;
        this->out_librpa_basis_v1(ucell,
                                  exx_full_coulomb.get(),
                                  use_shrink ? "aux_basis_s.txt" : "aux_basis.txt",
                                  use_shrink ? "wfc_basis_s.txt" : "basis_map.txt");
        this->out_coulomb_k_v1(ucell, this->Vs_period, "V_full_", exx_full_coulomb.get());
    }
    else
    {
        this->out_coulomb_k(ucell, this->Vs_period, "coulomb_mat_", exx_full_coulomb.get());
    }
    Vs_period.clear();

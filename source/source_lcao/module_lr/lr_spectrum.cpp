#include "lr_spectrum.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "source_hamilt/module_gint/gint_interface.h"

template <typename T>
elecstate::DensityMatrix<T, T> LR::LR_Spectrum<T>::cal_transition_density_matrix(const int istate, const T* X_in, const bool need_R)
{
    const T* const X = X_in == nullptr ? this->X : X_in;
    const int offset_b = istate * ldim;    //start index of band istate
    elecstate::DensityMatrix<T, T> DM_trans(&this->pmat, this->nspin_x, this->kv.kvec_d, this->nk);
    for (int is = 0;is < this->nspin_x; ++is)
    {
        const int offset_x = offset_b + is * nk * this->pX[0].get_local_size();
        //1. transition density 
#ifdef __MPI
        std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X + offset_x, this->pX[is], psi_ks_vec[is], this->pc, this->naos, this->nocc[is], this->nvirt[is], this->pmat, (T)1.0 / (T)nk);
        // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
        std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(X + offset_x, this->psi_ks_vec[is], this->nocc[is], this->nvirt[is], (T)1.0 / (T)nk);
        // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
        for (int ik = 0;ik < this->nk;++ik) { DM_trans.set_DMK_pointer(ik + is * nk, dm_trans_2d[ik].data<T>()); }
    }
    if (need_R)
    {
        LR_Util::initialize_DMR(DM_trans, this->pmat, this->ucell, this->gd_, this->orb_cutoff_);
        DM_trans.cal_DMR();
    }
    return DM_trans;
}

inline void check_sum_rule(const double& osc_tot)
{
    GlobalV::ofs_running << "Total oscillator strength = " << osc_tot << std::endl;
    std::cout << "Total oscillator strength = " << osc_tot << std::endl;
    if (std::abs(osc_tot - 1.0) > 1e-3) {
        GlobalV::ofs_running << "The sum rule is not well satisfied, try more nvirt and nstates if needed." << std::endl;
        std::cout << "The sum rule is not well satisfied, try more nvirt and nstates if needed." << std::endl;
    }
}

template<>
ModuleBase::Vector3<double> LR::LR_Spectrum<double>::cal_transition_dipole_istate_length(const int istate)
{
    ModuleBase::Vector3<double> trans_dipole(0.0, 0.0, 0.0);
    // 1. transition density matrix
    const elecstate::DensityMatrix<double, double> DM_trans = this->cal_transition_density_matrix(istate);
    for (int is = 0;is < this->nspin_x;++is)
    {
        // 2. transition density
        double** rho_trans = nullptr;
        LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, this->rho_basis.nrxx);
        ModuleBase::GlobalFunc::ZEROS(rho_trans[0], this->rho_basis.nrxx);
        ModuleGint::cal_gint_rho({ DM_trans.get_DMR_vector().at(is) }, 1, rho_trans, false);

        // 3. transition dipole moment
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            int i = ir / (rho_basis.ny * rho_basis.nplane);
            int j = ir / rho_basis.nplane - i * rho_basis.ny;
            int k = ir % rho_basis.nplane + rho_basis.startz_current;
            ModuleBase::Vector3<double> rd(static_cast<double>(i) / rho_basis.nx, static_cast<double>(j) / rho_basis.ny, static_cast<double>(k) / rho_basis.nz);  //+1/2 better?
            rd -= ModuleBase::Vector3<double>(0.5, 0.5, 0.5);   //shift to the center of the grid (need ?)
            ModuleBase::Vector3<double> rc = rd * ucell.latvec * ucell.lat0; // real coordinate
            trans_dipole += rc * rho_trans[0][ir];
        }
        LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);
    }
    trans_dipole *= (ucell.omega / static_cast<double>(rho_basis.nxyz));   // dv
    trans_dipole *= static_cast<double>(this->nk);  // nk is divided inside DM_trans, now recover it
    if (this->nspin_x == 1) { trans_dipole *= sqrt(2.0); } // *2 for 2 spins, /sqrt(2) for the halfed dimension of X in the normalizaiton
    Parallel_Reduce::reduce_all(trans_dipole.x);
    Parallel_Reduce::reduce_all(trans_dipole.y);
    Parallel_Reduce::reduce_all(trans_dipole.z);
    return trans_dipole;
}

template<>
ModuleBase::Vector3<std::complex<double>> LR::LR_Spectrum<std::complex<double>>::cal_transition_dipole_istate_length(const int istate)
{

    //1. transition density matrix
    ModuleBase::Vector3<std::complex<double>> trans_dipole(0.0, 0.0, 0.0);
    const elecstate::DensityMatrix<std::complex<double>, std::complex<double>> DM_trans = this->cal_transition_density_matrix(istate);
    for (int is = 0;is < this->nspin_x;++is)
    {
        // 2. transition density
        double** rho_trans_real = nullptr;
        double** rho_trans_imag = nullptr;
        LR_Util::_allocate_2order_nested_ptr(rho_trans_real, 1, this->rho_basis.nrxx);
        LR_Util::_allocate_2order_nested_ptr(rho_trans_imag, 1, this->rho_basis.nrxx);

        elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(&this->pmat, 1, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(DM_trans_real_imag, this->pmat, this->ucell, this->gd_, this->orb_cutoff_);

        // real part
        LR_Util::get_DMR_real_imag_part(DM_trans, DM_trans_real_imag, ucell.nat, 'R');
        ModuleBase::GlobalFunc::ZEROS(rho_trans_real[0], this->rho_basis.nrxx);
        ModuleGint::cal_gint_rho(DM_trans_real_imag.get_DMR_vector(), 1, rho_trans_real, false);
        // LR_Util::print_grid_nonzero(rho_trans_real[0], this->rho_basis.nrxx, 10, "rho_trans");

        // imag part
        LR_Util::get_DMR_real_imag_part(DM_trans, DM_trans_real_imag, ucell.nat, 'I');
        ModuleBase::GlobalFunc::ZEROS(rho_trans_imag[0], this->rho_basis.nrxx);
        ModuleGint::cal_gint_rho(DM_trans_real_imag.get_DMR_vector(), 1, rho_trans_imag, false);
        // LR_Util::print_grid_nonzero(rho_trans_imag[0], this->rho_basis.nrxx, 10, "rho_trans");

        // 3. transition dipole moment
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            int i = ir / (rho_basis.ny * rho_basis.nplane);
            int j = ir / rho_basis.nplane - i * rho_basis.ny;
            int k = ir % rho_basis.nplane + rho_basis.startz_current;
            ModuleBase::Vector3<double> rd(static_cast<double>(i) / rho_basis.nx, static_cast<double>(j) / rho_basis.ny, static_cast<double>(k) / rho_basis.nz);  //+1/2 better?
            rd -= ModuleBase::Vector3<double>(0.5, 0.5, 0.5);   //shift to the center of the grid (need ?)
            ModuleBase::Vector3<double> rc = rd * ucell.latvec * ucell.lat0; // real coordinate
            ModuleBase::Vector3<std::complex<double>> rc_complex(rc.x, rc.y, rc.z);
            trans_dipole += rc_complex * std::complex<double>(rho_trans_real[0][ir], rho_trans_imag[0][ir]);
        }
        LR_Util::_deallocate_2order_nested_ptr(rho_trans_real, 1);
        LR_Util::_deallocate_2order_nested_ptr(rho_trans_imag, 1);
    }
    trans_dipole *= (ucell.omega / static_cast<double>(rho_basis.nxyz));   // dv
    trans_dipole *= static_cast<double>(this->nk);  // nk is divided inside DM_trans, now recover it
    if (this->nspin_x == 1) { trans_dipole *= sqrt(2.0); } // *2 for 2 spins, /sqrt(2) for the halfed dimension of X in the normalizaiton
    Parallel_Reduce::reduce_all(trans_dipole.x);
    Parallel_Reduce::reduce_all(trans_dipole.y);
    Parallel_Reduce::reduce_all(trans_dipole.z);
    return trans_dipole;
}

template<> double LR::LR_Spectrum<double>::cal_mean_squared_dipole(ModuleBase::Vector3<double> dipole)
{
    return dipole.norm2() / 3.;
}
template<> double LR::LR_Spectrum<std::complex<double>>::cal_mean_squared_dipole(ModuleBase::Vector3<std::complex<double>> dipole)
{
    // return dipole.norm2().real() / 3.;       // ModuleBase::Vector3::norm2 calculates x*x + y*y + z*z, but here we need x*x.conj() + y*y.conj() + z*z.conj()
    return (std::norm(dipole.x) + std::norm(dipole.y) + std::norm(dipole.z)) / 3.;
}

template<typename T>
void LR::LR_Spectrum<T>::cal_transition_dipoles_length()
{
    transition_dipole_.resize(nstate);
    this->mean_squared_transition_dipole_.resize(nstate);

    for (int istate = 0;istate < nstate;++istate)
    {
        transition_dipole_[istate] = cal_transition_dipole_istate_length(istate);
        mean_squared_transition_dipole_[istate] = cal_mean_squared_dipole(transition_dipole_[istate]);
    }
}

template<typename T>
void LR::LR_Spectrum<T>::oscillator_strength()
{
    ModuleBase::TITLE("LR::LR_Spectrum", "oscillator_strength");
    std::vector<double>& osc = this->oscillator_strength_;  // unit: Ry
    osc.resize(nstate, 0.0);
    double osc_tot = 0.0;
    for (int istate = 0;istate < nstate;++istate)
    {
        osc[istate] = this->mean_squared_transition_dipole_[istate] * this->omega[istate] * 2.;
        osc_tot += osc[istate] / 2.; //Ry to Hartree (1/2) 
    }
    const int nele = (this->nspin_x == 2) ? this->nocc[0] + this->nocc[1] : 2 * this->nocc[0];
    osc_tot /= this->nk * nele;
    check_sum_rule(osc_tot);
}

template<typename T>
void LR::LR_Spectrum<T>::optical_absorption_method1(const std::vector<double>& freq, const double eta)
{
    // ============test dipole================
    // this->cal_transition_dipoles_length();
    // this->write_transition_dipole(PARAM.globalv.global_out_dir + "dipole_length.dat");
    // this->cal_transition_dipoles_velocity();
    // this->write_transition_dipole(PARAM.globalv.global_out_dir + "dipole_velocity.dat");
    // exit(0);
    // ============test dipole================
    ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption");
    // 4*pi^2/V * mean_squared_dipole *delta(w-Omega_S)
    // = -8*pi*Omega_S/V * mean_squared_dipole * Im[1/[(w+i\eta)^2-\Omega_S^2]]
    // = -4*pi/V * oscilator_strength * Im[1/[(w+i\eta)^2-\Omega_S^2]]
    std::vector<double>& osc = this->oscillator_strength_;
    std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption.dat");

	if (GlobalV::MY_RANK == 0) 
	{ 
		ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; 
	}

    double FourPI_div_c = ModuleBase::FOUR_PI / 137.036;
    double fac = 4 * M_PI / ucell.omega * ModuleBase::e2 / this->nk;   // e2 for Ry to Hartree in the denominator

    for (int f = 0;f < freq.size();++f)
    {
        std::complex<double> f_complex = std::complex<double>(freq[f], eta);
        double abs = 0.0;
        // for (int i = 0;i < osc.size();++i) { abs += (osc[i] / (f_complex * f_complex - omega[i] * omega[i])).imag() * freq[f] * FourPI_div_c; }
        for (int i = 0;i < osc.size();++i) { abs += (osc[i] / (f_complex * f_complex - omega[i] * omega[i])).imag() * fac; }
        if (GlobalV::MY_RANK == 0) { ofs << freq[f] * ModuleBase::Ry_to_eV << "\t" << 91.126664 / freq[f] << "\t" << std::abs(abs) << std::endl; }
    }
    ofs.close();
}

template<typename T>
void LR::LR_Spectrum<T>::transition_analysis(const std::string& spintype)
{
    ModuleBase::TITLE("LR::LR_Spectrum", "transition_analysis");
    ModuleBase::timer::start("LR_Spectrum", "transition_analysis");
    std::ofstream ofs;
    std::ofstream ofs_k;
    const int nbands = nocc[0] + nvirt[0];
    const bool use_td_weight = (this->vmo_ptr != nullptr && LR_Util::tolower(this->gauge) == "velocity");
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(PARAM.globalv.global_out_dir + "trans_analysis_" + spintype + ".dat");
        ofs_k.open(PARAM.globalv.global_out_dir + "trans_kweight_" + spintype + ".dat");
        ofs << "==================================================================== \n";
        ofs << std::setw(40) << spintype << '\n';
        ofs << "==================================================================== \n";
        ofs << std::setw(8) << "State" << std::setw(30) << "Excitation Energy (Ry, eV)" <<
            std::setw(90) << "Transition dipole x, y, z (a.u.)" << std::setw(30) << "Oscillator strength(a.u.)" << '\n';
        ofs << "------------------------------------------------------------------------------------ \n";
        for (int istate = 0;istate < nstate;++istate)
            ofs << std::setw(8) << istate << std::setw(15) << std::setprecision(6) << omega[istate]
            << std::setw(15) << omega[istate] * ModuleBase::Ry_to_eV
            << std::setprecision(4) << std::setw(30) << transition_dipole_[istate].x 
            << std::setw(30) << transition_dipole_[istate].y << std::setw(30) << transition_dipole_[istate].z
            << std::setprecision(6) << std::setw(30) << oscillator_strength_[istate] << std::endl;
        ofs << "------------------------------------------------------------------------------------ " << std::endl;
        ofs << std::setw(8) << "State" << std::setw(20) << "Occupied orbital"
            << std::setw(20) << "Virtual orbital" << std::setw(30) << "Excitation amplitude"
            << std::setw(30) << "Excitation rate"
            << std::setw(10) << "k-point" << '\n';
        ofs << "------------------------------------------------------------------------------------ \n";

        ofs_k << "# Sum of exciton contribution of k-point.\n";
        ofs_k << "# weight1(k) = sum_{state,spin,occ,virt} |X(state,spin,k,occ,virt)|^2+|Y(state,spin,k,occ,virt)|^2.\n";
        ofs_k << "# weight2(k) = sum_{state,spin,direction,occ,virt} |td * X|^2 + |td *Y|^2, where td index as (spin,dir,k,occ,virt).\n";
        if (!use_td_weight)
        {
            ofs_k << "# NOTE: td-weighted statistics require abs_gauge velocity and a valid velocity_mo pointer.\n";
        }
        else { ofs_k << "# \n"; }
        
        ofs_k << "k-point" << std::setw(10) << "kx" << std::setw(12) << "ky" << std::setw(12) << "kz"
              << std::setw(12) << "weight1" << std::setw(12) << "weight2" << '\n';
    }
    std::vector<double> local_k_weight(nk, 0.0);
    std::vector<double> local_k_td_weight(nk, 0.0);
    T amp_X(0.0), amp_Y(0.0);
    // Communicate only amplitudes that will be written, in batches of NCOMM states.
    constexpr int NCOMM = 256;
    for (int istart = 0; istart < nstate; istart += NCOMM)
    {
        const int iend = std::min(istart + NCOMM, nstate);
        const std::size_t ncount_c = std::size_t(iend - istart) * this->gdim;
        if (ncount_c > std::numeric_limits<int>::max())
        {
            throw std::overflow_error("in transition_analysis: overflow converting to int!");
        }

        std::vector<int> local_indices;
        std::vector<T> local_amplitudes;
        for (int istate = istart; istate < iend; ++istate)
        {
            const int loffset_b = istate * ldim;
            const int goffset_b = (istate - istart) * gdim;
            for (int is = 0;is < nspin_x;++is)
            {
                const int loffset_bs = loffset_b + is * nk * pX[0].get_local_size();
                const int goffset_bs = goffset_b + is * nk * nocc[0] * nvirt[0];
                const int goffset_v_s = is * 3 * nk * nbands * nbands;
                const int eig_ks_spin_offset = is * nk * nbands;
                for (int ik = 0;ik < nk;++ik)
                {
                    const int loffset_x = loffset_bs + ik * pX[is].get_local_size();
                    const int goffset_x = goffset_bs + ik * nocc[is] * nvirt[is];
                    const int goffset_v_k = goffset_v_s + ik * nbands * nbands;
                    const int eig_ks_k_offset = eig_ks_spin_offset + ik * nbands;
                    for (int io = 0; io < pX[is].get_col_size(); ++io)
                    {
                        const int iocc = pX[is].local2global_col(io);
                        for (int iv = 0; iv < pX[is].get_row_size(); ++iv)
                        {
                            const int ivirt = pX[is].local2global_row(iv);
                            amp_X = X[loffset_x + io * pX[is].get_row_size() + iv];
                            if (this->is_full && this->Y ) {
                                amp_Y = Y[loffset_x + io * pX[is].get_row_size() + iv];
                            }
                            local_k_weight[ik] += std::norm(amp_X) + std::norm(amp_Y);
                            if (use_td_weight)
                            {
                                const int goffset_v = goffset_v_k + (ivirt + nocc[is]) * nbands + iocc;
                                const double gap = (eig_ks[eig_ks_k_offset + nocc[is] + ivirt]
                                                  - eig_ks[eig_ks_k_offset + iocc]) / ModuleBase::e2;  // Ry to Hartree
                                for (int id = 0; id < 3; ++id)
                                {
                                    const int v_index = goffset_v + id * nk * nbands * nbands;;
                                    std::complex<double> td_X = ModuleBase::IMAG_UNIT * this->vmo_ptr[v_index] * amp_X / gap;
                                    std::complex<double> td_Y = -1.0 * ModuleBase::IMAG_UNIT * std::conj(this->vmo_ptr[v_index]) * amp_Y / gap;
                                    // |<ik|v|ak>X_{aik}/(Ea-Ei)|^2 + |<ak|v|ik>Y_{aik}/(Ei-Ea)|^2
                                    if (this->nspin_x == 1) { td_X *= sqrt(2.0); td_Y *= sqrt(2.0); }
                                    local_k_td_weight[ik] += std::norm(td_X) + std::norm(td_Y);
                                }
                            }
                            if (std::abs(amp_X) > ana_thr) // only X components temporarily
                            {
                                local_indices.push_back(goffset_x + iocc * nvirt[is] + ivirt);
                                local_amplitudes.push_back(amp_X);
                            }
                        }
                    }
                }
            }
        }

        std::vector<int> all_indices;
        std::vector<T> all_amplitudes;
    #ifdef __MPI
        const int local_count = static_cast<int>(local_indices.size());
        int comm_size = 0;
        MPI_Comm_size(this->pX[0].comm(), &comm_size);
        std::vector<int> recv_counts;
        std::vector<int> displs;
        if (GlobalV::MY_RANK == 0)
        {
            recv_counts.resize(comm_size);
        }
        MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, this->pX[0].comm());
        if (GlobalV::MY_RANK == 0)
        {
            displs.resize(comm_size, 0);
            for (int ip = 1; ip < comm_size; ++ip)
            {
                displs[ip] = displs[ip - 1] + recv_counts[ip - 1];
            }
            const std::size_t total_count = std::size_t(displs.back()) + recv_counts.back();
            all_indices.resize(total_count);
            all_amplitudes.resize(total_count);
        }
        MPI_Gatherv(local_indices.data(), local_count, MPI_INT,
                    all_indices.data(), recv_counts.data(), displs.data(), MPI_INT,
                    0, this->pX[0].comm());
        MPI_Gatherv(local_amplitudes.data(), local_count, LR_Util::MPIType<T>::value(),
                    all_amplitudes.data(), recv_counts.data(), displs.data(), LR_Util::MPIType<T>::value(),
                    0, this->pX[0].comm());
    #else
        all_indices = std::move(local_indices);
        all_amplitudes = std::move(local_amplitudes);
    #endif
        if (GlobalV::MY_RANK != 0) continue; // only rank 0 write the analysis file

        std::vector<std::vector<std::pair<int, T>>> contributions(iend - istart);
        for (std::size_t i = 0; i < all_indices.size(); ++i)
        {
            const int ibatch = all_indices[i] / gdim;
            contributions[ibatch].emplace_back(all_indices[i] - ibatch * gdim, all_amplitudes[i]);
        }

        for (int istate = istart; istate < iend; ++istate)
        {
            auto& state_contributions = contributions[istate - istart];
            std::sort(state_contributions.begin(), state_contributions.end(),
                [](const std::pair<int, T>& l, const std::pair<int, T>& r) { return std::abs(l.second) > std::abs(r.second); });

            for (auto it = state_contributions.cbegin(); it != state_contributions.cend(); ++it)
            {
                auto pair_info = get_pair_info(it->first);
                const int& is = pair_info["ispin"];
                const std::string s = nspin_x == 2 ? (is == 0 ? "a" : "b") : "";
                ofs << std::setw(8) << (it == state_contributions.cbegin() ? std::to_string(istate) : " ")
                    << std::setw(20) << std::to_string(pair_info["iocc"] + 1) + s << std::setw(20) << std::to_string(pair_info["ivirt"] + nocc[is] + 1) + s// iocc and ivirt
                    << std::setw(30) << it->second
                    << std::setw(30) << std::norm(it->second)
                    << std::setw(10) << pair_info["ik"] + 1 << '\n';
            }
        }
    }
    std::vector<double> k_weight(nk, 0.0);
    std::vector<double> k_td_weight(nk, 0.0);
#ifdef __MPI
    MPI_Reduce(local_k_weight.data(), k_weight.data(), nk, MPI_DOUBLE, MPI_SUM, 0, this->pX[0].comm());
    MPI_Reduce(local_k_td_weight.data(), k_td_weight.data(), nk, MPI_DOUBLE, MPI_SUM, 0, this->pX[0].comm());
#else
    k_weight = std::move(local_k_weight);
    k_td_weight = std::move(local_k_td_weight);
#endif
    if (GlobalV::MY_RANK == 0)
    {
        ofs_k << std::fixed << std::setprecision(5);
        for (int ik = 0; ik < nk; ++ik)
        {
            ofs_k << std::setw(5) << ik + 1 << std::setw(12) << kv.kvec_d[ik].x
                  << std::setw(12) << kv.kvec_d[ik].y << std::setw(12) << kv.kvec_d[ik].z
                  << std::setw(12) << k_weight[ik] << std::setw(12) << k_td_weight[ik] << '\n';
        }
        ofs.close();
        ofs_k.close();
    }
    ModuleBase::timer::end("LR_Spectrum", "transition_analysis");
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LR::LR_Spectrum::transition_analysis");
}

template<typename T>
std::map<std::string, int> LR::LR_Spectrum<T>::get_pair_info(const int i)
{
    assert(i >= 0 && i < gdim);
    const int dim_spin0 = nk * nocc[0] * nvirt[0];
    const int ispin = (nspin_x == 2 && i >= dim_spin0) ? 1 : 0;
    const int ik = (i - ispin*dim_spin0) / (nocc[ispin] * nvirt[ispin]);
    const int ipair = (i - ispin*dim_spin0) - ik * nocc[ispin] * nvirt[ispin];
    const int iocc = ipair / nvirt[ispin];
    const int ivirt = ipair % nvirt[ispin];
    return  { {"ispin", ispin}, {"ik", ik}, {"iocc", iocc}, {"ivirt", ivirt} };
}

template<typename T>
void LR::LR_Spectrum<T>::write_transition_dipole(const std::string& filename)
{
    std::ofstream ofs(filename);
    ofs << "Transition dipole moment (a.u.)" << std::endl;
    ofs << std::setw(6) << "State" << std::setw(13) << "Energy (eV)" << std::setw(15) << "x" << std::setw(23) << "|x|^2" << std::setw(19) << "y" << std::setw(23) <<"|y|^2" << std::setw(19) << "z" << std::setw(23) <<"|z|^2" << std::setw(13) << "average" << std::endl;
    for (int istate = 0;istate < nstate;++istate)
    {
        ofs << std::setw(4) << istate << std::setw(13) << std::setprecision(6) << omega[istate] * ModuleBase::Ry_to_eV
        << std::setw(29) << transition_dipole_[istate].x << std::setw(13) << std::norm(transition_dipole_[istate].x)
        << std::setw(29) << transition_dipole_[istate].y << std::setw(13) << std::norm(transition_dipole_[istate].y)
        << std::setw(29) << transition_dipole_[istate].z << std::setw(13) << std::norm(transition_dipole_[istate].z)
        << std::setw(13) << mean_squared_transition_dipole_[istate] << std::endl;
    }
    ofs.close();
    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "LR::LR_Spectrum::write_transition_dipole " + filename);
}

template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;

#include "lr_spectrum.h"
#include "module_lr/utils/lr_util.h"
#include "module_parameter/parameter.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_base/parallel_reduce.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "module_lr/utils/lr_util_print.h"

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
        std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X + offset_x, this->pX[is], psi_ks[is], this->pc, this->naos, this->nocc[is], this->nvirt[is], this->pmat, (T)1.0 / (T)nk);
        // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
        std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(X + offset_x, this->psi_ks[is], this->nocc[is], this->nvirt[is], (T)1.0 / (T)nk);
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

template<typename T>
void LR::LR_Spectrum<T>::cal_gint_rho(double** rho, const int& nrxx)
{
    ModuleBase::GlobalFunc::ZEROS(rho[0], nrxx);
    Gint_inout inout_rho(rho, Gint_Tools::job_type::rho, 1, false);
    this->gint->cal_gint(&inout_rho);
}

inline void check_sum_rule(const double& osc_tot)
{
    if (std::abs(osc_tot - 1.0) > 1e-3) {
        GlobalV::ofs_running << "Warning: in LR_Spectrum::oscillator_strength, \
        the sum rule is not satisfied, try more nstates if needed.\n \
        Total oscillator strength = " + std::to_string(osc_tot) + "\n";
}
}

template<>
ModuleBase::Vector3<double> LR::LR_Spectrum<double>::cal_transition_dipole_istate_length(const int istate)
{
    ModuleBase::Vector3<double> trans_dipole(0.0, 0.0, 0.0);
    // 1. transition density matrix
    const elecstate::DensityMatrix<double, double>& DM_trans = this->cal_transition_density_matrix(istate);
    for (int is = 0;is < this->nspin_x;++is)
    {
        this->gint->transfer_DM2DtoGrid({ DM_trans.get_DMR_vector().at(is) });

        // 2. transition density
        double** rho_trans;
        LR_Util::_allocate_2order_nested_ptr(rho_trans, 1, this->rho_basis.nrxx);
        this->cal_gint_rho(rho_trans, this->rho_basis.nrxx);

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
    trans_dipole *= (ucell.omega / static_cast<double>(gint->get_ncxyz()));   // dv
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
    const elecstate::DensityMatrix<std::complex<double>, std::complex<double>>& DM_trans = this->cal_transition_density_matrix(istate);
    for (int is = 0;is < this->nspin_x;++is)
    {
        // 2. transition density
        double** rho_trans_real;
        double** rho_trans_imag;
        LR_Util::_allocate_2order_nested_ptr(rho_trans_real, 1, this->rho_basis.nrxx);
        LR_Util::_allocate_2order_nested_ptr(rho_trans_imag, 1, this->rho_basis.nrxx);

        elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(&this->pmat, 1, this->kv.kvec_d, this->nk);
        LR_Util::initialize_DMR(DM_trans_real_imag, this->pmat, this->ucell, this->gd_, this->orb_cutoff_);

        // real part
        LR_Util::get_DMR_real_imag_part(DM_trans, DM_trans_real_imag, ucell.nat, 'R');
        this->gint->transfer_DM2DtoGrid(DM_trans_real_imag.get_DMR_vector());
        this->cal_gint_rho(rho_trans_real, this->rho_basis.nrxx);
        // LR_Util::print_grid_nonzero(rho_trans_real[0], this->rho_basis.nrxx, 10, "rho_trans");

        // imag part
        LR_Util::get_DMR_real_imag_part(DM_trans, DM_trans_real_imag, ucell.nat, 'I');
        this->gint->transfer_DM2DtoGrid(DM_trans_real_imag.get_DMR_vector());
        this->cal_gint_rho(rho_trans_imag, this->rho_basis.nrxx);
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
    trans_dipole *= (ucell.omega / static_cast<double>(gint->get_ncxyz()));   // dv
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
        osc[istate] = this->mean_squared_transition_dipole_[istate] * this->eig[istate] * 2.;
        osc_tot += osc[istate] / 2.; //Ry to Hartree (1/2) 
    }
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
    if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
    double FourPI_div_c = ModuleBase::FOUR_PI / 137.036;
    double fac = 4 * M_PI / ucell.omega * ModuleBase::e2 / this->nk;   // e2 for Ry to Hartree in the denominator
    for (int f = 0;f < freq.size();++f)
    {
        std::complex<double> f_complex = std::complex<double>(freq[f], eta);
        double abs = 0.0;
        // for (int i = 0;i < osc.size();++i) { abs += (osc[i] / (f_complex * f_complex - eig[i] * eig[i])).imag() * freq[f] * FourPI_div_c; }
        for (int i = 0;i < osc.size();++i) { abs += (osc[i] / (f_complex * f_complex - eig[i] * eig[i])).imag() * fac; }
        if (GlobalV::MY_RANK == 0) { ofs << freq[f] * ModuleBase::Ry_to_eV << "\t" << 91.126664 / freq[f] << "\t" << std::abs(abs) << std::endl; }
    }
    ofs.close();
}

template<typename T>
void LR::LR_Spectrum<T>::transition_analysis(const std::string& spintype)
{
    ModuleBase::TITLE("LR::LR_Spectrum", "transition_analysis");
    std::ofstream& ofs = GlobalV::ofs_running;
    ofs << "==================================================================== " << std::endl;
    ofs << std::setw(40) << spintype << std::endl;
    ofs << "==================================================================== " << std::endl;
    ofs << std::setw(8) << "State" << std::setw(30) << "Excitation Energy (Ry, eV)" <<
        std::setw(90) << "Transition dipole x, y, z (a.u.)" << std::setw(30) << "Oscillator strength(a.u.)" << std::endl;
    ofs << "------------------------------------------------------------------------------------ " << std::endl;
    for (int istate = 0;istate < nstate;++istate)
        ofs << std::setw(8) << istate << std::setw(15) << std::setprecision(6) << eig[istate] << std::setw(15) << eig[istate] * ModuleBase::Ry_to_eV
        << std::setprecision(4) << std::setw(30) << transition_dipole_[istate].x << std::setw(30) << transition_dipole_[istate].y << std::setw(30) << transition_dipole_[istate].z
        << std::setprecision(6) << std::setw(30) << oscillator_strength_[istate] << std::endl;
    ofs << "------------------------------------------------------------------------------------ " << std::endl;
    ofs << std::setw(8) << "State" << std::setw(20) << "Occupied orbital"
        << std::setw(20) << "Virtual orbital" << std::setw(30) << "Excitation amplitude"
        << std::setw(30) << "Excitation rate"
        << std::setw(10) << "k-point" << std::endl;
    ofs << "------------------------------------------------------------------------------------ " << std::endl;
    for (int istate = 0;istate < nstate;++istate)
    {
        /// find the main contributions (> 0.5)
        const int loffset_b = istate * ldim;
        std::vector<T> X_full(gdim, T(0));// one-band, global
        for (int is = 0;is < nspin_x;++is)
        {
            const int loffset_bs = loffset_b + is * nk * pX[0].get_local_size();
            const int goffset_s = is * nk * nocc[0] * nvirt[0];
            for (int ik = 0;ik < nk;++ik)
            {
                const int loffset_x = loffset_bs + ik * pX[is].get_local_size();
                const int goffset_x = goffset_s + ik * nocc[is] * nvirt[is];
#ifdef __MPI
                LR_Util::gather_2d_to_full(this->pX[is], X + loffset_x, X_full.data() + goffset_x, false, nvirt[is], nocc[is]);
#endif
            }
        }
        std::map<double, int, std::greater<double>> abs_order;
        for (int i = 0;i < gdim;++i) { double abs = std::abs(X_full.at(i));if (abs > ana_thr) { abs_order[abs] = i; } }
        if (abs_order.size() > 0) {
            for (auto it = abs_order.cbegin();it != abs_order.cend();++it)
            {
                auto pair_info = get_pair_info(it->second);
                const int& is = pair_info["ispin"];
                const std::string s = nspin_x == 2 ? (is == 0 ? "a" : "b") : "";
                ofs << std::setw(8) << (it == abs_order.cbegin() ? std::to_string(istate) : " ")
                    << std::setw(20) << std::to_string(pair_info["iocc"] + 1) + s << std::setw(20) << std::to_string(pair_info["ivirt"] + nocc[is] + 1) + s// iocc and ivirt
                    << std::setw(30) << X_full.at(it->second)
                    << std::setw(30) << std::norm(X_full.at(it->second))
                    << std::setw(10) << pair_info["ik"] + 1 << std::endl;
            }
        }
    }
    ofs << "==================================================================== " << std::endl;
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
    ofs << std::setw(20) << "State" << std::setw(20) << "x" << std::setw(20) << "y" << std::setw(20) << "z" << std::setw(20) << "average" << std::endl;
    for (int istate = 0;istate < nstate;++istate)
    {
        ofs << std::setw(20) << istate << std::setw(20) << transition_dipole_[istate].x << std::setw(20)
            << transition_dipole_[istate].y << std::setw(20)
            << transition_dipole_[istate].z << std::setw(20)
            << mean_squared_transition_dipole_[istate] << std::endl;
    }
    ofs.close();
}

template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;
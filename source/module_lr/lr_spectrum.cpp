#include "lr_spectrum.h"
#include "module_lr/utils/lr_util.h"
#include "module_parameter/parameter.h"
#include "module_lr/dm_trans/dm_trans.h"
#include "module_base/parallel_reduce.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_hcontainer.h"
#include "module_lr/utils/lr_util_print.h"
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
void LR::LR_Spectrum<double>::oscillator_strength(Grid_Driver& gd, const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("LR::LR_Spectrum", "oscillator_strength");
    std::vector<double>& osc = this->oscillator_strength_;  // unit: Ry
    osc.resize(nstate, 0.0);
    double osc_tot = 0.0;
    elecstate::DensityMatrix<double, double> DM_trans(&this->pmat, 1, this->kv.kvec_d, this->nk);
    LR_Util::initialize_DMR(DM_trans, this->pmat, this->ucell, gd, orb_cutoff);
    this->transition_dipole_.resize(nstate, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    for (int istate = 0;istate < nstate;++istate)
    {
        const int offset_b = istate * ldim;    //start index of band istate
        for (int is = 0;is < this->nspin_x;++is)
        {
            const int offset_x = offset_b + is * nk * pX[0].get_local_size();
            //1. transition density 
#ifdef __MPI
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X + offset_x, this->pX[is], psi_ks[is], this->pc, this->naos, this->nocc[is], this->nvirt[is], this->pmat);
            // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(X + offset_x, this->psi_ks[is], this->nocc[is], this->nvirt[is]);
            // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
            for (int ik = 0;ik < this->nk;++ik) { DM_trans.set_DMK_pointer(ik, dm_trans_2d[ik].data<double>()); }
            DM_trans.cal_DMR();
            this->gint->transfer_DM2DtoGrid(DM_trans.get_DMR_vector());

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
                transition_dipole_[istate] += rc * rho_trans[0][ir];
            }
            LR_Util::_deallocate_2order_nested_ptr(rho_trans, 1);
        }
        transition_dipole_[istate] *= (ucell.omega / static_cast<double>(gint->get_ncxyz()));   // dv
        Parallel_Reduce::reduce_all(transition_dipole_[istate].x);
        Parallel_Reduce::reduce_all(transition_dipole_[istate].y);
        Parallel_Reduce::reduce_all(transition_dipole_[istate].z);
        osc[istate] = transition_dipole_[istate].norm2() * eig[istate] * 2. / 3.;
        osc_tot += osc[istate] / 2.; //Ry to Hartree (1/2) 
    }
    check_sum_rule(osc_tot);
}

template<>
void LR::LR_Spectrum<std::complex<double>>::oscillator_strength(Grid_Driver& gd, const std::vector<double>& orb_cutoff)
{
    ModuleBase::TITLE("LR::LR_Spectrum", "oscillator_strength");
    std::vector<double>& osc = this->oscillator_strength_;  // unit: Ry
    osc.resize(nstate, 0.0);
    double osc_tot = 0.0;
    elecstate::DensityMatrix<std::complex<double>, std::complex<double>> DM_trans(&this->pmat, 1, this->kv.kvec_d, this->nk);
    LR_Util::initialize_DMR(DM_trans, this->pmat, this->ucell, gd, orb_cutoff);
    elecstate::DensityMatrix<std::complex<double>, double> DM_trans_real_imag(&this->pmat, 1, this->kv.kvec_d, this->nk);
    LR_Util::initialize_DMR(DM_trans_real_imag, this->pmat, this->ucell, gd, orb_cutoff);

    this->transition_dipole_.resize(nstate, ModuleBase::Vector3<std::complex<double>>(0.0, 0.0, 0.0));
    for (int istate = 0;istate < nstate;++istate)
    {
        const int offset_b = istate * ldim;    //start index of band istate
        for (int is = 0;is < this->nspin_x;++is)
        {
            const int offset_x = offset_b + is * nk * pX[0].get_local_size();
            //1. transition density 
#ifdef __MPI
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X + offset_x, this->pX[is], psi_ks[is], this->pc, this->naos, this->nocc[is], this->nvirt[is], this->pmat, /*renorm_k=*/false, 1);
            // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
            std::vector<container::Tensor>  dm_trans_2d = cal_dm_trans_blas(X + offset_x, psi_ks[is], this->nocc[is], this->nvirt[is],/*renorm_k=*/false, 1);
            // if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
            for (int ik = 0;ik < this->nk;++ik) { DM_trans.set_DMK_pointer(ik, dm_trans_2d[ik].data<std::complex<double>>()); }
            // for (int ik = 0;ik < this->nk;++ik)
            //     LR_Util::print_tensor<std::complex<double>>(dm_trans_2d[ik], "1.DMK[ik=" + std::to_string(ik) + "]", dynamic_cast<Parallel_2D*>(&this->pmat));
            DM_trans.cal_DMR();

            // 2. transition density
            double** rho_trans_real;
            double** rho_trans_imag;
            LR_Util::_allocate_2order_nested_ptr(rho_trans_real, 1, this->rho_basis.nrxx);
            LR_Util::_allocate_2order_nested_ptr(rho_trans_imag, 1, this->rho_basis.nrxx);
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
                transition_dipole_[istate] += rc_complex * std::complex<double>(rho_trans_real[0][ir], rho_trans_imag[0][ir]);
            }
            LR_Util::_deallocate_2order_nested_ptr(rho_trans_real, 1);
            LR_Util::_deallocate_2order_nested_ptr(rho_trans_imag, 1);
        }
        transition_dipole_[istate] *= (ucell.omega / static_cast<double>(gint->get_ncxyz()));   // dv
        Parallel_Reduce::reduce_all(transition_dipole_[istate].x);
        Parallel_Reduce::reduce_all(transition_dipole_[istate].y);
        Parallel_Reduce::reduce_all(transition_dipole_[istate].z);
        auto norm2 = [](const ModuleBase::Vector3<std::complex<double>>& v) -> double
            {
                return v.x.real() * v.x.real() + v.y.real() * v.y.real() + v.z.real() * v.z.real()
                    + v.x.imag() * v.x.imag() + v.y.imag() * v.y.imag() + v.z.imag() * v.z.imag();
            };
        osc[istate] = norm2(transition_dipole_[istate]) * eig[istate] * 2. / 3.;
        osc_tot += osc[istate] / 2.;   // Ry to Hartree (1/2)
    }
    check_sum_rule(osc_tot);
}
template<typename T>
void LR::LR_Spectrum<T>::optical_absorption(const std::vector<double>& freq, const double eta, const std::string& spintype)
{
    ModuleBase::TITLE("LR::LR_Spectrum", "optical_absorption");
    std::vector<double>& osc = this->oscillator_strength_;
    std::ofstream ofs(PARAM.globalv.global_out_dir + "absorption_" + spintype + ".dat");
    if (GlobalV::MY_RANK == 0) { ofs << "Frequency (eV) | wave length(nm) | Absorption (a.u.)" << std::endl; }
    double FourPI_div_c = ModuleBase::FOUR_PI / 137.036;
    for (int f = 0;f < freq.size();++f)
    {
        std::complex<double> f_complex = std::complex<double>(freq[f], eta);
        double abs = 0.0;
        for (int i = 0;i < osc.size();++i) { abs += (osc[i] / (f_complex * f_complex - eig[i] * eig[i])).imag() * freq[f] * FourPI_div_c; }
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
        std::setw(45) << "Transition dipole x, y, z (a.u.)" << std::setw(30) << "Oscillator strength(a.u.)" << std::endl;
    ofs << "------------------------------------------------------------------------------------ " << std::endl;
    for (int istate = 0;istate < nstate;++istate)
        ofs << std::setw(8) << istate << std::setw(15) << std::setprecision(6) << eig[istate] << std::setw(15) << eig[istate] * ModuleBase::Ry_to_eV
        << std::setw(15) << transition_dipole_[istate].x << std::setw(15) << transition_dipole_[istate].y << std::setw(15) << transition_dipole_[istate].z
        << std::setw(30) << oscillator_strength_[istate] << std::endl;
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

template class LR::LR_Spectrum<double>;
template class LR::LR_Spectrum<std::complex<double>>;
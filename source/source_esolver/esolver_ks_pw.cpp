#include "esolver_ks_pw.h"
#include "source_io/module_wf/exx_source_io.h"
#include "source_io/module_wf/read_wfc_pw.h"
#include "source_base/module_out/filename.h"

#include <cerrno>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "source_estate/elecstate_pw.h"
#include "source_estate/module_charge/chg_symm.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_hsolver/diago_params.h"
#include "source_hamilt/hamilt_hs_adapter.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/force_pw.h"
#include "source_pw/module_pwdft/hamilt_pw.h"
#include "source_pw/module_pwdft/stress_pw.h"

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#endif

#include "source_estate/module_charge/chg_routine.h" // use charge mixing, mohan add 20251006
#include "source_estate/setup_estate_pw.h"         // mohan add 20251005
#include "source_hamilt/module_xc/general_exx_info.h" // for General_Exx_Info type used via general_exx_info_
#include "source_io/module_ctrl/ctrl_output_pw.h"  // mohan add 20250927
#include "source_pw/module_pwdft/deltaspin_pw.h"   // mohan add 20250309
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/setup_pot.h"      // mohan add 20250929
#include "source_pw/module_pwdft/update_cell_pw.h" // mohan add 20250309
#include "source_pw/module_pwdft/setup_dftu_pw.h"  // mohan add 20250309

namespace ModuleESolver
{
namespace
{
std::string exx_restart_configuration(const Input_para& inp, const General_Exx_Info& info)
{
    std::ostringstream out;
    out << std::setprecision(17) << inp.dft_functional << ':' << inp.ecutwfc << ':'
        << info.ecut_exx << ':' << info.hybrid_alpha << ':' << info.gamma_extrapolation;
    for (const auto& kernel : info.coulomb_param)
    {
        out << ':' << static_cast<int>(kernel.first);
        for (const auto& term : kernel.second)
        {
            for (const auto& parameter : term)
            {
                out << ':' << parameter.first << '=' << parameter.second;
            }
        }
    }
    return out.str();
}
}


template <typename T, typename Device>
ESolver_KS_PW<T, Device>::ESolver_KS_PW()
{
    this->classname = "ESolver_KS_PW";
    this->basisname = "PW";
    // PW basis: the DFT+U object is the base class (no LCAO orbitals).
    this->dftu_.reset(new Plus_U_Base());
}

template <typename T, typename Device>
ESolver_KS_PW<T, Device>::~ESolver_KS_PW()
{
    //****************************************************
    // do not add any codes in this deconstructor funcion
    //****************************************************
    // delete Hamilt
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }

    // delete exx_helper
    if (this->exx_helper != nullptr)
    {
        delete this->exx_helper;
        this->exx_helper = nullptr;
    }

    // mohan add 2025-10-12
    this->stp.clean();
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::allocate_hamilt(const UnitCell& ucell)
{
    this->p_hamilt = new hamilt::HamiltPW<T, Device>(this->pelec->pot,
                                                     this->pw_wfc,
                                                     &this->kv,
                                                     &this->ppcell,
                                                     this->dftu_.get(),
                                                     &ucell,
                                                     &this->general_exx_info_);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_all_runners(BaseCell& basecell, const Input_para& inp)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ESolver_KS::before_all_runners(ucell, inp);

    if (inp.calculation == "nscf" && this->general_exx_info_.cal_exx)
    {
        const auto fock = this->general_exx_info_.coulomb_param.find(Conv_Coulomb_Pot_K::Coulomb_Type::Fock);
        const bool has_fock = fock != this->general_exx_info_.coulomb_param.end() && !fock->second.empty();
        if (inp.device != "cpu" || inp.kpar != 1 || inp.bndpar != 1 || inp.nspin == 4 || inp.exxace
            || inp.cal_force || inp.cal_stress || inp.exx_gamma_extrapolation || inp.symmetry == "1" || inp.init_wfc == "file" || has_fock)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "Hybrid NSCF currently requires CPU, kpar=bndpar=1, nspin=1/2, symmetry=-1/0, exxace=false, exx_gamma_extrapolation=false, screened exchange, cal_force=cal_stress=false and fresh target wavefunctions");
        }
    }

    //! setup and allocation for pelec, potentials, etc.
    elecstate::setup_estate_pw(ucell,
                               this->kv,
                               this->sf,
                               this->pelec,
                               this->chr,
                               this->locpp,
                               this->ppcell,
                               this->vsep_cell,
                               this->pw_wfc,
                               this->pw_rho,
                               this->pw_rhod,
                               this->pw_big,
                               this->solvent,
                               inp);

    if (inp.calculation == "scf" && this->general_exx_info_.cal_exx && inp.out_wfc_pw == 2
        && inp.kpar == 1 && inp.bndpar == 1 && this->pw_wfc->poolrank == 0)
    {
        const std::string checkpoint = PARAM.globalv.global_out_dir + "EXX_SOURCE";
        if (std::remove(checkpoint.c_str()) != 0 && errno != ENOENT)
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_PW", "Cannot invalidate the previous EXX source checkpoint");
        }
    }

    this->stp.before_runner(ucell, this->kv, this->sf, *this->pw_wfc, this->ppcell.lmaxkb, *this->inp_);

    ModuleBase::GlobalFunc::DONE(GlobalV::ofs_running, "INIT BASIS");

    //! Create exx_helper based on device and precision
    const bool is_gpu = (inp.device == "gpu");
    const bool is_single = (inp.precision == "single");

#if ((defined __CUDA) || (defined __ROCM))
    if (is_gpu)
    {
        if (is_single)
        {
            this->exx_helper = new Exx_Helper<std::complex<float>, base_device::DEVICE_GPU>();
        }
        else
        {
            this->exx_helper = new Exx_Helper<std::complex<double>, base_device::DEVICE_GPU>();
        }
    }
    else
#endif
    {
        if (is_single)
        {
            this->exx_helper = new Exx_Helper<std::complex<float>, base_device::DEVICE_CPU>();
        }
        else
        {
            this->exx_helper = new Exx_Helper<std::complex<double>, base_device::DEVICE_CPU>();
        }
    }

    //! Initialize exx pw
    this->exx_helper->init(ucell, inp, this->pelec->wg, this->general_exx_info_);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::before_scf(UnitCell& ucell, const int istep)
{
    ModuleBase::TITLE("ESolver_KS_PW", "before_scf");
    ModuleBase::timer::start("ESolver_KS_PW", "before_scf");

    ESolver_KS::before_scf(ucell, istep);

    //! Init variables (once the cell has changed)
    pw::update_cell_pw(ucell, this->ppcell, this->kv, this->pw_wfc, *this->inp_);

    if (ucell.cell_parameter_updated)
    {
        this->stp.p_psi_init->prepare_init(this->inp_->pw_seed, istep);
    }

    //! Init Hamiltonian (cell changed)
    //! Operators in HamiltPW should be reallocated once cell changed
    //! delete Hamilt if not first scf
    if (this->p_hamilt != nullptr)
    {
        delete this->p_hamilt;
        this->p_hamilt = nullptr;
    }

    //! Allocate HamiltPW
    this->allocate_hamilt(ucell);

    //! Setup potentials (local, non-local, sc, +U, DFT-1/2)
    // note: init DFT+U is done here for pw basis for every scf iteration, however,
    // init DFT+U is done in "before_all_runners" in LCAO basis. This should be refactored, mohan note 2025-11-06
    pw::setup_pot(istep,
                  ucell,
                  this->kv,
                  this->sf,
                  this->pelec,
                  this->Pgrid,
                  this->chr,
                  this->locpp,
                  this->ppcell,
                  *this->dftu_,
                  this->vsep_cell,
                  this->stp.template get_psi_t<T, Device>(),
                  this->p_hamilt,
                  this->pw_wfc,
                  this->pw_rhod,
                  PARAM.globalv.global_out_dir,
                  *this->inp_);

    // setup psi (electronic wave functions)
    this->stp.init(this->p_hamilt);

    //! Setup EXX helper for Hamiltonian and psi
    exx_helper->before_scf(this->p_hamilt, this->stp.template get_psi_t<T, Device>(), *this->inp_, this->general_exx_info_);

    if (this->inp_->calculation == "nscf" && this->general_exx_info_.cal_exx)
    {
        this->prepare_exx_nscf(ucell, PARAM.globalv.global_readin_dir);
    }

    ModuleBase::timer::end("ESolver_KS_PW", "before_scf");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::prepare_exx_nscf(const UnitCell& ucell, const std::string& readin_dir)
{
    ModuleBase::timer::start("ESolver_KS_PW", "prepare_exx_nscf");
    const Input_para& inp = *this->inp_;
    const std::string configuration = exx_restart_configuration(inp, this->general_exx_info_);
    const std::string filename = readin_dir + "EXX_SOURCE";
    std::ifstream in(filename);
    try
    {
        ModuleIO::read_exx_source(in, exx_source_points_, exx_source_weights_, inp.nspin,
                                 this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz, configuration);
    }
    catch (const std::exception& error)
    {
        ModuleBase::WARNING_QUIT("ESolver_KS_PW", error.what());
    }
    const int nqs = exx_source_points_.get_nks();
    int nq = exx_source_points_.get_nkstot_nospin();
    exx_source_points_.para_k.kinfo(nq, 1, 0, this->pw_wfc->poolrank, this->pw_wfc->poolnproc, inp.nspin);
    exx_source_basis_.reset(new ModulePW::PW_Basis_K(inp.device, inp.precision));
#ifdef __MPI
    exx_source_basis_->initmpi(this->pw_wfc->poolnproc, this->pw_wfc->poolrank, this->pw_wfc->pool_world);
#endif
    exx_source_basis_->initgrids(ucell.lat0, ucell.latvec,
                               this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz);
    exx_source_basis_->initparameters(false, inp.ecutwfc, nqs, exx_source_points_.kvec_d.data());
    exx_source_basis_->fft_bundle.initfftmode(inp.fft_mode);
    exx_source_basis_->setuptransform();
    exx_source_basis_->collect_local_pw(inp.erf_ecut, inp.erf_height, inp.erf_sigma);
    for (int iq = 0; iq < nqs; ++iq)
    {
        exx_source_points_.ngk[iq] = exx_source_basis_->npwk[iq];
        exx_source_points_.kvec_c[iq] = exx_source_points_.kvec_d[iq] * ucell.G;
    }
    const int nbands = exx_source_weights_.nc;
    const int nbasis = exx_source_basis_->npwk_max;
    exx_source_psi_.reset(new psi::Psi<T, Device>(nqs, nbands, nbasis, exx_source_points_.ngk, true));
    for (int iq = 0; iq < nqs; ++iq)
    {
        const std::string wfc_filename = ModuleIO::filename_output(readin_dir, "wf", "pw", iq,
            exx_source_points_.ik2iktot, inp.nspin, nqs, 2, false, false, -1);
        ModuleBase::ComplexMatrix wfc(nbands, nbasis);
        ModuleIO::read_wfc_pw(wfc_filename, exx_source_basis_.get(), this->pw_wfc->poolrank,
                             this->pw_wfc->poolnproc, nbands, 1, iq, iq, nqs, wfc);
        exx_source_psi_->fix_k(iq);
        std::vector<T> coefficients(nbands * nbasis);
        for (int ib = 0; ib < nbands; ++ib)
        {
            for (int ig = 0; ig < nbasis; ++ig)
            {
                coefficients[ib * nbasis + ig] = static_cast<T>(wfc(ib, ig));
            }
        }
        base_device::memory::synchronize_memory_op<T, Device, base_device::DEVICE_CPU>()(
            exx_source_psi_->get_pointer(), coefficients.data(), coefficients.size());
    }
    auto* helper = static_cast<Exx_Helper<T, Device>*>(this->exx_helper);
    helper->op_exx->set_source(*exx_source_basis_, exx_source_points_, *exx_source_psi_, exx_source_weights_);
    ModuleBase::timer::end("ESolver_KS_PW", "prepare_exx_nscf");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_init(UnitCell& ucell, const int istep, const int iter)
{
    ESolver_KS::iter_init(ucell, istep, iter);

    // query DeltaSpin convergence for U-ramping; treat as converged when disabled
    bool mag_converged = true;
    if (this->inp_->sc_mag_switch)
    {
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        mag_converged = sc.mag_converged();
    }

    module_charge::chgmixing_ks_pw(iter, this->p_chgmix, *this->dftu_, mag_converged, *this->inp_);

    // mohan move harris functional here, 2012-06-05
    // use 'rho(in)' and 'v_h and v_xc'(in)
    this->pelec->f_en.deband_harris = this->pelec->cal_delta_eband(ucell);

    // update local occupations for DFT+U
    // should before lambda loop in DeltaSpin
    DFTU_BASE::iter_init_dftu_pw(iter,
                          istep,
                          *this->dftu_,
                          this->stp.template get_psi_t<T, Device>(),
                          this->pelec->wg,
                          ucell,
                          this->p_chgmix,
                          this->kv.isk.data());

    // mohan add 2025-11: push DFT+U energy from Plus_U instance to ElecState
    if (this->inp_->dft_plus_u)
    {
        this->pelec->set_dftu_energy(this->dftu_->get_energy());
    }
}

// Temporary, it should be replaced by hsolver later.
template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr)
{
    ModuleBase::timer::start("ESolver_KS_PW", "hamilt2rho_single");

    // reset energy
    this->pelec->f_en.eband = 0.0;
    this->pelec->f_en.demet = 0.0;

    // setup diagonalization parameters
    hsolver::setup_diago_params_pw<T, Device>(istep, iter, ethr, *this->inp_);

    bool skip_charge = this->inp_->calculation == "nscf" ? true : false;

    // run the inner lambda loop to contrain atomic moments with the DeltaSpin method
    bool skip_solve = pw::run_deltaspin_lambda_loop(iter - 1, this->drho, *this->inp_, GlobalV::ofs_running);
    if (skip_solve)
    {
        // Fetch the most recent DeltaSpin RMS for display in the SCF iteration table.
        spinconstrain::SpinConstrain<std::complex<double>>& sc
            = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
        this->ds_rms_ = sc.get_last_rms_error();
    }

    if (!skip_solve)
    {
        hsolver::HSolverPW<T, Device> hsolver_pw_obj(this->pw_wfc,
                                                     this->inp_->calculation,
                                                     this->inp_->basis_type,
                                                     this->inp_->ks_solver,
                                                     PARAM.globalv.use_uspp,
                                                     this->inp_->nspin,
                                                     hsolver::DiagoIterAssist<T, Device>::SCF_ITER,
                                                     hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX,
                                                     hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                                     hsolver::DiagoIterAssist<T, Device>::need_subspace,
                                                     this->inp_->nbands,
                                                     this->inp_->diago_smooth_ethr,
                                                     this->inp_->pw_diag_ndim,
                                                     this->inp_->diag_subspace,
                                                     this->inp_->nb2d,
                                                     this->inp_->use_k_continuity);

        // the iterative eigensolvers see the Hamiltonian only through this operator
        hamilt::HamiltHSOperator<T, Device> op(static_cast<hamilt::Hamilt<T, Device>*>(this->p_hamilt), this->pw_wfc);
        hsolver_pw_obj.solve(op,
                             *this->stp.template get_psi_t<T, Device>(),
                             this->pelec,
                             this->pelec->ekb.c,
                             GlobalV::RANK_IN_POOL,
                             GlobalV::NPROC_IN_POOL,
                             GlobalV::ofs_running,
                             skip_charge);
    }

    // symmetrize the charge density
    module_charge::symmetrize_rho(this->inp_->nspin, this->chr, this->pw_rhod, ucell.symm);

    ModuleBase::timer::end("ESolver_KS_PW", "hamilt2rho_single");
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver)
{
    // Related to EXX
    bool cal_exx = general_exx_info_.cal_exx;
    double hybrid_alpha = general_exx_info_.hybrid_alpha;
    if (cal_exx && this->inp_->calculation != "nscf" && !exx_helper->get_op_first_iter())
    {
        this->pelec->set_exx(exx_helper->cal_exx_energy(this->stp.template get_psi_t<T, Device>()),
                             cal_exx,
                             hybrid_alpha);
    }

    // deband is calculated from "output" charge density
    this->pelec->f_en.deband = this->pelec->cal_delta_eband(ucell);

    // Call iter_finish() of ESolver_KS
    ESolver_KS::iter_finish(ucell, istep, iter, conv_esolver);

    // D in USPP needs vloc, thus needs update when veff updated
    // calculate the effective coefficient matrix for non-local
    // pp projectors, liuyu 2023-10-24
    if (PARAM.globalv.use_uspp)
    {
        ModuleBase::matrix veff = this->pelec->pot->get_eff_v();
        this->ppcell.cal_effective_D(veff, this->pw_rhod, ucell);
    }

    // Handle EXX-related operations after SCF iteration
    exx_helper->iter_finish(this->pelec,
                            &this->chr,
                            this->stp.template get_psi_t<T, Device>(),
                            ucell,
                            *this->inp_,
                            conv_esolver,
                            iter);

    // check if oscillate for delta_spin method
    pw::check_deltaspin_oscillation(iter, this->drho, this->p_chgmix, *this->inp_);

    // the output quantities
    ModuleIO::ctrl_iter_pw(istep, iter, conv_esolver, this->stp, this->kv, this->pw_wfc, *this->inp_);
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_PW", "after_scf");
    ModuleBase::timer::start("ESolver_KS_PW", "after_scf");

    // Calculate kinetic energy density tau for ELF if needed
    if (this->inp_->out_elf[0] > 0)
    {
        auto* elec_pw = static_cast<elecstate::ElecStatePW<T, Device>*>(this->pelec);
        auto& psi = *this->stp.template get_psi_t<T, Device>();
        elec_pw->cal_tau(psi);
    }

    ESolver_KS::after_scf(ucell, istep, conv_esolver);

    // Output quantities
    ModuleIO::ctrl_scf_pw<T, Device>(istep,
                                     ucell,
                                     this->pelec,
                                     this->chr,
                                     this->ppcell,
                                     this->kv,
                                     this->pw_wfc,
                                     this->pw_rho,
                                     this->pw_rhod,
                                     this->pw_big,
                                     this->stp,
                                     this->Pgrid,
                                     *this->inp_);

    const bool exchange_converged = !this->exx_helper->get_op_first_iter()
                                    && this->exx_helper->iteration_count() < this->inp_->exx_hybrid_step;
    if (conv_esolver && exchange_converged && this->inp_->calculation == "scf" && this->general_exx_info_.cal_exx
        && this->inp_->out_wfc_pw == 2 && this->inp_->out_freq_ion == 0 && this->inp_->kpar == 1 && this->inp_->bndpar == 1
        && this->inp_->nspin != 4 && this->inp_->symmetry != "1")
    {
        if (this->pw_wfc->poolrank == 0)
        {
            const std::string filename = PARAM.globalv.global_out_dir + "EXX_SOURCE";
            const std::string configuration = exx_restart_configuration(*this->inp_, this->general_exx_info_);
            std::ofstream out(filename);
            ModuleIO::write_exx_source(out, this->kv, this->pelec->wg, this->inp_->nspin,
                                      this->pw_wfc->nx, this->pw_wfc->ny, this->pw_wfc->nz, configuration);
        }
    }

    ModuleBase::timer::end("ESolver_KS_PW", "after_scf");
}

template <typename T, typename Device>
double ESolver_KS_PW<T, Device>::cal_energy()
{
    return this->pelec->f_en.etot;
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_force(BaseCell& basecell, ModuleBase::matrix& force)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    Forces<double, Device> ff(ucell.nat);

    // mohan add 2025-10-12
    this->stp.update_psi_d();

    // Calculate forces
    ff.cal_force(this->inp_->nspin, PARAM.globalv.domag, PARAM.globalv.domag_z, this->inp_->gga_grad,
                 this->inp_->dft_plus_u || this->inp_->sc_mag_switch,
                 ucell,
                 force,
                 this->get_vdw_result(),
                 *this->pelec,
                 this->pw_rhod,
                 &ucell.symm,
                 &this->sf,
                 this->solvent,
                 this->dftu_.get(),
                 &this->locpp,
                 &this->ppcell,
                 &this->kv,
                 this->pw_wfc,
                 this->stp.template get_psi_d<T, Device>());
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::cal_stress(BaseCell& basecell, ModuleBase::matrix& stress)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    Stress_PW<double, Device> ss(this->pelec);

    // mohan add 2025-10-12
    this->stp.update_psi_d();

    ss.cal_stress(this->inp_->nspin, PARAM.globalv.domag, PARAM.globalv.domag_z,
                       this->inp_->gga_grad, PARAM.globalv.gamma_only_pw,
                       stress,
                  ucell,
                  this->get_vdw_result(),
                  *this->dftu_,
                  this->locpp,
                  this->ppcell,
                  this->pw_rhod,
                  &ucell.symm,
                  &this->sf,
                  &this->kv,
                  this->pw_wfc,
                  this->general_exx_info_,
                  this->stp.template get_psi_d<T, Device>());

    // external stress
    double unit_transform = 0.0;
    unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    double external_stress[3] = {this->inp_->press1, this->inp_->press2, this->inp_->press3};
    for (int i = 0; i < 3; i++)
    {
        stress(i, i) -= external_stress[i] / unit_transform;
    }
}

template <typename T, typename Device>
void ESolver_KS_PW<T, Device>::after_all_runners(BaseCell& basecell)
{
    basecell.require_kind(BaseCell::Kind::unitcell, __FUNCTION__);
    UnitCell& ucell = static_cast<UnitCell&>(basecell);

    ESolver_KS::after_all_runners(ucell);

    ModuleIO::ctrl_runner_pw<T, Device>(ucell,
                                        this->pelec,
                                        this->pw_wfc,
                                        this->pw_rho,
                                        this->pw_rhod,
                                        this->chr,
                                        this->kv,
                                        this->stp,
                                        this->sf,
                                        this->ppcell,
                                        this->solvent,
                                        this->Pgrid,
                                        *this->inp_);

    elecstate::teardown_estate_pw(this->pelec, this->vsep_cell);
}

template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_CPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class ESolver_KS_PW<std::complex<float>, base_device::DEVICE_GPU>;
template class ESolver_KS_PW<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace ModuleESolver

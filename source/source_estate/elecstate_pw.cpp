#include "elecstate_pw.h"

#include "source_base/constants.h"
#include "source_base/libm/libm.h"
#include "source_base/module_device/device.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_estate/uspp_density.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/vnl_pw.h"

namespace elecstate
{

template <typename T, typename Device>
ElecStatePW<T, Device>::ElecStatePW(ModulePW::PW_Basis_K* wfc_basis_in,
                                    Charge* chr_in,
                                    K_Vectors* pkv_in,
                                    UnitCell* ucell_in,
                                    pseudopot_cell_vnl* ppcell_in,
                                    ModulePW::PW_Basis* rhopw_in,
                                    ModulePW::PW_Basis_Big* bigpw_in)
    : basis(wfc_basis_in)
{
    this->classname = "ElecStatePW";
    this->rhopw_smooth = rhopw_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->init_ks(chr_in, pkv_in, pkv_in->get_nks(), bigpw_in);
}

template <typename T, typename Device>
ElecStatePW<T, Device>::~ElecStatePW()
{
    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        delmem_var_op()(this->rho_data);
        delete[] this->rho;

        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            delmem_complex_op()(this->rhog_data);
            delete[] this->rhog;
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            delmem_var_op()(this->kin_r_data);
            delete[] this->kin_r;
        }
    }
    delmem_complex_op()(this->wfcr);
    delmem_complex_op()(this->wfcr_another_spin);
}

template <typename T, typename Device>
double ElecStatePW<T, Device>::get_spin_constrain_energy()
{
    spinconstrain::SpinConstrain<std::complex<double>>& sc = spinconstrain::SpinConstrain<std::complex<double>>::getScInstance();
    return sc.cal_escon();
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::init_rho_data()
{
    if (this->init_rho)
    {
        return;
    }

    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        this->rho = new Real*[this->charge->nspin];
        resmem_var_op()(this->rho_data, this->charge->nspin * this->charge->nrxx);
        for (int ii = 0; ii < this->charge->nspin; ii++)
        {
            this->rho[ii] = this->rho_data + ii * this->charge->nrxx;
        }
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            this->rhog = new T*[this->charge->nspin];
            resmem_complex_op()(this->rhog_data, this->charge->nspin * this->charge->rhopw->npw);
            for (int ii = 0; ii < this->charge->nspin; ii++)
            {
                this->rhog[ii] = this->rhog_data + ii * this->charge->rhopw->npw;
            }
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            this->kin_r = new Real*[this->charge->nspin];
            resmem_var_op()(this->kin_r_data, this->charge->nspin * this->charge->nrxx);
            for (int ii = 0; ii < this->charge->nspin; ii++)
            {
                this->kin_r[ii] = this->kin_r_data + ii * this->charge->nrxx;
            }
        }
    }
    else
    {
        this->rho = reinterpret_cast<Real**>(this->charge->rho);
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            this->rhog = reinterpret_cast<T**>(this->charge->rhog);
        }
        if (XC_Functional::get_func_type() == 3 || PARAM.inp.out_elf[0] > 0)
        {
            this->kin_r = reinterpret_cast<Real**>(this->charge->kin_r);
        }
    }
    resmem_complex_op()(this->wfcr, this->basis->nmaxgr, "ElecSPW::wfcr");
    resmem_complex_op()(this->wfcr_another_spin, this->basis->nrxx, "ElecSPW::wfcr_a");
    this->init_rho = true;
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::psiToRho(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "psiToRho");
    ModuleBase::timer::start("ElecStatePW", "psiToRho");

    this->init_rho_data();

    for (int is = 0; is < PARAM.inp.nspin; is++)
    {
        // denghui replaced at 20221110
        // ModuleBase::GlobalFunc::ZEROS(this->rho[is], this->charge->nrxx);
        setmem_var_op()(this->rho[is], 0, this->charge->nrxx);
        if (XC_Functional::get_func_type() == 3)
        {
            // ModuleBase::GlobalFunc::ZEROS(this->charge->kin_r[is], this->charge->nrxx);
            setmem_var_op()(this->kin_r[is], 0, this->charge->nrxx);
        }
        if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
        {
            setmem_complex_op()(this->rhog[is], 0, this->charge->rhopw->npw);
            std::fill(this->charge->rhog[is], this->charge->rhog[is] + this->charge->rhopw->npw, std::complex<double>(0, 0));
        }
    }

    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        this->updateRhoK(psi);
    }

    if (PARAM.inp.device == "gpu" || PARAM.inp.precision == "single")
    {
        for (int ii = 0; ii < PARAM.inp.nspin; ii++)
        {
            castmem_var_d2h_op()(this->charge->rho[ii], this->rho[ii], this->charge->nrxx);
            if (XC_Functional::get_func_type() == 3)
            {
                castmem_var_d2h_op()(this->charge->kin_r[ii], this->kin_r[ii], this->charge->nrxx);
            }
        }
    }

    this->add_usrho(psi);
    this->parallelK();
    ModuleBase::timer::end("ElecStatePW", "psiToRho");
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::updateRhoK(const psi::Psi<T, Device>& psi)
{
    this->rhoBandK(psi);
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::parallelK()
{
#ifdef __MPI
    this->charge->rho_mpi();
#endif
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::rhoBandK(const psi::Psi<T, Device>& psi)
{
    ModuleBase::TITLE("ElecStatePW", "rhoBandK");

    // moved by denghui to constructor at 20221110
    // used for plane wavefunction FFT3D to real space
    // static std::vector<T> wfcr;
    // wfcr.resize(this->basis->nmaxgr);
    // used for plane wavefunction FFT3D to real space, non-collinear spin case
    // static std::vector<std::complex<double>> wfcr_another_spin;
    // if (PARAM.inp.nspin == 4)
    //     wfcr_another_spin.resize(this->charge->nrxx);

    this->init_rho_data();
    int ik = psi.get_current_k();
    int npw = psi.get_current_ngk();
    int current_spin = 0;
    if (PARAM.inp.nspin == 2)
    {
        current_spin = this->klist->isk[ik];
    }
    int nbands = psi.get_nbands();
    //  here we compute the band energy: the sum of the eigenvalues
    if (PARAM.inp.nspin == 4)
    {
        int npwx = npw / 2;
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            /// be care of when smearing_sigma is large, wg would less than 0
            ///

            this->basis->recip_to_real(this->ctx, &psi(ibnd, 0), this->wfcr, ik);

            this->basis->recip_to_real(this->ctx, &psi(ibnd, npwx), this->wfcr_another_spin, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / ucell->omega);

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx,
                                  PARAM.globalv.domag,
                                  PARAM.globalv.domag_z,
                                  this->basis->nrxx,
                                  this->charge->nrxx,
                                  w1,
                                  this->rho,
                                  this->wfcr,
                                  this->wfcr_another_spin);
            }
        }
    }
    else
    {
        for (int ibnd = 0; ibnd < nbands; ibnd++)
        {
            ///
            /// only occupied band should be calculated.
            ///

            this->basis->recip_to_real(this->ctx, &psi(ibnd, 0), this->wfcr, ik);

            const auto w1 = static_cast<Real>(this->wg(ik, ibnd) / ucell->omega);

            if (w1 != 0.0)
            {
                // replaced by denghui at 20221110
                elecstate_pw_op()(this->ctx, current_spin, this->basis->nrxx, this->charge->nrxx, w1, this->rho, this->wfcr);
            }

            // kinetic energy density
            if (XC_Functional::get_func_type() == 3)
            {
                for (int j = 0; j < 3; j++)
                {
                    setmem_complex_op()(this->wfcr, 0, this->charge->nrxx);

                    meta_op()(this->ctx,
                              ik,
                              j,
                              npw,
                              this->basis->npwk_max,
                              static_cast<Real>(ucell->tpiba),
                              this->basis->template get_gcar_data<Real>(),
                              this->basis->template get_kvec_c_data<Real>(),
                              &psi(ibnd, 0),
                              this->wfcr);

                    this->basis->recip_to_real(this->ctx, this->wfcr, this->wfcr, ik);

                    elecstate_pw_op()(this->ctx, current_spin, this->charge->nrxx, this->charge->nrxx, w1, this->kin_r, this->wfcr);
                }
            }
        }
    }
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::cal_becsum(const psi::Psi<T, Device>& psi)
{
    this->becsum_.assign(uspp_becsum_size(*ucell, *ppcell, PARAM.inp.nspin), 0.0);
    UsppProjector<T, Device> projector(*ucell, *ppcell, psi.get_nbands());
    std::vector<double> weights(psi.get_nbands());
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        psi.fix_k(ik);
        // SCF supplies occupations including k weights; output callers supply their own state weights.
        for (int ib = 0; ib < psi.get_nbands(); ++ib)
        {
            weights[ib] = this->wg(ik, ib);
        }
        projector.accumulate(ik, psi.get_pointer(), psi.get_nbasis(), psi.get_current_ngk(), this->klist->isk[ik], weights, &this->becsum_);
    }
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::add_usrho(const psi::Psi<T, Device>& psi)
{
    if (PARAM.globalv.use_uspp)
    {
        this->cal_becsum(psi);
    }

    // transform soft charge to recip space using smooth grids
    if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            this->rhopw_smooth->real2recip(this->charge->rho[is], this->charge->rhog[is]);
        }
    }

    // \sum_lm Q_lm(r) \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
    // add to the charge density in reciprocal space the part which is due to the US augmentation.
    if (PARAM.globalv.use_uspp)
    {
        this->addusdens_g(this->charge->rhog);
    }
    // transform back to real space using dense grids
    if (PARAM.globalv.double_grid || PARAM.globalv.use_uspp)
    {
        for (int is = 0; is < PARAM.inp.nspin; is++)
        {
            this->charge->rhopw->recip2real(this->charge->rhog[is], this->charge->rho[is]);
        }
    }
}

template <typename T, typename Device>
void ElecStatePW<T, Device>::addusdens_g(std::complex<double>** rhog)
{
    add_uspp_density(*ucell, *ppcell, *this->charge->rhopw, PARAM.inp.nspin, this->becsum_, rhog);
}

// Taoni add 2026-09-02
// Added to fix USPP single force/stress reading the former float becsum as double.
// The double-only drivers receive a base ElecState, while becsum belongs to the precision-templated ElecStatePW.
// Refactor this bridge for true float force/stress.
template <typename Device>
const std::vector<double>* get_becsum(const ElecState& elec)
{
    const ElecStatePW<std::complex<double>, Device>* double_elec = dynamic_cast<const ElecStatePW<std::complex<double>, Device>*>(&elec);
    if (double_elec != nullptr)
    {
        return &double_elec->get_becsum();
    }

    const ElecStatePW<std::complex<float>, Device>* single_elec = dynamic_cast<const ElecStatePW<std::complex<float>, Device>*>(&elec);
    if (single_elec != nullptr)
    {
        return &single_elec->get_becsum();
    }

    return nullptr;
}

template class ElecStatePW<std::complex<float>, base_device::DEVICE_CPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_CPU>;
template const std::vector<double>* get_becsum<base_device::DEVICE_CPU>(const ElecState& elec);
#if ((defined __CUDA) || (defined __ROCM))
template class ElecStatePW<std::complex<float>, base_device::DEVICE_GPU>;
template class ElecStatePW<std::complex<double>, base_device::DEVICE_GPU>;
template const std::vector<double>* get_becsum<base_device::DEVICE_GPU>(const ElecState& elec);
#endif

} // namespace elecstate

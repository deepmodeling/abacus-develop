#include "hamilt_pw.h"

#include "op_pw_ekin.h"
#include "op_pw_ekin_td.h"
#include "op_pw_exx.h"
#include "op_pw_meta.h"
#include "op_pw_nl.h"
#include "op_pw_nl_td.h"
#include "op_pw_proj.h"
#include "op_pw_veff.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_xc/general_exx_info.h" // for General_Exx_Info type
#include "source_io/module_parameter/parameter.h"
#include "source_estate/module_pot/pot_pcc.h"

namespace hamilt
{

template <typename T, typename Device>
HamiltPW<T, Device>::HamiltPW(elecstate::Potential* pot_in,
                              ModulePW::PW_Basis_K* wfc_basis,
                              K_Vectors* pkv,
                              pseudopot_cell_vnl* nlpp,
                              Plus_U_Base* p_dftu,
                              const UnitCell* ucell,
                              const General_Exx_Info* exx_info)
    : ucell(ucell)
{
    this->classname = "HamiltPW";
    this->ppcell = nlpp;
    this->qq_nt = this->ppcell->template get_qq_nt_data<Real>();
    this->qq_so = this->ppcell->template get_qq_so_data<Real>();
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    const auto tpiba2 = static_cast<Real>(ucell->tpiba2);
    const auto tpiba = static_cast<Real>(ucell->tpiba);
    const int* isk = pkv->isk.data();
    const Real* gk2 = wfc_basis->get_gk2_data<Real>();
    const bool is_tddft = PARAM.inp.esolver_type == "tddft";
    const bool velocity_gauge = is_tddft && PARAM.inp.td_stype == 1;

    if (PARAM.inp.t_in_h)
    {
        // Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>
        Operator<T, Device>* ekinetic
            = new Ekinetic<OperatorPW<T, Device>>(tpiba2, gk2, wfc_basis->nks, wfc_basis->npwk_max);
        if (this->ops == nullptr)
        {
            this->ops = ekinetic;
        }
        else
        {
            this->ops->add(ekinetic);
        }

        // Add the velocity-gauge kinetic correction 2*A*p + A^2.
        if (velocity_gauge)
        {
            Operator<T, Device>* td_ekinetic = new TDEkineticPW<T, Device>(wfc_basis);

            if (this->ops == nullptr)
            {
                this->ops = td_ekinetic;
            }
            else
            {
                this->ops->add(td_ekinetic);
            }
        }
    }
    const Input_para& input = PARAM.inp;
    if (input.vl_in_h)
    {
        std::vector<std::string> pot_register_in;
        if (input.vion_in_h)
        {
            pot_register_in.push_back("local");
        }
        if (input.vh_in_h)
        {
            pot_register_in.push_back("hartree");
        }
        // no variable can choose xc, maybe it is necessary
        pot_register_in.push_back("xc");
        if (input.imp_sol == 1)
        {
            pot_register_in.push_back("surchem");
        }
        if (input.imp_sol == 2) { pot_register_in.push_back("sccs"); }
        if (input.assume_isolated == "pcc_0d")
        {
            pot_register_in.push_back("pcc");
        }
        if (input.assume_isolated == "pcc_2d")
        {
            const int count = pkv->get_nks();
            elecstate::PotPcc::validate_kpoints(pkv->kvec_d, count, input.pcc_2d_axis);
            pot_register_in.push_back("pcc_2d");
        }
        if (input.efield_flag)
        {
            pot_register_in.push_back("efield");
        }
        if (input.gate_flag)
        {
            pot_register_in.push_back("gatefield");
        }
        if (is_tddft)
        {
            pot_register_in.push_back("tddft");
        }
        if (PARAM.inp.ml_exx) // sunliang
        {
            pot_register_in.push_back("ml_exx");
        }
        // DFT-1/2
        if (PARAM.inp.dfthalf_type == 1)
        {
            pot_register_in.push_back("dfthalf");
        }
        // only Potential is not empty, Veff and Meta are available
        if (pot_register_in.size() > 0)
        {
            // register Potential by gathered operator
            pot_in->pot_register(pot_register_in, input);
            Operator<T, Device>* veff = new Veff<OperatorPW<T, Device>>(isk,
                                                                        pot_in->get_veff_smooth_data<Real>(),
                                                                        pot_in->get_veff_smooth().nr,
                                                                        pot_in->get_veff_smooth().nc,
                                                                        wfc_basis);
            if (this->ops == nullptr)
            {
                this->ops = veff;
            }
            else
            {
                this->ops->add(veff);
            }
            Operator<T, Device>* meta = new Meta<OperatorPW<T, Device>>(tpiba,
                                                                        isk,
                                                                        pot_in->get_vofk_smooth_data<Real>(),
                                                                        pot_in->get_vofk_smooth().nr,
                                                                        pot_in->get_vofk_smooth().nc,
                                                                        wfc_basis);
            this->ops->add(meta);
        }
    }
    if (PARAM.inp.vnl_in_h)
    {
        Operator<T, Device>* nonlocal = nullptr;

        // Use vector-potential-shifted projectors in the velocity gauge.
        if (velocity_gauge)
        {
            nonlocal = new TDNonlocalPW<T, Device>(isk, this->ppcell, ucell, wfc_basis);
        }
        else
        {
            // Ground-state and length-gauge calculations use the standard operator.
            nonlocal = new Nonlocal<OperatorPW<T, Device>>(isk, this->ppcell, ucell, wfc_basis);
        }

        if (this->ops == nullptr)
        {
            this->ops = nonlocal;
        }
        else
        {
            this->ops->add(nonlocal);
        }
    }
    if (PARAM.inp.sc_mag_switch || PARAM.inp.dft_plus_u)
    {
        Operator<T, Device>* onsite_proj = new OnsiteProj<OperatorPW<T, Device>>(isk,
                                                                                 ucell,
                                                                                 p_dftu,
                                                                                 PARAM.inp.sc_mag_switch,
                                                                                 (PARAM.inp.dft_plus_u > 0));
        this->ops->add(onsite_proj);
    }
    if (exx_info && exx_info->cal_exx)
    {
        auto exx = new OperatorEXXPW<T, Device>(isk,
                                                wfc_basis,
                                                pot_in->get_rho_basis(),
                                                pkv,
                                                ucell,
                                                *exx_info,
                                                PARAM.inp.nspin,
                                                GlobalV::KPAR,
                                                GlobalV::MY_RANK,
                                                GlobalV::MY_POOL);
        if (this->ops == nullptr)
        {
            this->ops = exx;
        }
        else
        {
            this->ops->add(exx);
            // exx->set_psi(&this->psi);
        }
    }
    return;
}

template <typename T, typename Device>
HamiltPW<T, Device>::~HamiltPW()
{
    if (this->ops != nullptr)
    {
        delete this->ops;
    }
}

template <typename T, typename Device>
void HamiltPW<T, Device>::bind_td_state(const Real* veff, const Real* vofk, const ModuleBase::Vector3<double>& A_ha)
{
    Operator<T, Device>* op = this->ops;
    while (op != nullptr)
    {
        if (op->get_cal_type() == calculation_type::pw_veff)
        {
            static_cast<Veff<OperatorPW<T, Device>>*>(op)->set_veff(veff);
        }
        else if (op->get_cal_type() == calculation_type::pw_meta)
        {
            static_cast<Meta<OperatorPW<T, Device>>*>(op)->set_vk(vofk);
        }
        else if (op->get_cal_type() == calculation_type::pw_ekinetic_td)
        {
            static_cast<TDEkineticPW<T, Device>*>(op)->set_A_ha(A_ha);
        }
        else if (op->get_cal_type() == calculation_type::pw_nonlocal)
        {
            TDNonlocalPW<T, Device>* td_nonlocal = dynamic_cast<TDNonlocalPW<T, Device>*>(op);
            if (td_nonlocal != nullptr) td_nonlocal->set_A_ha(A_ha);
        }
        op = op->next_op;
    }
}

template <typename T, typename Device>
void HamiltPW<T, Device>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW", "updateHk");
    this->ops->init(ik);
    ModuleBase::TITLE("HamiltPW", "updateHk");
}

// This routine applies the S matrix to m wavefunctions psi and puts
// the results in spsi.
// Requires the products of psi with all beta functions in array
// becp(nkb,m) (maybe calculated in hPsi).
template <typename T, typename Device>
void HamiltPW<T, Device>::sPsi(const T* psi_in, // psi
                               T* spsi,         // spsi
                               const int nrow,  // dimension of spsi: nbands * nrow
                               const int npw,   // number of plane waves
                               const int nbands // number of bands
) const
{
    ModuleBase::TITLE("HamiltPW", "sPsi");

    const T one{1, 0};
    const T zero{0, 0};

    syncmem_op()(spsi, psi_in, static_cast<size_t>(nbands * nrow));
    if (PARAM.globalv.use_uspp)
    {
        T* becp = nullptr;
        T* ps = nullptr;
        // psi updated, thus update <beta|psi>
        if (this->ppcell->nkb > 0)
        {
            resmem_complex_op()(becp, nbands * this->ppcell->nkb, "Hamilt<PW>::becp");
            char transa = 'C';
            char transb = 'N';
            if (nbands == 1)
            {
                int inc = 1;
                gemv_op()(transa,
                          npw,
                          this->ppcell->nkb,
                          &one,
                          this->vkb,
                          this->ppcell->vkbnc,
                          psi_in,
                          inc,
                          &zero,
                          becp,
                          inc);
            }
            else
            {
                gemm_op()(transa,
                          transb,
                          this->ppcell->nkb,
                          nbands,
                          npw,
                          &one,
                          this->vkb,
                          this->ppcell->vkbnc,
                          psi_in,
                          nrow,
                          &zero,
                          becp,
                          this->ppcell->nkb);
            }

            Parallel_Reduce::reduce_pool(becp, this->ppcell->nkb * nbands);
        }

        resmem_complex_op()(ps, this->ppcell->nkb * nbands, "Hamilt<PW>::ps");
        setmem_complex_op()(ps, 0, this->ppcell->nkb * nbands);

        // spsi = psi + sum qq <beta|psi> |beta>
        // qq <beta|psi>
        char transa = 'N';
        char transb = 'N';
        for (int it = 0; it < ucell->ntype; it++)
        {
            Atom* atoms = &ucell->atoms[it];
            if (atoms->ncpp.tvanp)
            {
                const int nh = atoms->ncpp.nh;
                T* qqc = nullptr;
                resmem_complex_op()(qqc, nh * nh, "Hamilt<PW>::qqc");
                std::vector<T> qqc_host(nh*nh);
                const double* qq_now_host = &this->ppcell->qq_nt.ptr[it * this->ppcell->nhm * this->ppcell->nhm];

                for (int i = 0; i < nh; i++)
                {
                    for (int j = 0; j < nh; j++)
                    {
                        const int source_index = i * this->ppcell->nhm + j;
                        const int target_index = i * nh + j;

                        qqc_host[target_index] = static_cast<Real>(qq_now_host[source_index]) * one;
                    }
                }

                syncmem_complex_h2d_op()(qqc, qqc_host.data(), qqc_host.size());

                for (int ia = 0; ia < atoms->na; ia++)
                {
                    const int iat = ucell->itia2iat(it, ia);
                    gemm_op()(transa,
                              transb,
                              nh,
                              nbands,
                              nh,
                              &one,
                              qqc,
                              nh,
                              &becp[this->ppcell->indv_ijkb0[iat]],
                              this->ppcell->nkb,
                              &zero,
                              &ps[this->ppcell->indv_ijkb0[iat]],
                              this->ppcell->nkb);
                }
                delmem_complex_op()(qqc);
            }
        }

        if (nbands == 1)
        {
            const int inc = 1;
            gemv_op()(transa,
                      npw,
                      this->ppcell->nkb,
                      &one,
                      this->vkb,
                      this->ppcell->vkbnc,
                      ps,
                      inc,
                      &one,
                      spsi,
                      inc);
        }
        else
        {
            gemm_op()(transa,
                      transb,
                      npw,
                      nbands,
                      this->ppcell->nkb,
                      &one,
                      this->vkb,
                      this->ppcell->vkbnc,
                      ps,
                      this->ppcell->nkb,
                      &one,
                      spsi,
                      nrow);
        }
        delmem_complex_op()(ps);
        delmem_complex_op()(becp);
    }
}

template <typename T, typename Device>
void HamiltPW<T, Device>::set_exx_helper(Exx_Helper<T, Device>& exx_helper)
{
    auto op = this->ops;
    while (op != nullptr)
    {
        if (op->get_cal_type() == calculation_type::pw_exx)
        {
            exx_helper.op_exx = reinterpret_cast<OperatorEXXPW<T, Device>*>(op);
            exx_helper.set_op();
        }
        op = op->next_op;
    }
}

template class HamiltPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HamiltPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HamiltPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HamiltPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hamilt

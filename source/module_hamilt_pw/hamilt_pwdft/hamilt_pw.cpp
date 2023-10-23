#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

#include "operator_pw/veff_pw.h"
#include "operator_pw/ekinetic_pw.h"
#include "operator_pw/meta_pw.h"
#include "operator_pw/nonlocal_pw.h"

namespace hamilt
{

template<typename T, typename Device>
HamiltPW<T, Device>::HamiltPW(elecstate::Potential* pot_in, ModulePW::PW_Basis_K* wfc_basis, K_Vectors* pkv)
{
    this->classname = "HamiltPW";
    this->ppcell = &GlobalC::ppcell;
    this->qq_nt = this->ppcell->template get_qq_nt_data<Real>();
    this->qq_so = this->ppcell->template get_qq_so_data<Real>();
    this->vkb = this->ppcell->template get_vkb_data<Real>();
    const auto tpiba2 = static_cast<Real>(GlobalC::ucell.tpiba2);
    const auto tpiba = static_cast<Real>(GlobalC::ucell.tpiba);
    const int* isk = pkv->isk.data();
    const Real* gk2 = wfc_basis->get_gk2_data<Real>();

    if (GlobalV::T_IN_H)
    {
        // Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>
        Operator<T, Device>* ekinetic
            = new Ekinetic<OperatorPW<T, Device>>(tpiba2, gk2, wfc_basis->nks, wfc_basis->npwk_max);
        if(this->ops == nullptr)
        {
            this->ops = ekinetic;
        }
        else
        {
            this->ops->add(ekinetic);
        }
    }
    if (GlobalV::VL_IN_H)
    {
        std::vector<std::string> pot_register_in;
        if (GlobalV::VION_IN_H)
        {
            pot_register_in.push_back("local");
        }
        if (GlobalV::VH_IN_H)
        {
            pot_register_in.push_back("hartree");
        }
        //no variable can choose xc, maybe it is necessary
        pot_register_in.push_back("xc");
        if (GlobalV::imp_sol)
        {
            pot_register_in.push_back("surchem");
        }
        if (GlobalV::EFIELD_FLAG)
        {
            pot_register_in.push_back("efield");
        }
        if (GlobalV::GATE_FLAG)
        {
            pot_register_in.push_back("gatefield");
        }
        //only Potential is not empty, Veff and Meta are available
        if(pot_register_in.size()>0)
        {
            //register Potential by gathered operator
            pot_in->pot_register(pot_register_in);
            Operator<T, Device>* veff = new Veff<OperatorPW<T, Device>>(isk,
                                                                        pot_in->get_veff_smooth_data<Real>(),
                                                                        pot_in->get_veff_smooth().nr,
                                                                        pot_in->get_veff_smooth().nc,
                                                                        wfc_basis);
            if(this->ops == nullptr)
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
    if (GlobalV::VNL_IN_H)
    {
        Operator<T, Device>* nonlocal
            = new Nonlocal<OperatorPW<T, Device>>(isk, &GlobalC::ppcell, &GlobalC::ucell, wfc_basis);
        if(this->ops == nullptr)
        {
            this->ops = nonlocal;
        }
        else
        {
            this->ops->add(nonlocal);
        }
    }
    return;
}

template<typename T, typename Device>
HamiltPW<T, Device>::~HamiltPW()
{
    if(this->ops!= nullptr)
    {
        delete this->ops;
    }
}

template<typename T, typename Device>
void HamiltPW<T, Device>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");
    this->ops->init(ik);
    ModuleBase::TITLE("HamiltPW","updateHk");
}

template<typename T, typename Device>
template<typename T_in, typename Device_in>
HamiltPW<T, Device>::HamiltPW(const HamiltPW<T_in, Device_in> *hamilt)
{
    this->classname = hamilt->classname;
    this->ppcell = hamilt->ppcell;
    this->qq_nt = hamilt->qq_nt;
    this->qq_so = hamilt->qq_so;
    this->vkb = hamilt->vkb;
    this->becp = hamilt->becp;
    OperatorPW<std::complex<T_in>, Device_in> * node =
            reinterpret_cast<OperatorPW<std::complex<T_in>, Device_in> *>(hamilt->ops);

    while(node != nullptr) {
        if (node->classname == "Ekinetic") {
            Operator<T, Device>* ekinetic =
                    new Ekinetic<OperatorPW<T, Device>>(
                            reinterpret_cast<const Ekinetic<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = ekinetic;
            }
            else {
                this->ops->add(ekinetic);
            }
            // this->ops = reinterpret_cast<Operator<T, Device>*>(node);
        }
        else if (node->classname == "Nonlocal") {
            Operator<T, Device>* nonlocal =
                    new Nonlocal<OperatorPW<T, Device>>(
                            reinterpret_cast<const Nonlocal<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = nonlocal;
            }
            else {
                this->ops->add(nonlocal);
            }
        }
        else if (node->classname == "Veff") {
            Operator<T, Device>* veff =
                    new Veff<OperatorPW<T, Device>>(
                            reinterpret_cast<const Veff<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = veff;
            }
            else {
                this->ops->add(veff);
            }
        }
        else if (node->classname == "Meta") {
            Operator<T, Device>* meta =
                    new Meta<OperatorPW<T, Device>>(
                            reinterpret_cast<const Meta<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = meta;
            }
            else {
                this->ops->add(meta);
            }
        }
        else {
            ModuleBase::WARNING_QUIT("HamiltPW", "Unrecognized Operator type!");
        }
        node = reinterpret_cast<OperatorPW<std::complex<T_in>, Device_in> *>(node->next_op);
    }
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

    syncmem_op()(this->ctx, this->ctx, spsi, psi_in, static_cast<size_t>(nbands * nrow));
    if (GlobalV::use_uspp)
    {
        // psi updated, thus update <beta|psi>
        if (this->ppcell->nkb > 0)
        {
            resmem_complex_op()(this->ctx, this->becp, nbands * this->ppcell->nkb, "Hamilt<PW>::becp");
            char transa = 'C';
            char transb = 'N';
            if (nbands == 1)
            {
                int inc = 1;
                gemv_op()(this->ctx,
                          transa,
                          npw,
                          this->ppcell->nkb,
                          &this->one,
                          this->vkb,
                          this->ppcell->vkb.nc,
                          psi_in,
                          inc,
                          &this->zero,
                          this->becp,
                          inc);
            }
            else
            {
                gemm_op()(this->ctx,
                          transa,
                          transb,
                          this->ppcell->nkb,
                          nbands,
                          npw,
                          &this->one,
                          this->vkb,
                          this->ppcell->vkb.nc,
                          psi_in,
                          nrow,
                          &this->zero,
                          this->becp,
                          this->ppcell->nkb);
            }

            Parallel_Reduce::reduce_pool(this->becp, this->ppcell->nkb * nbands);
        }

        resmem_complex_op()(this->ctx, ps, this->ppcell->nkb * nbands, "Hamilt<PW>::ps");
        setmem_complex_op()(this->ctx, ps, 0, this->ppcell->nkb * nbands);

        // spsi = psi + sum qq <beta|psi> |beta>
        int iat = 0;
        if (GlobalV::NONCOLIN)
        {
            // spsi_nc
            std::cout << " noncolinear in uspp is not implemented yet " << std::endl;
            exit(0);
        }
        else
        {
            // qq <beta|psi>
            char transa = 'N';
            char transb = 'N';
            for (int it = 0; it < GlobalC::ucell.ntype; it++)
            {
                Atom* atoms = &GlobalC::ucell.atoms[it];
                if (atoms->ncpp.tvanp)
                {
                    const int nh = atoms->ncpp.nh;
                    T* qqc = nullptr;
                    resmem_complex_op()(this->ctx, qqc, nh * nh, "Hamilt<PW>::qqc");
                    Real* qq_now = &qq_nt[it * this->ppcell->nhm * this->ppcell->nhm];
                    for (int i = 0; i < nh; i++)
                    {
                        for (int j = 0; j < nh; j++)
                        {
                            int index = i * this->ppcell->nhm + j;
                            qqc[i * nh + j] = qq_now[index] * this->one;
                        }
                    }
                    for (int ia = 0; ia < atoms->na; ia++)
                    {
                        gemm_op()(this->ctx,
                                  transa,
                                  transb,
                                  nh,
                                  nbands,
                                  nh,
                                  &this->one,
                                  qqc,
                                  nh,
                                  &becp[this->ppcell->indv_ijkb0[iat]],
                                  this->ppcell->nkb,
                                  &this->zero,
                                  &ps[this->ppcell->indv_ijkb0[iat]],
                                  this->ppcell->nkb);
                        iat++;
                    }
                    delmem_complex_op()(ctx, qqc);
                }
                else
                {
                    iat += atoms->na;
                }
            }

            if (nbands == 1)
            {
                const int inc = 1;
                gemv_op()(this->ctx,
                          transa,
                          npw,
                          this->ppcell->nkb,
                          &this->one,
                          this->vkb,
                          this->ppcell->vkb.nc,
                          ps,
                          inc,
                          &this->one,
                          spsi,
                          inc);
            }
            else
            {
                gemm_op()(this->ctx,
                          transa,
                          transb,
                          npw,
                          nbands,
                          this->ppcell->nkb,
                          &this->one,
                          this->vkb,
                          this->ppcell->vkb.nc,
                          ps,
                          this->ppcell->nkb,
                          &this->one,
                          spsi,
                          nrow);
            }
        }
    }

    delmem_complex_op()(ctx, ps);
    delmem_complex_op()(ctx, becp);

    ModuleBase::TITLE("HamiltPW", "sPsi");
}

template class HamiltPW<std::complex<float>, psi::DEVICE_CPU>;
template class HamiltPW<std::complex<double>, psi::DEVICE_CPU>;
// template HamiltPW<std::complex<double>, psi::DEVICE_CPU>::HamiltPW(const HamiltPW<std::complex<double>, psi::DEVICE_CPU> *hamilt);
#if ((defined __CUDA) || (defined __ROCM))
template class HamiltPW<std::complex<float>, psi::DEVICE_GPU>;
template class HamiltPW<std::complex<double>, psi::DEVICE_GPU>;
// template HamiltPW<std::complex<double>, psi::DEVICE_GPU>::HamiltPW(const HamiltPW<std::complex<double>, psi::DEVICE_GPU> *hamilt);
#endif

} // namespace hamilt
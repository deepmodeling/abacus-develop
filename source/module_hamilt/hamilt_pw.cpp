#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "src_pw/global.h"

#include "ks_pw/veff_pw.h"
#include "ks_pw/ekinetic_pw.h"
#include "ks_pw/meta_pw.h"
#include "ks_pw/nonlocal_pw.h"

namespace hamilt
{

template<typename FPTYPE, typename Device>
HamiltPW<FPTYPE, Device>::HamiltPW()
{
    this->classname = "HamiltPW";
    const double tpiba2 = GlobalC::ucell.tpiba2;
    const double tpiba = GlobalC::ucell.tpiba;
    const int* isk = GlobalC::kv.isk.data();
    const double* gk2 = GlobalC::wfcpw->gk2;

    if (GlobalV::T_IN_H)
    {
        // Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>
        Operator<std::complex<FPTYPE>, Device>* ekinetic = new Ekinetic<OperatorPW<FPTYPE, Device>>(
            tpiba2, 
            gk2,
            GlobalC::wfcpw->nks,
            GlobalC::wfcpw->npwk_max
        );
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
        Operator<std::complex<FPTYPE>, Device>* veff = new Veff<OperatorPW<FPTYPE, Device>>(
            isk,
            &(GlobalC::pot.vr_eff),
            GlobalC::wfcpw
        );
        if(this->ops == nullptr)
        {
            this->ops = veff;
        }
        else
        {
            this->ops->add(veff);
        }
    }
    if (GlobalV::VNL_IN_H)
    {
        Operator<std::complex<FPTYPE>, Device>* nonlocal = new Nonlocal<OperatorPW<FPTYPE, Device>>(
            isk,
            &GlobalC::ppcell,
            &GlobalC::ucell
        );
        if(this->ops == nullptr)
        {
            this->ops = nonlocal;
        }
        else
        {
            this->ops->add(nonlocal);
        }
    }
    Operator<std::complex<FPTYPE>, Device>* meta = new Meta<OperatorPW<FPTYPE, Device>>(
        tpiba,
        isk,
        &GlobalC::pot.vofk,
        GlobalC::wfcpw
    );
    if(this->ops == nullptr)
    {
        this->ops = meta;
    }
    else
    {
        this->ops->add(meta);
    }
}

template<typename FPTYPE, typename Device>
HamiltPW<FPTYPE, Device>::~HamiltPW()
{
    if(this->ops!= nullptr)
    {
        delete this->ops;
    }
}

template<typename FPTYPE, typename Device>
void HamiltPW<FPTYPE, Device>::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");
    this->ops->init(ik);
    ModuleBase::TITLE("HamiltPW","updateHk");
}

template<typename FPTYPE, typename Device>
void HamiltPW<FPTYPE, Device>::sPsi
(
    const std::complex<double> *psi,
    std::complex<double> *spsi,
    size_t size
) const
{
    // ModuleBase::GlobalFunc::COPYARRAY(psi, spsi, size);
    // denghui replaced at 2022.11.04
    syncmem_complex_op()(this->ctx, this->ctx, spsi, psi, size);
}

template<typename FPTYPE, typename Device>
template<typename T_in, typename Device_in>
HamiltPW<FPTYPE, Device>::HamiltPW(const HamiltPW<T_in, Device_in> *hamilt)
{
    this->classname = hamilt->classname;
    OperatorPW<std::complex<T_in>, Device_in> * node =
            reinterpret_cast<OperatorPW<std::complex<T_in>, Device_in> *>(hamilt->ops);

    while(node != nullptr) {
        if (node->classname == "Ekinetic") {
            Operator<std::complex<FPTYPE>, Device>* ekinetic =
                    new Ekinetic<OperatorPW<FPTYPE, Device>>(
                            reinterpret_cast<const Ekinetic<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = ekinetic;
            }
            else {
                this->ops->add(ekinetic);
            }
            // this->ops = reinterpret_cast<Operator<std::complex<FPTYPE>, Device>*>(node);
        }
        else if (node->classname == "Nonlocal") {
            Operator<std::complex<FPTYPE>, Device>* nonlocal =
                    new Nonlocal<OperatorPW<FPTYPE, Device>>(
                            reinterpret_cast<const Nonlocal<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = nonlocal;
            }
            else {
                this->ops->add(nonlocal);
            }
        }
        else if (node->classname == "Veff") {
            Operator<std::complex<FPTYPE>, Device>* veff =
                    new Veff<OperatorPW<FPTYPE, Device>>(
                            reinterpret_cast<const Veff<OperatorPW<T_in, Device_in>>*>(node));
            if(this->ops == nullptr) {
                this->ops = veff;
            }
            else {
                this->ops->add(veff);
            }
        }
        else if (node->classname == "Meta") {
            Operator<std::complex<FPTYPE>, Device>* meta =
                    new Meta<OperatorPW<FPTYPE, Device>>(
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

template class HamiltPW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HamiltPW<double, psi::DEVICE_GPU>;
template HamiltPW<double, psi::DEVICE_CPU>::HamiltPW(const HamiltPW<double, psi::DEVICE_GPU> *hamilt);
template HamiltPW<double, psi::DEVICE_GPU>::HamiltPW(const HamiltPW<double, psi::DEVICE_CPU> *hamilt);
#endif

} // namespace hamilt
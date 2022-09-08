#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"

#include "ks_pw/veff_pw.h"
#include "ks_pw/ekinetic_pw.h"
#include "ks_pw/meta_pw.h"
#include "ks_pw/nonlocal_pw.h"

namespace hamilt
{

HamiltPW::HamiltPW()
{
    this->classname = "HamiltPW";
    const double tpiba2 = GlobalC::ucell.tpiba2;
    const double tpiba = GlobalC::ucell.tpiba;
    const int* isk = GlobalC::kv.isk.data();
    const double* gk2 = GlobalC::wfcpw->gk2;

    if (GlobalV::T_IN_H)
    {
        Operator<std::complex<double>>* ekinetic = new Ekinetic<OperatorPW>( 
            tpiba2, 
            gk2,
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
        Operator<std::complex<double>>* veff = new Veff<OperatorPW>(
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
        Operator<std::complex<double>>* nonlocal = new Nonlocal<OperatorPW>(
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
    Operator<std::complex<double>>* meta = new Meta<OperatorPW>(
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

HamiltPW::~HamiltPW()
{
    int index = 0;
    delete this->ops;
}

void HamiltPW::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");

    this->ops->init(ik);

    return;
}

void HamiltPW::sPsi
(
    const std::complex<double> *psi,
    std::complex<double> *spsi,
    size_t size
) const
{
    ModuleBase::GlobalFunc::COPYARRAY(psi, spsi, size);
    return;
}

} // namespace hamilt
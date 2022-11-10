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

HamiltPW::HamiltPW(elecstate::Potential* pot_in)
{
    this->classname = "HamiltPW";
    const double tpiba2 = GlobalC::ucell.tpiba2;
    const double tpiba = GlobalC::ucell.tpiba;
    const int* isk = GlobalC::kv.isk.data();
    const double* gk2 = GlobalC::wfcpw->gk2;

    if (GlobalV::T_IN_H)
    {
        // Operator<double>* ekinetic = new Ekinetic<OperatorLCAO<double>>
        Operator<std::complex<double>>* ekinetic = new Ekinetic<OperatorPW<double>>(
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
        std::vector<string> pot_register_in;
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
            pot_register_in.push_back("gate");
        }
        //only Potential is not empty, Veff and Meta are available
        if(pot_register_in.size()>0)
        {
            //register Potential by gathered operator
            pot_in->pot_register(pot_register_in);
            Operator<std::complex<double>>* veff = new Veff<OperatorPW<double>>(
                isk,
                &(pot_in->get_effective_v()),
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
            Operator<std::complex<double>>* meta = new Meta<OperatorPW<double>>(
                tpiba,
                isk,
                &(pot_in->get_effective_vofk()),
                GlobalC::wfcpw
            );
            this->ops->add(meta);
        }
    }
    if (GlobalV::VNL_IN_H)
    {
        Operator<std::complex<double>>* nonlocal = new Nonlocal<OperatorPW<double>>(
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
    return;
}

HamiltPW::~HamiltPW()
{
    if(this->ops!= nullptr)
    {
        delete this->ops;
    }
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
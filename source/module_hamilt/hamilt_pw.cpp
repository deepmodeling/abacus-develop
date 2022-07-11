#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"

#include "veff_pw.h"
#include "ekinetic_pw.h"
#include "meta_pw.h"
#include "nonlocal_pw.h"

namespace hamilt
{

HamiltPW::HamiltPW()
{
    this->classname = "HamiltPW";
    this->ops.resize(0);
    const int npwx = GlobalC::wf.npwx;
    const int npol = GlobalV::NPOL;
    const double tpiba2 = GlobalC::ucell.tpiba2;
    const int* ngk = GlobalC::kv.ngk.data();
    const int* isk = GlobalC::kv.isk.data();
    const double* gk2 = GlobalC::wfcpw->gk2;

    if (GlobalV::T_IN_H)
    {
        Operator* ekinetic = new EkineticPW(
            npwx, 
            npol, 
            tpiba2, 
            ngk,
            gk2
        );
        this->ops.push_back(ekinetic);
    }
    if (GlobalV::VL_IN_H)
    {
        Operator* veff = new VeffPW(
            npwx,
            npol,
            ngk,
            isk,
            &(GlobalC::pot.vr_eff),
            GlobalC::wfcpw
        );
        this->ops.push_back(veff);
    }
    if (GlobalV::VNL_IN_H)
    {
        Operator* nonlocal = new NonlocalPW(
            npwx,
            npol,
            ngk,
            isk,
            &GlobalC::ppcell,
            &GlobalC::ucell
        );
        this->ops.push_back(nonlocal);
    }
    Operator* meta = new MetaPW(
        npwx,
        npol,
        tpiba2,
        ngk,
        isk,
        &GlobalC::pot.vofk,
        GlobalC::wfcpw
    );
    this->ops.push_back(meta);
}

HamiltPW::~HamiltPW()
{
    int index = 0;
    if (GlobalV::T_IN_H)
    {
        delete (EkineticPW*)this->ops[index++];
    }
    if (GlobalV::VL_IN_H)
    {
        delete (VeffPW*)this->ops[index++];
    }
    if (GlobalV::VNL_IN_H)
    {
        delete (NonlocalPW*)this->ops[index++];
    }
    delete (MetaPW*)this->ops[index++];
}

void HamiltPW::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");

    for(int iter = 0; iter < this->ops.size(); ++iter)
    {
        this->ops[iter]->init(ik);
    }

    return;
}

void HamiltPW::hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    ModuleBase::timer::tick("HamiltPW", "h_psi");
    for(int iter = 0; iter < this->ops.size(); ++iter)
    {
        this->ops[iter]->act(psi_in, hpsi, size);
    }
    ModuleBase::timer::tick("HamiltPW", "h_psi");
    return;
}

void HamiltPW::sPsi
(
    const std::complex<double> *psi,
    std::complex<double> *spsi,
    size_t size
) const
{
    for (size_t i=0; i<size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

} // namespace hamilt
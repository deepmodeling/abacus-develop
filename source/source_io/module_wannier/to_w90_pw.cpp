#include "to_w90_pw.h"

#include "source_io/module_parameter/parameter.h"

toW90_PW::toW90_PW(
    const bool &out_wannier_mmn, 
    const bool &out_wannier_amn, 
    const bool &out_wannier_unk, 
    const bool &out_wannier_eig,
    const bool &out_wannier_wvfn_formatted, 
    const std::string &nnkpfile,
    const std::string &wannier_spin,
    const int &nspin,
    const int &nbands,
    const int &nqx,
    const double &dq,
    const int &npol
):toW90(out_wannier_mmn, out_wannier_amn, out_wannier_unk, out_wannier_eig, out_wannier_wvfn_formatted, nnkpfile, wannier_spin, nspin, nbands, nqx, dq, npol)
{

}

toW90_PW::~toW90_PW()
{
    
}

void toW90_PW::calculate(
    const UnitCell& ucell,
    const ModuleBase::matrix& ekb,
    const ModulePW::PW_Basis_K* wfcpw,
    const ModulePW::PW_Basis_Big* bigpw,
    const K_Vectors& kv,
    const psi::Psi<std::complex<double>>* psi
)
{
    read_nnkp(ucell,kv);

    if (nspin_ == 2)
    {
        if (wannier_spin == "up")
        {
            start_k_index = 0;
        }
        else if (wannier_spin == "down")
        {
            start_k_index = num_kpts / 2;
        }
        else
        {
            ModuleBase::WARNING_QUIT("toW90::calculate", "Error wannier_spin set,is not \"up\" or \"down\" ");
        }
    }

    if (out_wannier_eig)
    {
        out_eig(ekb);
    }

    if (out_wannier_mmn)
    {
        cal_Mmn(*psi, wfcpw);
    }

    if (out_wannier_amn)
    {
        cal_Amn(*psi, wfcpw);
    }

    if (out_wannier_unk)
    {
        out_unk(*psi, wfcpw, bigpw);
    }

}
void toW90_PW::set_tpiba_omega(const double& tpiba, const double& omega)
{
    this->tpiba = &tpiba;
    this->omega = &omega;
    
}

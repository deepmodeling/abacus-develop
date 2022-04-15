#include "hsolverpw.h"
#include "module_base/tool_quit.h"
#include "module_elecstate/elecstatepw.h"
#include "diagocg.h"
#include "diagodavid.h"

#include "src_pw/global.h"

namespace ModuleHSolver
{


/*void HSolverPW::init(const PW_Basis* pbas_in)
{
    this->pbas = pbas_in;
    return;
}

void HSolverPW::update()
{
    return;
}*/

void HSolverPW::solve
(
    ModuleHamilt::Hamilt* pHamilt, 
    ModulePsi::Psi<std::complex<double>>& psi, 
    ModuleElecS::ElecState* pes
)
{
    //prepare for the precondition of diagonalization
    std::vector<double> precondition(psi.get_nbasis());

    // select the method of diagonalization
    if(this->method == "cg") pdiagh = new DiagoCG(&(GlobalC::hm.hpw), precondition.data());
    else if(this->method == "david") pdiagh = new DiagoDavid(&(GlobalC::hm.hpw), precondition.data());
    else ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");

    ///Loop over k points for solve Hamiltonian to charge density 
    for(int ik=0;ik<psi.get_nk();++ik)
    {
        ///update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);

        //template add precondition calculating here
        update_precondition(precondition, psi.get_current_nbas(), GlobalC::wf.g2kin);

        ///solve eigenvector and eigenvalue for H(k)
        double *p_eigenvalues = GlobalC::wf.ekb[ik];
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
        ///calculate the contribution of Psi for charge density rho
        pes->updateRhoK(psi);
    }
}

void HSolverPW::hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

void HSolverPW::update_precondition(std::vector<double> h_diag, const int npw, const double* g2kin)
{
    int precondition_type=2;
    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type==1)
    {
        for (int ig = 0;ig < npw; ig++)
        {
            h_diag[ig] = std::max(1.0, g2kin[ig]);
        }
    }
    else if (precondition_type==2)
    {
        for (int ig = 0;ig < npw; ig++)
        {
            h_diag[ig] = 1 + g2kin[ig] + sqrt( 1 + (g2kin[ig] - 1) * (g2kin[ig] - 1));
        }
    }
}

}
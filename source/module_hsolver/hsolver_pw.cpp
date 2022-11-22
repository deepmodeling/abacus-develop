#include "hsolver_pw.h"

#include "diago_cg.h"
#include "diago_david.h"
#include "diago_iter_assist.h"
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
#include "module_hamilt/hamilt_pw.h"
#include "module_elecstate/elecstate_pw.h"
#include "src_pw/global.h"
#include <algorithm>

namespace hsolver
{
HSolverPW::HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in)
{
    this->wfc_basis = wfc_basis_in;
    this->classname = "HSolverPW";
    this->diag_ethr = GlobalV::PW_DIAG_THR;
    /*this->init(pbas_in);*/
}
/*void HSolverPW::init(const PW_Basis* pbas_in)
{
    this->pbas = pbas_in;
    return;
}

void HSolverPW::update()
{
    return;
}*/
void HSolverPW::initDiagh()
{
    if (this->method == "cg")
    {
        if(pdiagh!=nullptr)
        {
            if(pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoCG<double>(precondition.data());
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoCG<double>(precondition.data());
            pdiagh->method = this->method;
            // temperary added for debugging!
            #if defined(__CUDA) || defined(__ROCM)
            gpu_diagh = new DiagoCG<double, psi::DEVICE_GPU>(precondition.data());
            gpu_diagh->method = this->method;
            #endif
        }
    }
    else if (this->method == "dav")
    {
        DiagoDavid<double>::PW_DIAG_NDIM = GlobalV::PW_DIAG_NDIM;
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoDavid<double>( precondition.data());
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoDavid<double>( precondition.data());
            pdiagh->method = this->method;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");
    }
}

void HSolverPW::solve(hamilt::Hamilt<double>* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");
    // prepare for the precondition of diagonalization
    this->precondition.resize(psi.get_nbasis());

    // select the method of diagonalization
    this->method = method_in;
    this->initDiagh();

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        this->updatePsiK(pHamilt, psi, ik);

        // template add precondition calculating here
        update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
        /// calculate the contribution of Psi for charge density rho
    }

    this->endDiagh();

    if(skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    pes->psiToRho(psi);

    ModuleBase::timer::tick("HSolverPW", "solve");
    return;
}

void HSolverPW::endDiagh()
{
    // DiagoCG would keep 9*nbasis memory in cache during loop-k
    // it should be deleted before calculating charge
    if(this->method == "cg")
    {
        delete (DiagoCG<double>*)pdiagh;
        pdiagh = nullptr;
        #if defined(__CUDA) || defined(__ROCM)
        delete (DiagoCG<double, psi::DEVICE_GPU>*)gpu_diagh;
        gpu_diagh = nullptr;
        #endif
    }
    if(this->method == "dav")
    {
        delete (DiagoDavid<double>*)pdiagh;
        pdiagh = nullptr;
    }

    //in PW base, average iteration steps for each band and k-point should be printing
    GlobalV::ofs_running<< "Average iterative diagonalization steps: "<<DiagoIterAssist<double>::avg_iter / this->wfc_basis->nks
        <<" ; where current threshold is: "<<DiagoIterAssist<double>::PW_DIAG_THR<<" . "<<std::endl;
    //reset avg_iter
    DiagoIterAssist<double>::avg_iter = 0.0;

    //psi only should be initialed once for PW
    if(!this->initialed_psi)
    {
        this->initialed_psi = true;
    }
}

void HSolverPW::updatePsiK(hamilt::Hamilt<double>* pHamilt, psi::Psi<std::complex<double>>& psi, const int ik)
{
    psi.fix_k(ik);
    if(!this->initialed_psi)
    {
        if(GlobalV::BASIS_TYPE=="pw")
        {
            // generate PAOs first, then diagonalize to get
            // inital wavefunctions.
            GlobalC::wf.diago_PAO_in_pw_k2(ik, psi, pHamilt);
        }
        else
        {
            ModuleBase::WARNING_QUIT("HSolverPW::updatePsiK", "lcao_in_pw is not supported now.");
        }
        return;
    }
}

void HSolverPW::hamiltSolvePsiK(hamilt::Hamilt<double>* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    /// a huge victory here!
    /// psi::Psi<std::complex<double>, psi::DEVICE_GPU> gpu_psi = psi;
    /// psi::Psi<std::complex<double>, psi::DEVICE_CPU> cpu_psi = gpu_psi;
    /// pdiagh->diag(hm, psi, eigenvalue);
    /// psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_CPU>()(
    ///     psi.get_device(),
    ///     cpu_psi.get_device(),
    ///     psi.get_pointer() - psi.get_psi_bias(),
    ///     cpu_psi.get_pointer() - cpu_psi.get_psi_bias(),
    ///     psi.size());


    /// hamilt::Hamilt<double, psi::DEVICE_CPU>* h_phm_in =
    ///         new hamilt::HamiltPW<double, psi::DEVICE_CPU>(
    ///                 reinterpret_cast<hamilt::HamiltPW<double, psi::DEVICE_GPU>*>(d_phm_in));
    ///
    /// pdiagh->diag(h_phm_in, psi, eigenvalue);
    /// new era
    /// if this works, then everything done!

#if defined(__CUDA) || defined(__ROCM)
    psi::Psi<std::complex<double>, psi::DEVICE_GPU> gpu_psi = psi;
    hamilt::Hamilt<double, psi::DEVICE_GPU>* d_phm_in =
            new hamilt::HamiltPW<double, psi::DEVICE_GPU>(
                    reinterpret_cast<hamilt::HamiltPW<double, psi::DEVICE_CPU>*>(hm));
    gpu_diagh->diag(d_phm_in, gpu_psi, eigenvalue);
    psi::memory::synchronize_memory_op<std::complex<double>, psi::DEVICE_CPU, psi::DEVICE_GPU>()(
         psi.get_device(),
         gpu_psi.get_device(),
         psi.get_pointer() - psi.get_psi_bias(),
         gpu_psi.get_pointer() - gpu_psi.get_psi_bias(),
         psi.size());
    delete reinterpret_cast<hamilt::HamiltPW<double, psi::DEVICE_GPU>*>(d_phm_in);
#else
    pdiagh->diag(hm, psi, eigenvalue);
#endif
}

void HSolverPW::update_precondition(std::vector<double> &h_diag, const int ik, const int npw)
{
    h_diag.assign(h_diag.size(), 1.0);
    int precondition_type = 2;
    const double tpiba2 = this->wfc_basis->tpiba2;
    
    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type == 1)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            double g2kin = this->wfc_basis->getgk2(ik,ig) * tpiba2;    
            h_diag[ig] = std::max(1.0, g2kin);
        }
    }
    else if (precondition_type == 2)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            double g2kin = this->wfc_basis->getgk2(ik,ig) * tpiba2;
            h_diag[ig] = 1 + g2kin + sqrt(1 + (g2kin - 1) * (g2kin - 1));
        }
    }
    if(GlobalV::NSPIN==4)
    {
        const int size = h_diag.size();
        for (int ig = 0; ig < npw; ig++)
        {
            h_diag[ig+size/2] = h_diag[ig];
        }
    }
}

double HSolverPW::cal_hsolerror()
{
    return this->diag_ethr * std::max(1.0, GlobalV::nelec);
}

double HSolverPW::set_diagethr(const int istep, const int iter, const double drho)
{
    //It is too complex now and should be modified.
    if (iter == 1)
    {
        if (abs(this->diag_ethr - 1.0e-2) < 1.0e-10)
        {
            if (GlobalV::init_chg == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict this->diag_ethr in the input file ()diago_the_init
                //======================================================
                this->diag_ethr = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                this->diag_ethr = 1.0e-2;
            }
        }
        // if (GlobalV::FINAL_SCF) this->diag_ethr = 1.0e-2;
        if (GlobalV::CALCULATION == "md" || GlobalV::CALCULATION == "relax" || GlobalV::CALCULATION == "cell-relax")
        {
            this->diag_ethr = std::max(this->diag_ethr, GlobalV::PW_DIAG_THR);
        }
    }
    else
    {
        if (iter == 2)
        {
            this->diag_ethr = 1.e-2;
        }
        this->diag_ethr = std::min(this->diag_ethr, 0.1 * drho / std::max(1.0, GlobalV::nelec));
    }
    return this->diag_ethr;
}

double HSolverPW::reset_diagethr(std::ofstream& ofs_running, const double hsover_error, const double drho)
{
    ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
    ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");
    ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
    ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
    this->diag_ethr = 0.1 * drho / GlobalV::nelec;
    ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
    return this->diag_ethr;
}


} // namespace hsolver
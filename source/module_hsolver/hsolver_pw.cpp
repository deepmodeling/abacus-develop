#include "hsolver_pw.h"

#include "diago_cg.h"
#include "diago_david.h"
#include "diago_iter_assist.h"
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include <algorithm>
#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#endif
namespace hsolver {

template <typename FPTYPE, typename Device>
HSolverPW<FPTYPE, Device>::HSolverPW(ModulePW::PW_Basis_K* wfc_basis_in, wavefunc* pwf_in)
{
    this->classname = "HSolverPW";
    this->wfc_basis = wfc_basis_in;
    this->pwf = pwf_in;
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
template<typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::initDiagh()
{
    if (this->method == "cg")
    {
        if(this->pdiagh!=nullptr)
        {
            if(this->pdiagh->method != this->method)
            {
                delete (DiagoCG<FPTYPE, Device>*)this->pdiagh;
                this->pdiagh = new DiagoCG<FPTYPE, Device>(precondition.data());
                this->pdiagh->method = this->method;
            }
        }
        else
        {
            this->pdiagh = new DiagoCG<FPTYPE, Device>(precondition.data());
            this->pdiagh->method = this->method;
        }
    }
    else if (this->method == "dav")
    {
        DiagoDavid<double>::PW_DIAG_NDIM = GlobalV::PW_DIAG_NDIM;
        if (this->pdiagh != nullptr)
        {
            if (this->pdiagh->method != this->method)
            {
                delete (DiagoDavid<FPTYPE, Device>*)this->pdiagh;
                this->pdiagh = new DiagoDavid<FPTYPE, Device>(precondition.data());
                this->pdiagh->method = this->method;
            }
        }
        else
        {
            this->pdiagh = new DiagoDavid<FPTYPE, Device>( precondition.data());
            this->pdiagh->method = this->method;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");
    }
}

template <typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::solve(hamilt::Hamilt<FPTYPE, Device>* pHamilt,
                                      psi::Psi<std::complex<FPTYPE>, Device>& psi,
                                      elecstate::ElecState* pes,
                                      const std::string method_in,
                                      const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");
    // prepare for the precondition of diagonalization
    this->precondition.resize(psi.get_nbasis());

    // select the method of diagonalization
    this->method = method_in;
    this->initDiagh();
    std::vector<FPTYPE> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);
    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        this->updatePsiK(pHamilt, psi, ik);

        // template add precondition calculating here
        update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);

        /// solve eigenvector and eigenvalue for H(k)
        this->hamiltSolvePsiK(pHamilt, psi, eigenvalues.data() + ik * pes->ekb.nc);
        if(skip_charge)
        {
            GlobalV::ofs_running<< "Average iterative diagonalization steps for k-points "<<ik<<" is: "<<DiagoIterAssist<FPTYPE, Device>::avg_iter
                <<" ; where current threshold is: "<<DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR<<" . "<<std::endl;
            DiagoIterAssist<FPTYPE, Device>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    castmem_2d_2h_op()(cpu_ctx, cpu_ctx, pes->ekb.c, eigenvalues.data(), pes->ekb.nr * pes->ekb.nc);

    this->endDiagh();

    if(skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<FPTYPE, Device>*>(pes)->psiToRho(psi);

#ifdef USE_PAW
    if(GlobalV::use_paw)
    {
        if(typeid(FPTYPE) != typeid(double))
        {
            ModuleBase::WARNING_QUIT("HSolverPW::solve", "PAW is only supported for double precision!");
        }

        GlobalC::paw_cell.reset_rhoij();
        for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
        {
            psi.fix_k(ik);
            int nbands = psi.get_nbands();
            for(int ib = 0; ib < nbands; ib ++)
            {
                GlobalC::paw_cell.accumulate_rhoij(reinterpret_cast<std::complex<double>*> (psi.get_pointer(ib)), pes->wg(ik,ib));
            }
        }

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;

        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for(int iat = 0; iat < GlobalC::ucell.nat; iat ++)
        {
            GlobalC::paw_cell.set_rhoij(iat,nrhoijsel[iat],rhoijp[iat].size(),rhoijselect[iat].data(),rhoijp[iat].data());
        }
    }
#endif

    ModuleBase::timer::tick("HSolverPW", "solve");
    return;
}

template<typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::endDiagh()
{
    // DiagoCG would keep 9*nbasis memory in cache during loop-k
    // it should be deleted before calculating charge
    if(this->method == "cg")
    {
        delete (DiagoCG<FPTYPE, Device>*)this->pdiagh;
        this->pdiagh = nullptr;
    }
    if(this->method == "dav")
    {
        delete (DiagoDavid<FPTYPE, Device>*)this->pdiagh;
        this->pdiagh = nullptr;
    }

    //in PW base, average iteration steps for each band and k-point should be printing
    if(DiagoIterAssist<FPTYPE, Device>::avg_iter > 0.0)
    {
        GlobalV::ofs_running<< "Average iterative diagonalization steps: "<<DiagoIterAssist<FPTYPE, Device>::avg_iter / this->wfc_basis->nks
            <<" ; where current threshold is: "<<DiagoIterAssist<FPTYPE, Device>::PW_DIAG_THR<<" . "<<std::endl;
        //reset avg_iter
        DiagoIterAssist<FPTYPE, Device>::avg_iter = 0.0;
    }
    //psi only should be initialed once for PW
    if(!this->initialed_psi)
    {
        this->initialed_psi = true;
    }
}

template <typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::updatePsiK(hamilt::Hamilt<FPTYPE, Device>* pHamilt,
                                           psi::Psi<std::complex<FPTYPE>, Device>& psi,
                                           const int ik)
{
    psi.fix_k(ik);
    if(!this->initialed_psi)
    {
        if(GlobalV::BASIS_TYPE=="pw")
        {
            // generate PAOs first, then diagonalize to get
            // inital wavefunctions.
            hamilt::diago_PAO_in_pw_k2(this->ctx, ik, psi, this->wfc_basis, this->pwf, pHamilt);
        }
        else
        {
            ModuleBase::WARNING_QUIT("HSolverPW::updatePsiK", "lcao_in_pw is not supported now.");
        }
    }
}

template<typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::hamiltSolvePsiK(hamilt::Hamilt<FPTYPE, Device>* hm, psi::Psi<std::complex<FPTYPE>, Device>& psi, FPTYPE* eigenvalue)
{
    this->pdiagh->diag(hm, psi, eigenvalue);
}

template<typename FPTYPE, typename Device>
void HSolverPW<FPTYPE, Device>::update_precondition(std::vector<FPTYPE> &h_diag, const int ik, const int npw)
{
    h_diag.assign(h_diag.size(), 1.0);
    int precondition_type = 2;
    const auto tpiba2 = static_cast<FPTYPE>(this->wfc_basis->tpiba2);

    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type == 1)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            FPTYPE g2kin = static_cast<FPTYPE>(this->wfc_basis->getgk2(ik,ig)) * tpiba2;
            h_diag[ig] = std::max(static_cast<FPTYPE>(1.0), g2kin);
        }
    }
    else if (precondition_type == 2)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            FPTYPE g2kin = static_cast<FPTYPE>(this->wfc_basis->getgk2(ik,ig)) * tpiba2;
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

template<typename FPTYPE, typename Device>
FPTYPE HSolverPW<FPTYPE, Device>::cal_hsolerror()
{
    return this->diag_ethr * static_cast<FPTYPE>(std::max(1.0, GlobalV::nelec));
}

template<typename FPTYPE, typename Device>
FPTYPE HSolverPW<FPTYPE, Device>::set_diagethr(const int istep, const int iter, const FPTYPE drho)
{
    //It is too complex now and should be modified.
    if (iter == 1)
    {
        if (std::abs(this->diag_ethr - 1.0e-2) < 1.0e-6)
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
            this->diag_ethr = std::max(this->diag_ethr, static_cast<FPTYPE>(GlobalV::PW_DIAG_THR));
        }
    }
    else
    {
        if (iter == 2)
        {
            this->diag_ethr = 1.e-2;
        }
        this->diag_ethr = std::min(this->diag_ethr, static_cast<FPTYPE>(0.1) * drho / std::max(static_cast<FPTYPE>(1.0), static_cast<FPTYPE>(GlobalV::nelec)));
    }
    // It is essential for single precision implementation to keep the diag_ethr value
    // less or equal to the single-precision limit of convergence(0.5e-4).
    // modified by denghuilu at 2023-05-15
    if (GlobalV::precision_flag == "single") {
        this->diag_ethr = std::max(this->diag_ethr, static_cast<FPTYPE>(0.5e-4));
    }
    return this->diag_ethr;
}

template<typename FPTYPE, typename Device>
FPTYPE HSolverPW<FPTYPE, Device>::reset_diagethr(std::ofstream& ofs_running, const FPTYPE hsover_error, const FPTYPE drho)
{
    ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
    ModuleBase::WARNING("scf", "Threshold on eigenvalues was too large.");
    ofs_running << " hsover_error=" << hsover_error << " > DRHO=" << drho << std::endl;
    ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
    this->diag_ethr = 0.1 * drho / GlobalV::nelec;
    ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
    return this->diag_ethr;
}

template class HSolverPW<float, psi::DEVICE_CPU>;
template class HSolverPW<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HSolverPW<float, psi::DEVICE_GPU>;
template class HSolverPW<double, psi::DEVICE_GPU>;
#endif

} // namespace hsolver
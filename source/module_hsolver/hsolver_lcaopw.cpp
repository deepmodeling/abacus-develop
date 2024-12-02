#include "hsolver_lcaopw.h"

#include "module_base/global_variable.h"
#include "module_base/parallel_global.h" // for MPI
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_psi/wavefunc.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#endif

#ifdef __EXX
#include "module_hamilt_pw/hamilt_pwdft/hamilt_lcaopw.h"
#endif

namespace hsolver
{

#ifdef USE_PAW
template <typename T>
void HSolverLIP<T>::paw_func_in_kloop(const int ik,const double tpiba)
{
    if (PARAM.inp.use_paw)
    {
        const int npw = this->wfc_basis->npwk[ik];
        ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
        for (int ig = 0; ig < npw; ig++)
        {
            _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
        }

        std::vector<double> kpt(3, 0);
        kpt[0] = this->wfc_basis->kvec_c[ik].x;
        kpt[1] = this->wfc_basis->kvec_c[ik].y;
        kpt[2] = this->wfc_basis->kvec_c[ik].z;

        double** kpg;
        double** gcar;
        kpg = new double*[npw];
        gcar = new double*[npw];
        for (int ipw = 0; ipw < npw; ipw++)
        {
            kpg[ipw] = new double[3];
            kpg[ipw][0] = _gk[ipw].x;
            kpg[ipw][1] = _gk[ipw].y;
            kpg[ipw][2] = _gk[ipw].z;

            gcar[ipw] = new double[3];
            gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
            gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
            gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
        }

        GlobalC::paw_cell.set_paw_k(npw,
                                    wfc_basis->npwk_max,
                                    kpt.data(),
                                    this->wfc_basis->get_ig2ix(ik).data(),
                                    this->wfc_basis->get_ig2iy(ik).data(),
                                    this->wfc_basis->get_ig2iz(ik).data(),
                                    (const double**)kpg,
                                    tpiba,
                                    (const double**)gcar);

        std::vector<double>().swap(kpt);
        for (int ipw = 0; ipw < npw; ipw++)
        {
            delete[] kpg[ipw];
            delete[] gcar[ipw];
        }
        delete[] kpg;
        delete[] gcar;

        GlobalC::paw_cell.get_vkb();

        GlobalC::paw_cell.set_currentk(ik);
    }
}

template <typename T>
void HSolverLIP<T>::paw_func_after_kloop(psi::Psi<T>& psi, 
                                         elecstate::ElecState* pes,
                                         const double tpiba,
                                         const int nat)
{
    if (PARAM.inp.use_paw)
    {
        if (typeid(Real) != typeid(double))
        {
            ModuleBase::WARNING_QUIT("HSolverLIP::solve", "PAW is only supported for double precision!");
        }

        GlobalC::paw_cell.reset_rhoij();
        for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
        {
            const int npw = this->wfc_basis->npwk[ik];
            ModuleBase::Vector3<double>* _gk = new ModuleBase::Vector3<double>[npw];
            for (int ig = 0; ig < npw; ig++)
            {
                _gk[ig] = this->wfc_basis->getgpluskcar(ik, ig);
            }

            std::vector<double> kpt(3, 0);
            kpt[0] = this->wfc_basis->kvec_c[ik].x;
            kpt[1] = this->wfc_basis->kvec_c[ik].y;
            kpt[2] = this->wfc_basis->kvec_c[ik].z;

            double** kpg;
            double** gcar;
            kpg = new double*[npw];
            gcar = new double*[npw];
            for (int ipw = 0; ipw < npw; ipw++)
            {
                kpg[ipw] = new double[3];
                kpg[ipw][0] = _gk[ipw].x;
                kpg[ipw][1] = _gk[ipw].y;
                kpg[ipw][2] = _gk[ipw].z;

                gcar[ipw] = new double[3];
                gcar[ipw][0] = this->wfc_basis->getgcar(ik, ipw).x;
                gcar[ipw][1] = this->wfc_basis->getgcar(ik, ipw).y;
                gcar[ipw][2] = this->wfc_basis->getgcar(ik, ipw).z;
            }

            GlobalC::paw_cell.set_paw_k(npw,
                                        wfc_basis->npwk_max,
                                        kpt.data(),
                                        this->wfc_basis->get_ig2ix(ik).data(),
                                        this->wfc_basis->get_ig2iy(ik).data(),
                                        this->wfc_basis->get_ig2iz(ik).data(),
                                        (const double**)kpg,
                                        tpiba,
                                        (const double**)gcar);

            std::vector<double>().swap(kpt);
            for (int ipw = 0; ipw < npw; ipw++)
            {
                delete[] kpg[ipw];
                delete[] gcar[ipw];
            }
            delete[] kpg;
            delete[] gcar;

            GlobalC::paw_cell.get_vkb();

            psi.fix_k(ik);
            GlobalC::paw_cell.set_currentk(ik);
            int nbands = psi.get_nbands();
            for (int ib = 0; ib < nbands; ib++)
            {
                GlobalC::paw_cell.accumulate_rhoij(reinterpret_cast<std::complex<double>*>(psi.get_pointer(ib)),
                                                   pes->wg(ik, ib));
            }
        }

        std::vector<std::vector<double>> rhoijp;
        std::vector<std::vector<int>> rhoijselect;
        std::vector<int> nrhoijsel;

#ifdef __MPI
        if (GlobalV::RANK_IN_POOL == 0)
        {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < nat; iat++)
            {
                GlobalC::paw_cell.set_rhoij(iat,
                                            nrhoijsel[iat],
                                            rhoijselect[iat].size(),
                                            rhoijselect[iat].data(),
                                            rhoijp[iat].data());
            }
        }
#else
        GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

        for (int iat = 0; iat < nat; iat++)
        {
            GlobalC::paw_cell.set_rhoij(iat,
                                        nrhoijsel[iat],
                                        rhoijselect[iat].size(),
                                        rhoijselect[iat].data(),
                                        rhoijp[iat].data());
        }

#endif
        double* nhatgr;
        GlobalC::paw_cell.get_nhat(pes->charge->nhat, nhatgr);
    }
}
#endif


/*
    lcao_in_pw
*/
template <typename T>
void HSolverLIP<T>::solve(hamilt::Hamilt<T>* pHamilt, // ESolver_KS_PW::p_hamilt
                          psi::Psi<T>& psi,           // ESolver_KS_PW::kspw_psi
                          elecstate::ElecState* pes,  // ESolver_KS_PW::pes
                          psi::Psi<T>& transform,
                          const bool skip_charge,
                          const double tpiba,
                          const int nat)
{
    ModuleBase::TITLE("HSolverLIP", "solve");
    ModuleBase::timer::tick("HSolverLIP", "solve");
    std::vector<Real> eigenvalues(pes->ekb.nr * pes->ekb.nc, 0);
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

#ifdef USE_PAW
        this->paw_func_in_kloop(ik,tpiba);
#endif

        psi.fix_k(ik);
        transform.fix_k(ik);

#ifdef __EXX
        auto& exx_lip = dynamic_cast<hamilt::HamiltLIP<T>*>(pHamilt)->exx_lip;
        auto add_exx_to_subspace_hamilt = [&ik, &exx_lip](T* hcc, const int naos) -> void {
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                for (int n = 0; n < naos; ++n)
                {
                    for (int m = 0; m < naos; ++m)
                    {
                        hcc[n * naos + m]
                            += (T)GlobalC::exx_info.info_global.hybrid_alpha * exx_lip.get_exx_matrix()[ik][m][n];
                    }
                }
            }
        };
        auto set_exxlip_lcaowfc = [&ik, &exx_lip](const T* const vcc, const int naos, const int nbands) -> void {
            if (GlobalC::exx_info.info_global.cal_exx)
            {
                exx_lip.set_hvec(ik, vcc, naos, nbands);
            }
        };
#endif
        /// solve eigenvector and eigenvalue for H(k)
        hsolver::DiagoIterAssist<T>::diagH_subspace_init(
            pHamilt,                 // interface to hamilt
            transform.get_pointer(), // transform matrix between lcao and pw
            transform.get_nbands(),
            transform.get_nbasis(),
            psi,                                  // psi in pw basis
            eigenvalues.data() + ik * pes->ekb.nc // eigenvalues
#ifdef __EXX
            ,
            add_exx_to_subspace_hamilt,
            set_exxlip_lcaowfc
#endif
        );

        if (skip_charge)
        {
            GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                 << " is: " << DiagoIterAssist<T>::avg_iter
                                 << " ; where current threshold is: " << DiagoIterAssist<T>::PW_DIAG_THR << " . "
                                 << std::endl;
            DiagoIterAssist<T>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, base_device::DEVICE_CPU>()(
        cpu_ctx,
        cpu_ctx,
        pes->ekb.c,
        eigenvalues.data(),
        pes->ekb.nr * pes->ekb.nc);

    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->calculate_weights();
    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->calEBand();
    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverLIP", "solve");
        return;
    }
    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->psiToRho(psi);

#ifdef USE_PAW
    this->paw_func_after_kloop(psi, pes,tpiba,nat);
#endif

    ModuleBase::timer::tick("HSolverLIP", "solve");
    return;
}

template class HSolverLIP<std::complex<float>>;
template class HSolverLIP<std::complex<double>>;

} // namespace hsolver

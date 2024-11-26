#include "hsolver_pw.h"

#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_elecstate/elecstate_pw.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hamilt_pw/hamilt_pwdft/wavefunc.h"
#include "module_hsolver/diag_comm_info.h"
#include "module_hsolver/diago_bpcg.h"
#include "module_hsolver/diago_cg.h"
#include "module_hsolver/diago_dav_subspace.h"
#include "module_hsolver/diago_david.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

#include <algorithm>
#include <vector>

#ifdef USE_PAW
#include "module_cell/module_paw/paw_cell.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
// #include "module_base/parallel_global.h" // for MPI
// #include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#endif
namespace hsolver
{

#ifdef USE_PAW
template <typename T, typename Device>
void HSolverPW<T, Device>::paw_func_in_kloop(const int ik)
{
    if (this->use_paw)
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
                                    GlobalC::ucell.tpiba,
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

template <typename T, typename Device>
void HSolverPW<T, Device>::call_paw_cell_set_currentk(const int ik)
{
    if (this->use_paw)
    {
        GlobalC::paw_cell.set_currentk(ik);
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::paw_func_after_kloop(psi::Psi<T, Device>& psi, elecstate::ElecState* pes)
{
    if (this->use_paw)
    {
        if (typeid(Real) != typeid(double))
        {
            ModuleBase::WARNING_QUIT("HSolverPW::solve", "PAW is only supported for double precision!");
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
                                        GlobalC::ucell.tpiba,
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
        if (this->rank_in_pool == 0)
        {
            GlobalC::paw_cell.get_rhoijp(rhoijp, rhoijselect, nrhoijsel);

            for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
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

        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
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

template <typename T, typename Device>
void HSolverPW<T, Device>::cal_ethr_band(const double& wk,
                                         const double* wg,
                                         const double& ethr,
                                         std::vector<double>& ethrs)
{
    // threshold for classifying occupied and unoccupied bands
    const double occ_threshold = 1e-2;
    // diagonalization threshold limitation for unoccupied bands
    const double ethr_limit = 1e-5;
    if (wk > 0.0)
    {
        // Note: the idea of threshold for unoccupied bands (1e-5) comes from QE
        // In ABACUS, We applied a smoothing process to this truncation to avoid abrupt changes in energy errors between
        // different bands.
        const double ethr_unocc = std::max(ethr_limit, ethr);
        for (int i = 0; i < ethrs.size(); i++)
        {
            double band_weight = wg[i] / wk;
            if (band_weight > occ_threshold)
            {
                ethrs[i] = ethr;
            }
            else if (band_weight > ethr_limit)
            { // similar energy difference for different bands when band_weight in range [1e-5, 1e-2]
                ethrs[i] = std::min(ethr_unocc, ethr / band_weight);
            }
            else
            {
                ethrs[i] = ethr_unocc;
            }
        }
    }
    else
    {
        for (int i = 0; i < ethrs.size(); i++)
        {
            ethrs[i] = ethr;
        }
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::solve(hamilt::Hamilt<T, Device>* pHamilt,
                                 psi::Psi<T, Device>& psi,
                                 elecstate::ElecState* pes,
                                 double* out_eigenvalues,
                                 const int rank_in_pool_in,
                                 const int nproc_in_pool_in,
                                 const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");

    this->rank_in_pool = rank_in_pool_in;
    this->nproc_in_pool = nproc_in_pool_in;

    // report if the specified diagonalization method is not supported
    const std::initializer_list<std::string> _methods = {"cg", "dav", "dav_subspace", "bpcg"};
    if (std::find(std::begin(_methods), std::end(_methods), this->method) == std::end(_methods))
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This type of eigensolver is not supported!");
    }

    // prepare for the precondition of diagonalization
    std::vector<Real> precondition(psi.get_nbasis(), 0.0);
    std::vector<Real> eigenvalues(this->wfc_basis->nks * psi.get_nbands(), 0.0);
    ethr_band.resize(psi.get_nbands(), DiagoIterAssist<T, Device>::PW_DIAG_THR);

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

#ifdef USE_PAW
        this->paw_func_in_kloop(ik);
#endif

        /// update psi pointer for each k point
        psi.fix_k(ik);

        // template add precondition calculating here
        update_precondition(precondition, ik, this->wfc_basis->npwk[ik], Real(pes->pot->get_vl_of_0()));

        // only dav_subspace method used smooth threshold for all bands now,
        // for other methods, this trick can be added in the future to accelerate calculation without accuracy loss.
        if (this->method == "dav_subspace")
        {
            this->cal_ethr_band(pes->klist->wk[ik],
                                &pes->wg(ik, 0),
                                DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                ethr_band);
        }

#ifdef USE_PAW
        this->call_paw_cell_set_currentk(ik);
#endif

        /// solve eigenvector and eigenvalue for H(k)
        this->hamiltSolvePsiK(pHamilt, psi, precondition, eigenvalues.data() + ik * psi.get_nbands());

        if (skip_charge)
        {
            GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                 << " is: " << DiagoIterAssist<T, Device>::avg_iter
                                 << " ; where current threshold is: " << this->diag_thr << " . " << std::endl;
            DiagoIterAssist<T, Device>::avg_iter = 0.0;
        }
        /// calculate the contribution of Psi for charge density rho
    }
    // END Loop over k points

    // copy eigenvalues to ekb in ElecState
    base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, base_device::DEVICE_CPU>()(
        cpu_ctx,
        cpu_ctx,
        // pes->ekb.c,
        out_eigenvalues,
        eigenvalues.data(),
        // pes->ekb.nr * pes->ekb.nc
        this->wfc_basis->nks * psi.get_nbands());

    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->calculate_weights();
    reinterpret_cast<elecstate::ElecStatePW<T>*>(pes)->calEBand();
    if (skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    else
    {
        reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pes)->psiToRho(psi);

#ifdef USE_PAW
        this->paw_func_after_kloop(psi, pes);
#endif

        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::hamiltSolvePsiK(hamilt::Hamilt<T, Device>* hm,
                                           psi::Psi<T, Device>& psi,
                                           std::vector<Real>& pre_condition,
                                           Real* eigenvalue)
{
#ifdef __MPI
    const diag_comm_info comm_info = {POOL_WORLD, this->rank_in_pool, this->nproc_in_pool};
#else
    const diag_comm_info comm_info = {this->rank_in_pool, this->nproc_in_pool};
#endif

    if (this->method == "cg")
    {
        // wrap the subspace_func into a lambda function
        auto ngk_pointer = psi.get_ngk_pointer();
        auto subspace_func = [hm, ngk_pointer](const ct::Tensor& psi_in, ct::Tensor& psi_out) {
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim == 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_in_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                      1,
                                                      psi_in.shape().dim_size(0),
                                                      psi_in.shape().dim_size(1),
                                                      ngk_pointer);
            auto psi_out_wrapper = psi::Psi<T, Device>(psi_out.data<T>(),
                                                       1,
                                                       psi_out.shape().dim_size(0),
                                                       psi_out.shape().dim_size(1),
                                                       ngk_pointer);
            auto eigen = ct::Tensor(ct::DataTypeToEnum<Real>::value,
                                    ct::DeviceType::CpuDevice,
                                    ct::TensorShape({psi_in.shape().dim_size(0)}));

            DiagoIterAssist<T, Device>::diagH_subspace(hm, psi_in_wrapper, psi_out_wrapper, eigen.data<Real>());
        };
        DiagoCG<T, Device> cg(this->basis_type,
                              this->calculation_type,
                              this->need_subspace,
                              subspace_func,
                              this->diag_thr,
                              this->diag_iter_max,
                              this->nproc_in_pool);

        // wrap the hpsi_func and spsi_func into a lambda function
        using ct_Device = typename ct::PsiToContainer<Device>::type;

        // wrap the hpsi_func and spsi_func into a lambda function
        auto hpsi_func = [hm, ngk_pointer](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                   1,
                                                   ndim == 1 ? 1 : psi_in.shape().dim_size(0),
                                                   ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                                                   ngk_pointer);
            psi::Range all_bands_range(true, psi_wrapper.get_current_k(), 0, psi_wrapper.get_nbands() - 1);
            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_wrapper, all_bands_range, hpsi_out.data<T>());
            hm->ops->hPsi(info);
            ModuleBase::timer::tick("DiagoCG_New", "hpsi_func");
        };
        auto spsi_func = [this, hm](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");

            if (this->use_uspp)
            {
                // Convert a Tensor object to a psi::Psi object
                hm->sPsi(psi_in.data<T>(),
                         spsi_out.data<T>(),
                         ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                         ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                         ndim == 1 ? 1 : psi_in.shape().dim_size(0));
            }
            else
            {
                base_device::memory::synchronize_memory_op<T, Device, Device>()(
                    this->ctx,
                    this->ctx,
                    spsi_out.data<T>(),
                    psi_in.data<T>(),
                    static_cast<size_t>((ndim == 1 ? 1 : psi_in.shape().dim_size(0))
                                        * (ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1))));
            }

            ModuleBase::timer::tick("DiagoCG_New", "spsi_func");
        };
        auto psi_tensor = ct::TensorMap(psi.get_pointer(),
                                        ct::DataTypeToEnum<T>::value,
                                        ct::DeviceTypeToEnum<ct_Device>::value,
                                        ct::TensorShape({psi.get_nbands(), psi.get_nbasis()}));
        auto eigen_tensor = ct::TensorMap(eigenvalue,
                                          ct::DataTypeToEnum<Real>::value,
                                          ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                          ct::TensorShape({psi.get_nbands()}));
        auto prec_tensor = ct::TensorMap(pre_condition.data(),
                                         ct::DataTypeToEnum<Real>::value,
                                         ct::DeviceTypeToEnum<ct::DEVICE_CPU>::value,
                                         ct::TensorShape({static_cast<int>(pre_condition.size())}))
                               .to_device<ct_Device>()
                               .slice({0}, {psi.get_current_nbas()});

        cg.diag(hpsi_func, spsi_func, psi_tensor, eigen_tensor, prec_tensor);
        // TODO: Double check tensormap's potential problem
        ct::TensorMap(psi.get_pointer(), psi_tensor, {psi.get_nbands(), psi.get_nbasis()}).sync(psi_tensor);
    }
    else if (this->method == "bpcg")
    {
        DiagoBPCG<T, Device> bpcg(pre_condition.data());
        bpcg.init_iter(psi);
        bpcg.diag(hm, psi, eigenvalue);
    }
    else if (this->method == "dav_subspace")
    {
        auto ngk_pointer = psi.get_ngk_pointer();
        // hpsi_func (X, HX, ld, nvec) -> HX = H(X), X and HX blockvectors of size ld x nvec
        auto hpsi_func = [hm, ngk_pointer](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
            ModuleBase::timer::tick("DavSubspace", "hpsi_func");

            // Convert "pointer data stucture" to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, ngk_pointer);

            psi::Range bands_range(true, 0, 0, nvec - 1);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);

            ModuleBase::timer::tick("DavSubspace", "hpsi_func");
        };
        bool scf = this->calculation_type == "nscf" ? false : true;

        Diago_DavSubspace<T, Device> dav_subspace(pre_condition,
                                                  psi.get_nbands(),
                                                  psi.get_k_first() ? psi.get_current_nbas()
                                                                    : psi.get_nk() * psi.get_nbasis(),
                                                  PARAM.inp.pw_diag_ndim,
                                                  this->diag_thr,
                                                  this->diag_iter_max,
                                                  this->need_subspace,
                                                  comm_info);

        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(
            dav_subspace.diag(hpsi_func, psi.get_pointer(), psi.get_nbasis(), eigenvalue, this->ethr_band.data(), scf));
    }
    else if (this->method == "dav")
    {
        // Davidson iter parameters

        /// Allow 5 tries at most. If ntry > ntry_max = 5, exit diag loop.
        const int ntry_max = 5;
        /// In non-self consistent calculation, do until totally converged. Else
        /// allow 5 eigenvecs to be NOT converged.
        const int notconv_max = ("nscf" == this->calculation_type) ? 0 : 5;
        /// convergence threshold
        const Real david_diag_thr = this->diag_thr;
        /// maximum iterations
        const int david_maxiter = this->diag_iter_max;

        // dimensions of matrix to be solved
        const int dim = psi.get_current_nbas(); /// dimension of matrix
        const int nband = psi.get_nbands();     /// number of eigenpairs sought
        const int ld_psi = psi.get_nbasis();    /// leading dimension of psi

        // Davidson matrix-blockvector functions

        auto ngk_pointer = psi.get_ngk_pointer();
        /// wrap hpsi into lambda function, Matrix \times blockvector
        // hpsi_func (X, HX, ld, nvec) -> HX = H(X), X and HX blockvectors of size ld x nvec
        auto hpsi_func = [hm, ngk_pointer](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {
            ModuleBase::timer::tick("David", "hpsi_func");

            // Convert pointer of psi_in to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, ngk_pointer);

            psi::Range bands_range(true, 0, 0, nvec - 1);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);

            ModuleBase::timer::tick("David", "hpsi_func");
        };

        /// wrap spsi into lambda function, Matrix \times blockvector
        /// spsi(X, SX, ld, nvec)
        /// ld is leading dimension of psi and spsi
        auto spsi_func = [hm](const T* psi_in,
                              T* spsi_out,
                              const int ld_psi, // Leading dimension of psi and spsi.
                              const int nvec    // Number of vectors(bands)
                         ) {
            ModuleBase::timer::tick("David", "spsi_func");
            // sPsi determines S=I or not by  PARAM.globalv.use_uspp inside
            // sPsi(psi, spsi, nrow, npw, nbands)
            hm->sPsi(psi_in, spsi_out, ld_psi, ld_psi, nvec);
            ModuleBase::timer::tick("David", "spsi_func");
        };

        DiagoDavid<T, Device> david(pre_condition.data(), nband, dim, PARAM.inp.pw_diag_ndim, this->use_paw, comm_info);
        // do diag and add davidson iteration counts up to avg_iter
        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(david.diag(hpsi_func,
                                                                               spsi_func,
                                                                               ld_psi,
                                                                               psi.get_pointer(),
                                                                               eigenvalue,
                                                                               david_diag_thr,
                                                                               david_maxiter,
                                                                               ntry_max,
                                                                               notconv_max));
    }
    return;
}

template <typename T, typename Device>
void HSolverPW<T, Device>::update_precondition(std::vector<Real>& h_diag,
                                               const int ik,
                                               const int npw,
                                               const Real vl_of_0)
{
    h_diag.assign(h_diag.size(), 1.0);
    int precondition_type = 2;
    const auto tpiba2 = static_cast<Real>(this->wfc_basis->tpiba2);

    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type == 1)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            Real g2kin = static_cast<Real>(this->wfc_basis->getgk2(ik, ig)) * tpiba2;
            h_diag[ig] = std::max(static_cast<Real>(1.0), g2kin);
        }
    }
    else if (precondition_type == 2)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            Real g2kin = static_cast<Real>(this->wfc_basis->getgk2(ik, ig)) * tpiba2;

            if (this->method == "dav_subspace")
            {
                h_diag[ig] = g2kin + vl_of_0;
            }
            else
            {
                h_diag[ig] = 1 + g2kin + sqrt(1 + (g2kin - 1) * (g2kin - 1));
            }
        }
    }
    if (this->nspin == 4)
    {
        const int size = h_diag.size();
        for (int ig = 0; ig < npw; ig++)
        {
            h_diag[ig + size / 2] = h_diag[ig];
        }
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::output_iterInfo()
{
    // in PW base, average iteration steps for each band and k-point should be printing
    if (DiagoIterAssist<T, Device>::avg_iter > 0.0)
    {
        GlobalV::ofs_running << "Average iterative diagonalization steps: "
                             << DiagoIterAssist<T, Device>::avg_iter / this->wfc_basis->nks
                             << " ; where current threshold is: " << this->diag_thr << " . " << std::endl;
        // reset avg_iter
        DiagoIterAssist<T, Device>::avg_iter = 0.0;
    }
}

template class HSolverPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HSolverPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver
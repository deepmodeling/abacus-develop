#include "hsolver_pw.h"

#include "source_base/global_variable.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_estate/elecstate_pw.h"
#include "source_hamilt/hamilt.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/diago_bpcg.h"
#include "source_hsolver/diago_cg.h"
#include "source_hsolver/diago_dav_subspace.h"
#include "source_hsolver/diago_david.h"
#include "source_hsolver/diago_iter_assist.h"
#include "source_io/module_parameter/parameter.h"
#include "source_psi/psi.h"
#include "source_estate/elecstate_tools.h"


#include <algorithm>
#include <vector>

namespace hsolver
{

template <typename T, typename Device>
void HSolverPW<T, Device>::cal_smooth_ethr(const double& wk,
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
                                 const bool skip_charge,
                                 const double tpiba,
                                 const int nat)
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
    ethr_band.resize(psi.get_nbands(), this->diag_thr);

    // Initialize k-point continuity if enabled
    static int count = 0;
    if (use_k_continuity) {
        build_k_neighbors();
    }

    // Loop over k points for solve Hamiltonian to charge density
    if (use_k_continuity) {
        // K-point continuity case
        for (int i = 0; i < this->wfc_basis->nks; ++i)
        {
            const int ik = k_order[i];
            
            // update H(k) for each k point
            pHamilt->updateHk(ik);

#ifdef USE_PAW
            this->paw_func_in_kloop(ik, tpiba);
#endif

            // update psi pointer for each k point
            psi.fix_k(ik);

            // If using k-point continuity and not first k-point, propagate from parent
            if (ik > 0 && count == 0 && k_parent.find(ik) != k_parent.end()) {
                propagate_psi(psi, k_parent[ik], ik);
            }

            // template add precondition calculating here
            update_precondition(precondition, ik, this->wfc_basis->npwk[ik], Real(pes->pot->get_vl_of_0()));

            // use smooth threshold for all iter methods
            if (PARAM.inp.diago_smooth_ethr == true)
            {
                this->cal_smooth_ethr(pes->klist->wk[ik],
                                    &pes->wg(ik, 0),
                                    DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                    ethr_band);
            }

#ifdef USE_PAW
            this->call_paw_cell_set_currentk(ik);
#endif

            // solve eigenvector and eigenvalue for H(k)
            this->hamiltSolvePsiK(pHamilt, psi, precondition, eigenvalues.data() + ik * psi.get_nbands(), this->wfc_basis->nks);

            if (skip_charge)
            {
                GlobalV::ofs_running << "Average iterative diagonalization steps for k-points " << ik
                                    << " is: " << DiagoIterAssist<T, Device>::avg_iter
                                    << " ; where current threshold is: " << this->diag_thr << " . " << std::endl;
                DiagoIterAssist<T, Device>::avg_iter = 0.0;
            }
        }
    }
    else {
        // Original code without k-point continuity
        for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
        {
            // update H(k) for each k point
            pHamilt->updateHk(ik);

#ifdef USE_PAW
            this->paw_func_in_kloop(ik, tpiba);
#endif

            // update psi pointer for each k point
            psi.fix_k(ik);

            // template add precondition calculating here
            update_precondition(precondition, ik, this->wfc_basis->npwk[ik], Real(pes->pot->get_vl_of_0()));

            // use smooth threshold for all iter methods
            if (PARAM.inp.diago_smooth_ethr == true)
            {
                this->cal_smooth_ethr(pes->klist->wk[ik],
                                    &pes->wg(ik, 0),
                                    DiagoIterAssist<T, Device>::PW_DIAG_THR,
                                    ethr_band);
            }

#ifdef USE_PAW
            this->call_paw_cell_set_currentk(ik);
#endif

            // solve eigenvector and eigenvalue for H(k)
            this->hamiltSolvePsiK(pHamilt, psi, precondition, eigenvalues.data() + ik * psi.get_nbands(), this->wfc_basis->nks);

            if (skip_charge)
            {
                GlobalV::ofs_running << " k(" << ik+1 << "/" << pes->klist->get_nkstot()
                                     << ") Iter steps (avg)=" << DiagoIterAssist<T, Device>::avg_iter
                                     << " threshold=" << this->diag_thr << std::endl;
                DiagoIterAssist<T, Device>::avg_iter = 0.0;
            }
            /// calculate the contribution of Psi for charge density rho
        }
    }
    
    count++;
    // END Loop over k points

    // copy eigenvalues to ekb in ElecState
    base_device::memory::cast_memory_op<double, Real, base_device::DEVICE_CPU, base_device::DEVICE_CPU>()(
        out_eigenvalues,
        eigenvalues.data(),
        this->wfc_basis->nks * psi.get_nbands());

    auto _pes_pw = reinterpret_cast<elecstate::ElecStatePW<T>*>(pes);
    elecstate::calculate_weights(_pes_pw->ekb,
                                 _pes_pw->wg,
                                 _pes_pw->klist,
                                 _pes_pw->eferm,
                                 _pes_pw->f_en,
                                 _pes_pw->nelec_spin,
                                 _pes_pw->skip_weights);

    elecstate::calEBand(_pes_pw->ekb,_pes_pw->wg,_pes_pw->f_en);
    if (skip_charge)
    {
        if (PARAM.globalv.use_uspp)
        {
            reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pes)->cal_becsum(psi);
        }
    }
    else
    {
        reinterpret_cast<elecstate::ElecStatePW<T, Device>*>(pes)->psiToRho(psi);
    }

	ModuleBase::timer::tick("HSolverPW", "solve");
	return;
}

template <typename T, typename Device>
void HSolverPW<T, Device>::hamiltSolvePsiK(hamilt::Hamilt<T, Device>* hm,
                                           psi::Psi<T, Device>& psi,
                                           std::vector<Real>& pre_condition,
                                           Real* eigenvalue,
                                           const int& nk_nums)
{
    ModuleBase::timer::tick("HSolverPW", "solve_psik");
#ifdef __MPI
    const diag_comm_info comm_info = {POOL_WORLD, this->rank_in_pool, this->nproc_in_pool};
#else
    const diag_comm_info comm_info = {this->rank_in_pool, this->nproc_in_pool};
#endif

    const int cur_nbasis = psi.get_current_nbas();

    if (this->method == "cg")
    {
        // wrap the subspace_func into a lambda function
        auto subspace_func = [hm, cur_nbasis](const ct::Tensor& psi_in, ct::Tensor& psi_out) {
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim == 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_in_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                      1,
                                                      psi_in.shape().dim_size(0),
                                                      psi_in.shape().dim_size(1),
                                                      cur_nbasis);
            auto psi_out_wrapper = psi::Psi<T, Device>(psi_out.data<T>(),
                                                       1,
                                                       psi_out.shape().dim_size(0),
                                                       psi_out.shape().dim_size(1),
                                                       cur_nbasis);
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
        auto hpsi_func = [hm, cur_nbasis](const ct::Tensor& psi_in, ct::Tensor& hpsi_out) {
            // psi_in should be a 2D tensor:
            // psi_in.shape() = [nbands, nbasis]
            const auto ndim = psi_in.shape().ndim();
            REQUIRES_OK(ndim <= 2, "dims of psi_in should be less than or equal to 2");
            // Convert a Tensor object to a psi::Psi object
            auto psi_wrapper = psi::Psi<T, Device>(psi_in.data<T>(),
                                                   1,
                                                   ndim == 1 ? 1 : psi_in.shape().dim_size(0),
                                                   ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1),
                                                   cur_nbasis);
            psi::Range all_bands_range(true, psi_wrapper.get_current_k(), 0, psi_wrapper.get_nbands() - 1);
            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_wrapper, all_bands_range, hpsi_out.data<T>());
            hm->ops->hPsi(info);
        };
        auto spsi_func = [this, hm](const ct::Tensor& psi_in, ct::Tensor& spsi_out) {
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
                    spsi_out.data<T>(),
                    psi_in.data<T>(),
                    static_cast<size_t>((ndim == 1 ? 1 : psi_in.shape().dim_size(0))
                                        * (ndim == 1 ? psi_in.NumElements() : psi_in.shape().dim_size(1))));
            }
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
                               .slice({0}, {psi.get_current_ngk()});

        cg.diag(hpsi_func, spsi_func, psi_tensor, eigen_tensor, this->ethr_band, prec_tensor);
        // TODO: Double check tensormap's potential problem
        // ct::TensorMap(psi.get_pointer(), psi_tensor, {psi.get_nbands(), psi.get_nbasis()}).sync(psi_tensor);
    }
    else if (this->method == "bpcg")
    {
        const int nband_l = psi.get_nbands();
        const int nbasis = psi.get_nbasis();
        const int ndim = psi.get_current_ngk();
        // hpsi_func (X, HX, ld, nvec) -> HX = H(X), X and HX blockvectors of size ld x nvec
        auto hpsi_func = [hm, cur_nbasis](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {

            // Convert "pointer data stucture" to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, cur_nbasis);

            psi::Range bands_range(true, 0, 0, nvec - 1);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);
        };
        DiagoBPCG<T, Device> bpcg(pre_condition.data());
        bpcg.init_iter(PARAM.inp.nbands, nband_l, nbasis, ndim);
        bpcg.diag(hpsi_func, psi.get_pointer(), eigenvalue, this->ethr_band);
    }
    else if (this->method == "dav_subspace")
    {
        // hpsi_func (X, HX, ld, nvec) -> HX = H(X), X and HX blockvectors of size ld x nvec
        auto hpsi_func = [hm, cur_nbasis](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {

            // Convert "pointer data stucture" to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, cur_nbasis);

            psi::Range bands_range(true, 0, 0, nvec - 1);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);
        };
        bool scf = this->calculation_type == "nscf" ? false : true;

        Diago_DavSubspace<T, Device> dav_subspace(pre_condition,
                                                  psi.get_nbands(),
                                                  psi.get_k_first() ? psi.get_current_ngk()
                                                                    : psi.get_nk() * psi.get_nbasis(),
                                                  PARAM.inp.pw_diag_ndim,
                                                  this->diag_thr,
                                                  this->diag_iter_max,
                                                  this->need_subspace,
                                                  comm_info,
                                                  PARAM.inp.diag_subspace,
                                                  PARAM.inp.nb2d);

        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(
            dav_subspace.diag(hpsi_func, psi.get_pointer(), psi.get_nbasis(), eigenvalue, this->ethr_band, scf));
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
        const int dim = psi.get_current_ngk(); /// dimension of matrix
        const int nband = psi.get_nbands();            /// number of eigenpairs sought
        const int ld_psi = psi.get_nbasis();           /// leading dimension of psi

        // Davidson matrix-blockvector functions
        /// wrap hpsi into lambda function, Matrix \times blockvector
        // hpsi_func (X, HX, ld, nvec) -> HX = H(X), X and HX blockvectors of size ld x nvec
        auto hpsi_func = [hm, cur_nbasis](T* psi_in, T* hpsi_out, const int ld_psi, const int nvec) {

            // Convert pointer of psi_in to a psi::Psi object
            auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, cur_nbasis);

            psi::Range bands_range(true, 0, 0, nvec - 1);

            using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
            hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
            hm->ops->hPsi(info);
        };

        /// wrap spsi into lambda function, Matrix \times blockvector
        /// spsi(X, SX, ld, nvec)
        /// ld is leading dimension of psi and spsi
        auto spsi_func = [hm](const T* psi_in,
                              T* spsi_out,
                              const int ld_psi, // Leading dimension of psi and spsi.
                              const int nvec    // Number of vectors(bands)
                         ) {
            // sPsi determines S=I or not by  PARAM.globalv.use_uspp inside
            // sPsi(psi, spsi, nrow, npw, nbands)
            hm->sPsi(psi_in, spsi_out, ld_psi, ld_psi, nvec);
        };

        DiagoDavid<T, Device> david(pre_condition.data(), nband, dim, PARAM.inp.pw_diag_ndim, this->use_paw, comm_info);
        // do diag and add davidson iteration counts up to avg_iter
        DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(david.diag(hpsi_func,
                                                                               spsi_func,
                                                                               ld_psi,
                                                                               psi.get_pointer(),
                                                                               eigenvalue,
                                                                               this->ethr_band,
                                                                               david_maxiter,
                                                                               ntry_max,
                                                                               notconv_max));
    }
    ModuleBase::timer::tick("HSolverPW", "solve_psik");
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

template <typename T, typename Device>
void HSolverPW<T, Device>::build_k_neighbors() {
    const int nk = this->wfc_basis->nks;
    kvecs_c.resize(nk);
    k_order.clear();
    k_order.reserve(nk);
    
    // Store k-points and corresponding indices
    struct KPoint {
        ModuleBase::Vector3<double> kvec;
        int index;
        double norm;
        
        KPoint(const ModuleBase::Vector3<double>& v, int i) : 
            kvec(v), index(i), norm(v.norm()) {}
    };
    
    // Build k-point list
    std::vector<KPoint> klist;
    for (int ik = 0; ik < nk; ++ik) {
        kvecs_c[ik] = this->wfc_basis->kvec_c[ik];
        klist.push_back(KPoint(kvecs_c[ik], ik));
    }
    
    // Sort k-points by distance from origin
    std::sort(klist.begin(), klist.end(),
        [](const KPoint& a, const KPoint& b) {
            return a.norm < b.norm;
        });
    
    // Build parent-child relationships
    k_order.push_back(klist[0].index);
    
    // Find nearest processed k-point as parent for each k-point
    for (int i = 1; i < nk; ++i) {
        int current_k = klist[i].index;
        double min_dist = 1e10;
        int parent = -1;
        
        // find the nearest k-point as parent
        for (int j = 0; j < k_order.size(); ++j) {
            int processed_k = k_order[j];
            double dist = (kvecs_c[current_k] - kvecs_c[processed_k]).norm2();
            if (dist < min_dist) {
                min_dist = dist;
                parent = processed_k;
            }
        }
        
        k_parent[current_k] = parent;
        k_order.push_back(current_k);
    }
}

template <typename T, typename Device>
void HSolverPW<T, Device>::propagate_psi(psi::Psi<T, Device>& psi, const int from_ik, const int to_ik) {
    const int nbands = psi.get_nbands();
    const int npwk = this->wfc_basis->npwk[to_ik];
    
    // Get k-point difference
    ModuleBase::Vector3<double> dk = kvecs_c[to_ik] - kvecs_c[from_ik];
    
    // Allocate porter locally
    T* porter = nullptr;
    resmem_complex_op()(porter, this->wfc_basis->nmaxgr, "HSolverPW::porter");
    
    // Process each band
    for (int ib = 0; ib < nbands; ib++)
    {
        // Fix current k-point and band
        // psi.fix_k(from_ik);
        
        // FFT to real space
        // this->wfc_basis->recip_to_real(this->ctx, psi.get_pointer(ib), porter, from_ik);
        this->wfc_basis->recip_to_real(this->ctx, &psi(from_ik, ib, 0), porter, from_ik);
        
        // Apply phase factor
        //     // TODO: Check how to get the r vector
        //     ModuleBase::Vector3<double> r = this->wfc_basis->get_ir2r(ir);
        //     double phase = this->wfc_basis->tpiba * (dk.x * r.x + dk.y * r.y + dk.z * r.z);
        //     psi_real[ir] *= std::exp(std::complex<double>(0.0, phase));
        // }
        
        // Fix k-point for target
        // psi.fix_k(to_ik);
        
        // FFT back to reciprocal space
        // this->wfc_basis->real_to_recip(this->ctx, porter, psi.get_pointer(ib), to_ik, true);
        this->wfc_basis->real_to_recip(this->ctx, porter, &psi(to_ik, ib, 0), to_ik);
    }

    // Clean up porter
    delmem_complex_op()(porter);
}

template class HSolverPW<std::complex<float>, base_device::DEVICE_CPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class HSolverPW<std::complex<float>, base_device::DEVICE_GPU>;
template class HSolverPW<std::complex<double>, base_device::DEVICE_GPU>;
#endif

} // namespace hsolver

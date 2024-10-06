#include "hsolver_lrtd.h"
#include "module_parameter/parameter.h"
#include "module_hsolver/diago_david.h"
#include "module_hsolver/diago_dav_subspace.h"
#include "module_hsolver/diago_cg.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_print.h"

namespace LR
{
    inline double square(double x) { return x * x; };
    inline double square(std::complex<double> x) { return x.real() * x.real() + x.imag() * x.imag(); };
    template<typename T>
    inline void print_eigs(const std::vector<T>& eigs, const std::string& label = "", const double factor = 1.0)
    {
        std::cout << label << std::endl;
        for (auto& e : eigs) {std::cout << e * factor << " ";
}
        std::cout << std::endl;
    }
    template<typename T, typename Device>
    void HSolverLR<T, Device>::solve(hamilt::Hamilt<T, Device>* pHamilt,
        psi::Psi<T, Device>& psi,
        elecstate::ElecState* pes,
        const std::string method_in,
        const bool hermitian)
    {
        ModuleBase::TITLE("HSolverLR", "solve");
        assert(psi.get_nk() == nk);
        const std::vector<std::string> spin_types = { "Spin Singlet", "Spin Triplet" };
        // note: if not TDA, the eigenvalues will be complex
        // then we will need a new constructor of DiagoDavid

        // 1. allocate precondition and eigenvalue
        std::vector<Real> precondition(psi.get_nk() * psi.get_nbasis());
        std::vector<Real> eigenvalue(psi.get_nbands());   //nstates
        // 2. select the method
        this->method = method_in;
#ifdef __MPI
        const hsolver::diag_comm_info comm_info = { POOL_WORLD, GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL };
#else
        const hsolver::diag_comm_info comm_info = { GlobalV::RANK_IN_POOL, GlobalV::NPROC_IN_POOL };
#endif

        if (this->method == "lapack")
        {
            std::vector<T> Amat_full = pHamilt->matrix();
            eigenvalue.resize(nk * npairs);
            if (hermitian) { LR_Util::diag_lapack(nk * npairs, Amat_full.data(), eigenvalue.data()); }
            else
            {
                std::vector<std::complex<double>> eig_complex(nk * npairs);
                LR_Util::diag_lapack_nh(nk * npairs, Amat_full.data(), eig_complex.data());
                print_eigs(eig_complex, "Right eigenvalues: of the non-Hermitian matrix: (Ry)");
                for (int i = 0; i < nk * npairs; i++) { eigenvalue[i] = eig_complex[i].real(); }
            }
            psi.fix_kb(0, 0);
            // copy eigenvectors
            for (int i = 0;i < psi.size();++i) { psi.get_pointer()[i] = Amat_full[i];
}
        }
        else
        {
            // 3. set precondition and diagethr
            for (int i = 0; i < psi.get_nk() * psi.get_nbasis(); ++i) {
                precondition[i] = static_cast<Real>(1.0);
}

            // wrap band-first psi as k1-first psi_k1_dav
            psi::Psi<T> psi_k1_dav = LR_Util::bfirst_to_k1_wrapper(psi);
            assert(psi_k1_dav.get_nbands() == psi.get_nbands());
            assert(psi_k1_dav.get_nbasis() == psi.get_nbasis() * psi.get_nk());

            const int david_maxiter = hsolver::DiagoIterAssist<T, Device>::PW_DIAG_NMAX;

            if (this->method == "dav")
            {
                // Allow 5 tries at most. If ntry > ntry_max = 5, exit diag loop.
                const int ntry_max = 5;
                // In non-self consistent calculation, do until totally converged. Else allow 5 eigenvecs to be NOT
                // converged.
                const int notconv_max = ("nscf" == PARAM.inp.calculation) ? 0 : 5;
                // do diag and add davidson iteration counts up to avg_iter

                auto hpsi_func = [pHamilt](
                    T* psi_in,
                    T* hpsi_out,
                    const int ld_psi,
                    const int nvec)
                    {
                        auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, nullptr);
                        psi::Range bands_range(true, 0, 0, nvec-1);
                        using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
                        hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
                        pHamilt->ops->hPsi(info);
                    };
                auto spsi_func = [pHamilt](const T* psi_in, T* spsi_out,
                               const int nrow, const int npw,  const int nbands){
                    // sPsi determines S=I or not by PARAM.globalv.use_uspp inside
                    pHamilt->sPsi(psi_in, spsi_out, nrow, npw, nbands);
                };

                const int& dim = psi_k1_dav.get_nbasis();   //equals to leading dimension here
                const int& nband = psi_k1_dav.get_nbands();
                hsolver::DiagoDavid<T, Device> david(precondition.data(), nband, dim, PARAM.inp.pw_diag_ndim, PARAM.inp.use_paw, comm_info);
                hsolver::DiagoIterAssist<T, Device>::avg_iter += static_cast<double>(david.diag(hpsi_func, spsi_func,
                    dim, psi_k1_dav.get_pointer(), eigenvalue.data(), this->diag_ethr, david_maxiter, ntry_max, 0));
            }
            else if (this->method == "dav_subspace") //need refactor
            {
                hsolver::Diago_DavSubspace<T, Device> dav_subspace(precondition,
                    psi_k1_dav.get_nbands(),
                    psi_k1_dav.get_nbasis(),
                    PARAM.inp.pw_diag_ndim,
                    this->diag_ethr,
                    david_maxiter,
                    false, //always do the subspace diag (check the implementation)
                    comm_info);

                auto hpsi_func = [pHamilt](
                    T* psi_in,
                    T* hpsi_out,
                    const int ld_psi,
                    const int nvec)
                    {
                        auto psi_iter_wrapper = psi::Psi<T, Device>(psi_in, 1, nvec, ld_psi, nullptr);
                        psi::Range bands_range(true, 0, 0, nvec-1);
                        using hpsi_info = typename hamilt::Operator<T, Device>::hpsi_info;
                        hpsi_info info(&psi_iter_wrapper, bands_range, hpsi_out);
                        pHamilt->ops->hPsi(info);
                    };
                auto subspace_func = [pHamilt](T* psi_out,
                    T* psi_in,
                    Real* eigenvalue_in_hsolver,
                    const int nband_in,
                    const int nbasis_max_in) {
                        // Convert "pointer data stucture" to a psi::Psi object
                        auto psi_in_wrapper = psi::Psi<T, Device>(psi_in, 1, nband_in, nbasis_max_in, nullptr);
                        auto psi_out_wrapper = psi::Psi<T, Device>(psi_out, 1, nband_in, nbasis_max_in, nullptr);

                        hsolver::DiagoIterAssist<T, Device>::diagH_subspace(pHamilt,
                            psi_in_wrapper,
                            psi_out_wrapper,
                            eigenvalue_in_hsolver,
                            nband_in);
                    };

                hsolver::DiagoIterAssist<T, Device>::avg_iter
                    += static_cast<double>(dav_subspace.diag(
                        hpsi_func, psi_k1_dav.get_pointer(),
                        psi_k1_dav.get_nbasis(),
                        eigenvalue.data(),
                        std::vector<bool>(psi_k1_dav.get_nbands(), true),
                        false /*scf*/));
            }
            else {throw std::runtime_error("HSolverLR::solve: method not implemented");}
        }

        // 5. copy eigenvalue to pes
        for (int ist = 0;ist < psi.get_nbands();++ist) { pes->ekb(ispin_solve, ist) = eigenvalue[ist];}


        // 6. output eigenvalues and eigenvectors
        print_eigs(eigenvalue, "eigenvalues: (Ry)");
        print_eigs(eigenvalue, "eigenvalues: (eV)", ModuleBase::Ry_to_eV);
        if (out_wfc_lr)
        {
            if (GlobalV::MY_RANK == 0)
            {
                std::ofstream ofs(PARAM.globalv.global_out_dir + "Excitation_Energy_" + spin_types[ispin_solve] + ".dat");
                ofs << std::setprecision(8) << std::scientific;
                for (auto& e : eigenvalue) {ofs << e << " ";}
                ofs.close();
            }
            LR_Util::write_psi_bandfirst(psi, PARAM.globalv.global_out_dir + "Excitation_Amplitude_" + spin_types[ispin_solve], GlobalV::MY_RANK);
        }

        // normalization is already satisfied
        // std::cout << "check normalization of eigenvectors:" << std::endl;
        // for (int ist = 0;ist < psi.get_nbands();++ist)
        // {
        //     double norm2 = 0;
        //     for (int ik = 0;ik < psi.get_nk();++ik)
        //     {
        //         for (int ib = 0;ib < psi.get_nbasis();++ib)
        //         {
        //             norm2 += square(psi(ist, ik, ib));
        //             // std::cout << "norm2_now=" << norm2 << std::endl;
        //         }
        //     }
        //     std::cout << "state " << ist << ", norm2=" << norm2 << std::endl;
        // }

        // output iters
        std::cout << "Average iterative diagonalization steps: " << hsolver::DiagoIterAssist<T, Device>::avg_iter
            << " ; where current threshold is: " << hsolver::DiagoIterAssist<T, Device>::PW_DIAG_THR << " . " << std::endl;
        // castmem_2d_2h_op()(cpu_ctx, cpu_ctx, pes->ekb.c, eigenvalues.data(), pes->ekb.nr * pes->ekb.nc);
    }
    template class HSolverLR<double>;
    template class HSolverLR<std::complex<double>>;
};
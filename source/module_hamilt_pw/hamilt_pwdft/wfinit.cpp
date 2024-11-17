#include "wfinit.h"

#include "module_base/macros.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi_initializer_atomic.h"
#include "module_psi/psi_initializer_atomic_random.h"
#include "module_psi/psi_initializer_nao.h"
#include "module_psi/psi_initializer_nao_random.h"
#include "module_psi/psi_initializer_random.h"
namespace psi
{

template <typename T, typename Device>
WFInit<T, Device>::WFInit(const std::string& init_wfc_in,
                          const std::string& ks_solver_in,
                          const std::string& basis_type_in,
                          const bool& use_psiinitializer_in,
                          wavefunc* p_wf_in,
                          ModulePW::PW_Basis_K* pw_wfc_in)
{
    this->init_wfc = init_wfc_in;
    this->ks_solver = ks_solver_in;
    this->basis_type = basis_type_in;
    this->use_psiinitializer = use_psiinitializer_in;
    this->p_wf = p_wf_in;
    this->pw_wfc = pw_wfc_in;
}

template <typename T, typename Device>
void WFInit<T, Device>::prepare_init(Structure_Factor* p_sf,
                                     UnitCell* p_ucell,
                                     const int& random_seed,
#ifdef __MPI
                                     Parallel_Kpoints* p_parak,
                                     const int& rank,
#endif
                                     pseudopot_cell_vnl* p_ppcell)
{
    if (!this->use_psiinitializer)
    {
        return;
    }
    // under restriction of C++11, std::unique_ptr can not be allocate via std::make_unique
    // use new instead, but will cause asymmetric allocation and deallocation, in literal aspect
    ModuleBase::timer::tick("WFInit", "prepare_init");
    if ((this->init_wfc.substr(0, 6) == "atomic") && (p_ucell->natomwfc == 0))
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_random<T, Device>());
    }
    else if (this->init_wfc == "atomic")
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_atomic<T, Device>());
    }
    else if (this->init_wfc == "random")
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_random<T, Device>());
    }
    else if (this->init_wfc == "nao")
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_nao<T, Device>());
    }
    else if (this->init_wfc == "atomic+random")
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_atomic_random<T, Device>());
    }
    else if (this->init_wfc == "nao+random")
    {
        this->psi_init = std::unique_ptr<psi_initializer<T, Device>>(new psi_initializer_nao_random<T, Device>());
    }
    else
    {
        ModuleBase::WARNING_QUIT("WFInit::prepare_init", "for new psi initializer, init_wfc type not supported");
    }

    //! function polymorphism is moved from constructor to function initialize.
    //! Two slightly different implementation are for MPI and serial case, respectively.
#ifdef __MPI
    this->psi_init->initialize(p_sf, pw_wfc, p_ucell, p_parak, random_seed, p_ppcell, rank);
#else
    this->psi_init->initialize(p_sf, pw_wfc, p_ucell, random_seed, p_ppcell);
#endif

    // always new->initialize->tabulate->allocate->proj_ao_onkG
    this->psi_init->tabulate();
    ModuleBase::timer::tick("WFInit", "prepare_init");
}

template <typename T, typename Device>
void WFInit<T, Device>::allocate_psi(Psi<std::complex<double>>*& psi,
                                     const int nkstot,
                                     const int nks,
                                     const int* ngk,
                                     const int npwx,
                                     Structure_Factor* p_sf)
{
    // allocate memory for std::complex<double> datatype psi
    // New psi initializer in ABACUS, Developer's note:
    // Because the calling relationship between WFInit and derived class is
    // complicated, up to upcoming of ABACUS 3.4, we only implement this new psi
    // initialization method for ksdft_pw, which means the routinely used dft theory.
    // For other theories like stochastic DFT, we still use the old method.

    // LCAOINPW also temporarily uses WFInit workflow, but in principle, it
    // should have its own ESolver. ESolver class is for controlling workflow for each
    // theory-basis combination, in the future it is also possible to seperate/decouple
    // the basis (representation) with operator (hamiltonian) and solver (diagonalization).
    // This feature requires feasible Linear Algebra library in-built in ABACUS, which
    // is not ready yet.
    if (this->use_psiinitializer) // new method
    {
        // psi_initializer drag initialization of pw wavefunction out of HSolver, make psi
        // initialization decoupled with HSolver (diagonalization) procedure.
        // However, due to EXX is hard to maintain, we still use the old method for EXX.
        // LCAOINPW in version >= 3.5.0 uses this new method.
        psi = this->psi_init->allocate();
    }
    else // old method
    {
        // old method explicitly requires variables such as total number of kpoints, number
        // of bands, number of G-vectors, and so on. Comparatively in new method, these
        // variables are imported in function called initialize.
        psi = this->p_wf->allocate(nkstot, nks, ngk, npwx);

        // however, init_at_1 does not actually initialize the psi, instead, it is a
        // function to calculate a interpolate table saving overlap intergral or say
        // Spherical Bessel Transform of atomic orbitals.
        this->p_wf->init_at_1(p_sf);
        // similarly, wfcinit not really initialize any wavefunction, instead, it initialize
        // the mapping from ixy, the 1d flattened index of point on fft grid (x, y) plane,
        // to the index of "stick", composed of grid points.
        this->p_wf->wfcinit(psi, pw_wfc);
    }
}

template <typename T, typename Device>
void WFInit<T, Device>::make_table(const int nks, Structure_Factor* p_sf)
{
    if (this->use_psiinitializer)
    {
    }    // do not need to do anything because the interpolate table is unchanged
    else // old initialization method, used in EXX calculation
    {
        this->p_wf->init_after_vc(nks); // reallocate wanf2, the planewave expansion of lcao
        this->p_wf->init_at_1(p_sf);    // re-calculate tab_at, the overlap matrix between atomic pswfc and jlq
    }
}

template <typename T, typename Device>
void WFInit<T, Device>::initialize_psi(Psi<std::complex<double>>* psi,
                                       psi::Psi<T, Device>* kspw_psi,
                                       hamilt::Hamilt<T, Device>* p_hamilt,
                                       std::ofstream& ofs_running)
{
    ModuleBase::timer::tick("WFInit", "initialize_psi");

    if (PARAM.inp.psi_initializer)
    {
        // if psig is not allocated before, allocate it
        if (!this->psi_init->psig_use_count())
        {
            this->psi_init->allocate(/*psig_only=*/true);
        }

        // loop over kpoints, make it possible to only allocate memory for psig at the only one kpt
        // like (1, nbands, npwx), in which npwx is the maximal npw of all kpoints
        for (int ik = 0; ik < this->pw_wfc->nks; ik++)
        {
            //! Fix the wavefunction to initialize at given kpoint
            psi->fix_k(ik);

            //! Update Hamiltonian from other kpoint to the given one
            p_hamilt->updateHk(ik);

            //! Project atomic orbitals on |k+G> planewave basis, where k is wavevector of kpoint
            //! and G is wavevector of the peroiodic part of the Bloch function
            this->psi_init->proj_ao_onkG(ik);

            //! psi_initializer manages memory of psig with shared pointer,
            //! its access to use is shared here via weak pointer
            //! therefore once the psi_initializer is destructed, psig will be destructed, too
            //! this way, we can avoid memory leak and undefined behavior
            std::weak_ptr<psi::Psi<T, Device>> psig = this->psi_init->share_psig();

            if (psig.expired())
            {
                ModuleBase::WARNING_QUIT("WFInit::initialize_psi", "psig lifetime is expired");
            }

            //! to use psig, we need to lock it to get a shared pointer version,
            //! then switch kpoint of psig to the given one
            auto psig_ = psig.lock();
            // CHANGE LOG: if not lcaoinpw, the psig will only be used in psi-initialization
            // so we can only allocate memory for one kpoint with the maximal number of pw
            // over all kpoints, then the memory space will be always enough. Then for each
            // kpoint, the psig is calculated in an overwrite manner.
            const int ik_psig = (psig_->get_nk() == 1) ? 0 : ik;
            psig_->fix_k(ik_psig);

            std::vector<typename GetTypeReal<T>::type> etatom(psig_->get_nbands(), 0.0);

            // then adjust dimension from psig to psi
            // either by matrix-multiplication or by copying-discarding
            if (this->psi_init->method() != "random")
            {
                // lcaoinpw and pw share the same esolver. In the future, we will have different esolver
                if (((this->ks_solver == "cg") || (this->ks_solver == "lapack")) && (this->basis_type == "pw"))
                {
                    // the following function is only run serially, to be improved
                    hsolver::DiagoIterAssist<T, Device>::diagH_subspace_init(p_hamilt,
                                                                             psig_->get_pointer(),
                                                                             psig_->get_nbands(),
                                                                             psig_->get_nbasis(),
                                                                             *(kspw_psi),
                                                                             etatom.data());
                    continue;
                }
                else if ((this->ks_solver == "lapack") && (this->basis_type == "lcao_in_pw"))
                {
                    if (ik == 0)
                    {
                        ofs_running << " START WAVEFUNCTION: LCAO_IN_PW, psi initialization skipped " << std::endl;
                    }
                    continue;
                }
                // else the case is davidson
            }
            else
            {
                if (this->ks_solver == "cg")
                {
                    hsolver::DiagoIterAssist<T, Device>::diagH_subspace(p_hamilt, *(psig_), *(kspw_psi), etatom.data());
                    continue;
                }
                // else the case is davidson
            }

            // for the Davidson method, we just copy the wavefunction (partially)
            for (int iband = 0; iband < kspw_psi->get_nbands(); iband++)
            {
                for (int ibasis = 0; ibasis < kspw_psi->get_nbasis(); ibasis++)
                {
                    (*(kspw_psi))(iband, ibasis) = (*psig_)(iband, ibasis);
                }
            }
        } // end k-point loop

        if (this->basis_type != "lcao_in_pw")
        {
            this->psi_init->deallocate_psig();
        }
    }
    else
    {
        // if (PARAM.inp.basis_type == "pw")
        // {
        //     for (int ik = 0; ik < this->pw_wfc->nks; ++ik)
        //     {
        //         //! Update Hamiltonian from other kpoint to the given one
        //         p_hamilt->updateHk(ik);

        //         //! Fix the wavefunction to initialize at given kpoint
        //         kspw_psi->fix_k(ik);

        //         /// for psi init guess!!!!
        //         hamilt::diago_PAO_in_pw_k2(this->ctx, ik, *kspw_psi, this->pw_wfc, this->p_wf, p_hamilt);
        //     }
        // }
    }

    ModuleBase::timer::tick("WFInit", "initialize_psi");
}

template class WFInit<std::complex<float>, base_device::DEVICE_CPU>;
template class WFInit<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class WFInit<std::complex<float>, base_device::DEVICE_GPU>;
template class WFInit<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace psi
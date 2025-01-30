#include "psi_init.h"

#include "module_base/macros.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_quit.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi_initializer_atomic.h"
#include "module_psi/psi_initializer_atomic_random.h"
#include "module_psi/psi_initializer_file.h"
#include "module_psi/psi_initializer_nao.h"
#include "module_psi/psi_initializer_nao_random.h"
#include "module_psi/psi_initializer_random.h"
namespace psi
{

template <typename T, typename Device>
PSIInit<T, Device>::PSIInit(const std::string& init_wfc_in,
                            const std::string& ks_solver_in,
                            const std::string& basis_type_in,
                            const int& rank_in,
                            const UnitCell& ucell_in,
                            const Structure_Factor& sf_in,
                            const K_Vectors& kv_in,
                            const pseudopot_cell_vnl& nlpp_in,
                            const ModulePW::PW_Basis_K& pw_wfc_in)
    : ucell(ucell_in), sf(sf_in), nlpp(nlpp_in), kv(kv_in), pw_wfc(pw_wfc_in), rank(rank_in)
{
    this->init_wfc = init_wfc_in;
    this->ks_solver = ks_solver_in;
    this->basis_type = basis_type_in;
}

template <typename T, typename Device>
void PSIInit<T, Device>::prepare_init(const int& random_seed)
{

    // under restriction of C++11, std::unique_ptr can not be allocate via std::make_unique
    // use new instead, but will cause asymmetric allocation and deallocation, in literal aspect
    ModuleBase::timer::tick("PSIInit", "prepare_init");
    this->psi_initer.reset();
    if (this->init_wfc == "random")
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_random<T>());
    }
    else if (this->init_wfc == "file")
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_file<T>());
    }
    else if ((this->init_wfc.substr(0, 6) == "atomic") && (this->ucell.natomwfc == 0))
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_random<T>());
    }
    else if (this->init_wfc == "atomic"
             || (this->init_wfc == "atomic+random" && this->ucell.natomwfc < PARAM.inp.nbands))
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_atomic<T>());
    }
    else if (this->init_wfc == "atomic+random")
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_atomic_random<T>());
    }
    else if (this->init_wfc == "nao")
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_nao<T>());
    }
    else if (this->init_wfc == "nao+random")
    {
        this->psi_initer = std::unique_ptr<psi_initializer<T>>(new psi_initializer_nao_random<T>());
    }
    else
    {
        ModuleBase::WARNING_QUIT("PSIInit::prepare_init", "for new psi initializer, init_wfc type not supported");
    }

    this->psi_initer->initialize(&sf, &pw_wfc, &ucell, &kv, random_seed, &nlpp, rank);
    this->psi_initer->tabulate();

    ModuleBase::timer::tick("PSIInit", "prepare_init");
}

template <typename T, typename Device>
void PSIInit<T, Device>::initialize_psi(Psi<std::complex<double>>* psi,
                                        psi::Psi<T, Device>* kspw_psi,
                                        hamilt::Hamilt<T, Device>* p_hamilt,
                                        std::ofstream& ofs_running)
{
    if (kspw_psi->get_nbands() == 0 || GlobalV::MY_STOGROUP != 0)
    {
        return;
    }
    if (this->basis_type == "lcao_in_pw")
    {
        return;
    }
    ModuleBase::timer::tick("PSIInit", "initialize_psi");

    const int nbands_start = this->psi_initer->nbands_start();
    const int nbands = psi->get_nbands();
    const int nbasis = psi->get_nbasis();
    const bool not_equal = (nbands_start != nbands);

    Psi<T>* psi_cpu = reinterpret_cast<psi::Psi<T>*>(psi);
    Psi<T, Device>* psi_device = kspw_psi;

    if (not_equal)
    {
        psi_cpu = new Psi<T>(1, nbands_start, nbasis, nbasis, true);
        psi_device = PARAM.inp.device == "gpu" ? new psi::Psi<T, Device>(psi_cpu[0])
                                               : reinterpret_cast<psi::Psi<T, Device>*>(psi_cpu);
    }
    else if (PARAM.inp.precision == "single")
    {
        if (PARAM.inp.device == "cpu")
        {
            psi_cpu = reinterpret_cast<psi::Psi<T>*>(kspw_psi);
            psi_device = kspw_psi;
        }
        else
        {
            psi_cpu = new Psi<T>(1, nbands_start, nbasis, nbasis, true);
            psi_device = kspw_psi;  
        }
    }

    // loop over kpoints, make it possible to only allocate memory for psig at the only one kpt
    // like (1, nbands, npwx), in which npwx is the maximal npw of all kpoints
    for (int ik = 0; ik < this->pw_wfc.nks; ik++)
    {
        //! Fix the wavefunction to initialize at given kpoint
        psi->fix_k(ik);
        kspw_psi->fix_k(ik);

        //! Update Hamiltonian from other kpoint to the given one
        p_hamilt->updateHk(ik);

        //! initialize psi_cpu
        this->psi_initer->init_psig(psi_cpu->get_pointer(), ik);
        if (psi_device->get_pointer() != psi_cpu->get_pointer())
        {
            syncmem_h2d_op()(psi_device->get_pointer(), psi_cpu->get_pointer(), nbands_start * nbasis);
        }

        std::vector<typename GetTypeReal<T>::type> etatom(nbands_start, 0.0);

        if (this->ks_solver == "cg")
        {
            if (not_equal)
            {
                // for diagH_subspace_init, psi_device->get_pointer() and kspw_psi->get_pointer() should be different
                hsolver::DiagoIterAssist<T, Device>::diagH_subspace_init(p_hamilt,
                                                                         psi_device->get_pointer(),
                                                                         nbands_start,
                                                                         nbasis,
                                                                         *(kspw_psi),
                                                                         etatom.data());
            }
            else
            {
                // for diagH_subspace, psi_device->get_pointer() and kspw_psi->get_pointer() can be the same
                hsolver::DiagoIterAssist<T, Device>::diagH_subspace(p_hamilt,
                                                                    *psi_device,
                                                                    *kspw_psi,
                                                                    etatom.data(),
                                                                    nbands_start);
            }
        }
        else // dav, bpcg
        {
            if (psi_device->get_pointer() != kspw_psi->get_pointer())
            {
                syncmem_complex_op()(kspw_psi->get_pointer(), psi_device->get_pointer(), nbands * nbasis);
            }
        }
    } // end k-point loop

    if (not_equal)
    {
        delete psi_cpu;
        if(PARAM.inp.device == "gpu")
        {
            delete psi_device;
        }
    }
    else if (PARAM.inp.precision == "single" && PARAM.inp.device == "gpu")
    {
        delete psi_cpu;
    }

    ModuleBase::timer::tick("PSIInit", "initialize_psi");
}

template <typename T, typename Device>
void PSIInit<T, Device>::initialize_lcao_in_pw(Psi<T>* psi_local, std::ofstream& ofs_running)
{
    ofs_running << " START WAVEFUNCTION: LCAO_IN_PW, psi initialization skipped " << std::endl;
    assert(this->psi_initer->method() == "nao");
    for (int ik = 0; ik < this->pw_wfc.nks; ik++)
    {
        psi_local->fix_k(ik);
        this->psi_initer->init_psig(psi_local->get_pointer(), ik);
    }
}

void allocate_psi(Psi<std::complex<double>>*& psi, const int& nks, const std::vector<int>& ngk, const int& nbands, const int& npwx)
{
    assert(npwx > 0);
    assert(nks > 0);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "npwx", npwx);

    delete psi;
    int nks2 = nks;
    if (PARAM.inp.calculation == "nscf" && PARAM.inp.mem_saver == 1)
    {
        nks2 = 1;
    }
    psi = new psi::Psi<std::complex<double>>(nks2, nbands, npwx * PARAM.globalv.npol, ngk, true);
    const size_t memory_cost = sizeof(std::complex<double>) * nks2 * nbands * (PARAM.globalv.npol * npwx);
    std::cout << " MEMORY FOR PSI (MB)  : " << static_cast<double>(memory_cost) / 1024.0 / 1024.0 << std::endl;
    ModuleBase::Memory::record("Psi_PW", memory_cost);
}

template class PSIInit<std::complex<float>, base_device::DEVICE_CPU>;
template class PSIInit<std::complex<double>, base_device::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class PSIInit<std::complex<float>, base_device::DEVICE_GPU>;
template class PSIInit<std::complex<double>, base_device::DEVICE_GPU>;
#endif
} // namespace psi
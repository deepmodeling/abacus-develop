#include "source_psi/setup_psi.h"
#include "source_io/module_parameter/parameter.h" // use parameter
#include "source_io/read_wfc_nao.h"

template <typename T>
Setup_Psi<T>::Setup_Psi(){}

template <typename T>
Setup_Psi<T>::~Setup_Psi(){}

// the size of psi may change during scf (in the future we can support such
// calculations: first SZP or DZ, and then DZP calculations.
// In that case, psi may change its size multiple times during SCF
template <typename T>
void Setup_Psi<T>::allocate_psi(
		psi::Psi<T>* &psi,
		ModuleBase::matrix &ekb,
		ModuleBase::matrix &wg,
		const K_Vectors &kv,
        const Parallel_Orbitals &para_orb,
		const Input_para &inp)
{
    // init electronic wave function psi
    if (psi == nullptr)
    {
        int nsk = 0;
        int ncol = 0;
        if (PARAM.globalv.gamma_only_local)
        {
            nsk = inp.nspin;
            ncol = para_orb.ncol_bands;
            if (inp.ks_solver == "genelpa" || inp.ks_solver == "elpa" || inp.ks_solver == "lapack"
                || inp.ks_solver == "pexsi" || inp.ks_solver == "cusolver"
                || inp.ks_solver == "cusolvermp")
            {
                ncol = para_orb.ncol;
            }
        }
        else
        {
            nsk = kv.get_nks();
#ifdef __MPI
            ncol = para_orb.ncol_bands;
#else
            ncol = inp.nbands;
#endif
        }
        psi = new psi::Psi<T>(nsk, ncol, para_orb.nrow, kv.ngk, true);
    }

    //------------------------------
    //! read psi from file
    //------------------------------
    if (inp.init_wfc == "file" && inp.esolver_type != "tddft")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
             para_orb, *psi, ekb, wg, kv.ik2iktot,
             kv.get_nkstot(), inp.nspin))
        {
            ModuleBase::WARNING_QUIT("Setup_Psi", "read electronic wave functions failed");
        }
    }
}


template <typename T>
void Setup_Psi<T>::deallocate_psi(psi::Psi<T>* &psi)
{
	if(psi!=nullptr)
	{
		delete psi;
	}
}

template class Setup_Psi<double>;
template class Setup_Psi<std::complex<double>>;

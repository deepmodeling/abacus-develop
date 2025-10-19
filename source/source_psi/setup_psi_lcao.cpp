#include "source_psi/setup_psi_lcao.h"

template <typename T>
Setup_Psi_LCAO<T>::Setup_Psi_LCAO(){}

template <typename T>
Setup_Psi_LCAO<T>::~Setup_Psi_LCAO(){}

template <typename T>
void Setup_Psi_LCAO<T>::before_runner(
		psi::Psi<T>* &psi,
		matrix &ekb,
		matrix &wg,
		const K_Vectors &kv,
        const Parallel_Orbitals &para_orb,
		const Input_para &inp)
{
    // 5) init electronic wave function psi
    if (this->psi == nullptr)
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
        this->psi = new psi::Psi<TK>(nsk, ncol, para_orb.nrow, kv.ngk, true);
    }

    // 6) read psi from file
    if (inp.init_wfc == "file" && inp.esolver_type != "tddft")
    {
        if (!ModuleIO::read_wfc_nao(PARAM.globalv.global_readin_dir,
             para_orb, *psi, ekb, wg, kv.ik2iktot,
             kv.get_nkstot(), inp.nspin))
        {
            ModuleBase::WARNING_QUIT("ESolver_KS_LCAO", "read electronic wave functions failed");
        }
    }
}


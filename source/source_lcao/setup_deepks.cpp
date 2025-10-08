#include "source_lcao/setup_deepks.h"
#include "source_lcao/LCAO_domain.h"

template <typename TK>
Setup_DeePKS<TK>::Setup_DeePKS(){}

template <typename TK>
Setup_DeePKS<TK>::~Setup_DeePKS(){}

template <typename TK>
void Setup_DeePKS<TK>::before_runner(const UnitCell& ucell, // unitcell
	const int nks, // number of k points
    const LCAO_Orbitals &orb, // orbital info
	Parallel_Orbitals &pv, // parallel orbitals
	const Input_para &inp)
{
#ifdef __MLALGO
    LCAO_domain::DeePKS_init(ucell, pv, nks, orb, this->ld, GlobalV::ofs_running);
    if (inp.deepks_scf)
    {
        // load the DeePKS model from deep neural network
        DeePKS_domain::load_model(inp.deepks_model, this->ld.model_deepks);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        DeePKS_domain::read_pdm((inp.init_chg == "file"), inp.deepks_equiv,
          this->ld.init_pdm, ucell.nat, orb.Alpha[0].getTotal_nchi() * ucell.nat,
          this->ld.lmaxd, this->ld.inl2l, *orb.Alpha, this->ld.pdm);
    }
#endif
}

template class Setup_DeePKS<double>;
template class Setup_DeePKS<std::complex<double>>;

#include "source_lcao/setup_deepks.h"

template <typename TK>
DeePKS<TK>::DeePKS(){}

template <typename TK>
DeePKS<TK>::~DeePKS(){}

template <typename TK>
void DeePKS<TK>::before_runner(UnitCell& ucell, // unitcell
	K_Vectors &kv, // k points
    const LCAO_Orbitals &orb, // orbital info
	const Parallel_Orbitals &pv, // parallel orbitals
	const Input_para& inp)
{
#ifdef __MLALGO
    LCAO_domain::DeePKS_init(ucell, pv, kv.get_nks(), orb_, this->ld, GlobalV::ofs_running);
    if (inp.deepks_scf)
    {
        // load the DeePKS model from deep neural network
        DeePKS_domain::load_model(inp.deepks_model, ld.model_deepks);
        // read pdm from file for NSCF or SCF-restart, do it only once in whole calculation
        DeePKS_domain::read_pdm((inp.init_chg == "file"), inp.deepks_equiv,
          ld.init_pdm, ucell.nat, orb_.Alpha[0].getTotal_nchi() * ucell.nat,
          ld.lmaxd, ld.inl2l, *orb_.Alpha, ld.pdm);
    }
#endif
}

template class DeePKS<double>;
template class DeePKS<std::complex<double>>;

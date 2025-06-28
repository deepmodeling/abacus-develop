#include "esolver_ks_lcao.h"
#include "source_io/ctrl_output_lcao.h"

namespace ModuleESolver
{

template <typename TK, typename TR>
void ESolver_KS_LCAO<TK, TR>::after_scf(UnitCell& ucell, const int istep, const bool conv_esolver)
{
    ModuleBase::TITLE("ESolver_KS_LCAO", "after_scf");
    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");

    //------------------------------------------------------------------
    //! 1) call after_scf() of ESolver_KS
    //------------------------------------------------------------------
    ESolver_KS<TK>::after_scf(ucell, istep, conv_esolver);

    //------------------------------------------------------------------
    //! 2) output of lcao every few ionic steps 
    //------------------------------------------------------------------
	auto* estate = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(this->pelec);
    auto* hamilt_lcao = dynamic_cast<hamilt::HamiltLCAO<TK, TR>*>(this->p_hamilt);

	if(!estate)
	{
		ModuleBase::WARNING_QUIT("ModuleIO::ctrl_output_lcao","pelec does not exist");
	}

    if(istep % PARAM.inp.out_interval == 0)
    {

        ModuleIO::ctrl_output_lcao<TK, TR>(ucell, 
				this->kv,
				estate, 
				this->pv, 
				this->gd,
				this->psi,
				hamilt_lcao,
				this->two_center_bundle_,
				this->GK,
				this->orb_,
				this->pw_wfc,
				this->pw_rho,
				this->GridT,
				this->pw_big,
				this->sf,
				this->rdmft_solver,
#ifdef __MLALGO
				this->ld,
#endif
#ifdef __EXX
				*this->exd,
				*this->exc,
#endif
				istep);
    }

    //------------------------------------------------------------------
    //! 3) Clean up RA, which is used to serach for adjacent atoms
    //------------------------------------------------------------------
    if (!PARAM.inp.cal_force && !PARAM.inp.cal_stress)
    {
        RA.delete_grid();
    }

    ModuleBase::timer::tick("ESolver_KS_LCAO", "after_scf");
}

template class ESolver_KS_LCAO<double, double>;
template class ESolver_KS_LCAO<std::complex<double>, double>;
template class ESolver_KS_LCAO<std::complex<double>, std::complex<double>>;
} // namespace ModuleESolver

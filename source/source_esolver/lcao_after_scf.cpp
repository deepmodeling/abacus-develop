#include "esolver_ks_lcao.h"
#include "source_io/ctrl_scf_lcao.h"

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
		ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","pelec does not exist");
	}

	if(!hamilt_lcao)
	{
		ModuleBase::WARNING_QUIT("ESolver_KS_LCAO::after_scf","p_hamilt does not exist");
	}


    //*****
    // if istep_in = -1, istep will not appear in file name
    // if iter_in = -1, iter will not appear in file name
	int istep_in = -1;
    int iter_in = -1;
	bool out_flag = false;
	if (PARAM.inp.out_freq_ion>0) // default value of out_freq_ion is 0
	{
		if (istep % PARAM.inp.out_freq_ion == 0)
		{
			istep_in = istep;
			out_flag = true;
		}
	}
	else if(conv_esolver || this->scf_nmax_flag) // mohan add scf_nmax_flag on 20250921
	{
		out_flag = true;
	}
    //*****

	if (out_flag)
	{
		ModuleIO::ctrl_scf_lcao<TK, TR>(ucell,
		  PARAM.inp, this->kv, estate, this->pv, 
		  this->gd, this->psi, hamilt_lcao,
          this->two_center_bundle_, this->GK,
          this->orb_, this->pw_wfc, this->pw_rho,
          this->GridT, this->pw_big, this->sf,
		  this->rdmft_solver,
          this->deepks,
		  this->exx_nao,
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

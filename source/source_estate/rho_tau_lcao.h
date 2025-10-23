#ifndef RHO_TAU_LCAO_H
#define RHO_TAU_LCAO_H

// generate charge density from different basis or methods
namespace elecstate
{
    // interface for HSolver to calculate rho from Psi
	void psi2rho_lcao(const psi::Psi<std::complex<double>>& psi, 
			std::vector<hamilt::HContainer<double>*> &dmr,
			const int nspin,
			Charge* chg);

  
}

#endif

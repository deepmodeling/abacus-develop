#ifndef RHO_TAU_LCAO_H
#define RHO_TAU_LCAO_H

#include "source_lcao/module_hcontainer/hcontainer.h"
#include "source_estate/module_charge/charge.h"

// generate charge density from different basis or methods
namespace LCAO_domain
{
	void dm2rho(std::vector<hamilt::HContainer<double>*> &dmr,
			const int nspin,
			Charge* chr);

	void dm2tau(std::vector<hamilt::HContainer<double>*> &dmr,
			const int nspin,
			Charge* chr);
}

#endif

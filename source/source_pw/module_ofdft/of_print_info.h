#ifndef OF_PRINT_INFO_H
#define OF_PRINT_INFO_H

#include "source_estate/elecstate.h" // electronic states
#include "source_pw/module_ofdft/kedf_manager.h"

#include "source_base/timer_wrapper.h"


namespace OFDFT
{

void print_info(const int iter,
	ModuleBase::TimePoint &iter_time,
	const double &energy_current,
	const double &energy_last,
	const double &normdLdphi,
	const elecstate::ElecState *pelec,
	KEDF_Manager *kedf_manager,
	const bool conv_esolver);

}

#endif



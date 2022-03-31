#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H
#include<vector>
#include<module_base/matrix.h>
#include<module_base/complexmatrix.h>
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/gint_gamma.h"

class IState_Charge
{
public:
    IState_Charge(std::vector<ModuleBase::matrix> &wfc_gamma_in,
        Local_Orbital_Charge &loc_in);
    ~IState_Charge();

	void begin(Gint_Gamma &gg);

private:

	int *bands_picked;

#ifdef __MPI
	void idmatrix(const int &ib);
#endif
    std::vector<ModuleBase::matrix>* wfc_gamma;
    Local_Orbital_Charge* loc;

};
#endif

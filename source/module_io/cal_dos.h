#ifndef DOS_H
#define DOS_H
#include<vector>
#include "module_base/matrix.h"

namespace ModuleIO
{
	bool calculate_dos(const int &is,		
		const std::string &fn,// file address for DOS.
		const std::string &fn1,// file address for DOS_smearing.
		const double &de_ev, // delta energy in ev.
		const double &emax_ev,// maximal energy in ev.
		const double &emin_ev,// minimal energy in ev.
		const double &bcoeff,
		const int &nks,//number of k points
		const int &nkstot,
		const std::vector<double> &wk,//weight of k points
		const std::vector<int> &isk,
		const int &nbands,// number of bands
		const ModuleBase::matrix &ekb, //store energy for each k point and each band
		const ModuleBase::matrix &wg); //weight of (kpoint,bands))
}

#endif 

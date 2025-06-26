#ifndef DOS_H
#define DOS_H

#include<vector>
#include "source_base/matrix.h"
#include "source_estate/fp_energy.h"

namespace ModuleIO
{

	void prepare_dos(std::ofstream& ofs_running,
			const elecstate::efermi &energy_fermi,
			const ModuleBase::matrix& ekb,
			const int nks,
			const int nbands,
			const double& dos_edelta_ev,
            const double& dos_scale,
			double &emax,
            double &emin);

	bool cal_dos(const int &is,		
		const std::string &fn,// file address for DOS.
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

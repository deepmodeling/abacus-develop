#ifndef DOS_TEST_H
#define DOS_TEST_H

#include <string>
#include <vector>

#include "source_base/constants.h"
#include "source_base/matrix.h"

#include "./dos_test_data.h"

class DosPrepare
{
public:
	DosPrepare(int is_in,
			std::string fa_in,
			double de_ev_in,
			double emax_ev_in,
			double emin_ev_in,
			double bcoeff_in,
			int nks_in,
			int nkstot_in,
			int nbands_in):
		is(is_in),fa(fa_in),
		de_ev(de_ev_in),emax_ev(emax_ev_in),
		emin_ev(emin_ev_in),bcoeff(bcoeff_in),
		nks(nks_in),nkstot(nkstot_in),nbands(nbands_in){}

	int is;
	std::string fa;
	double de_ev;
	double emax_ev;
	double emin_ev;
	double bcoeff;
	int nks;
	int nkstot;
	int nbands;
	std::vector<int> isk;
	std::vector<double> wk;
	ModuleBase::matrix ekb;
	ModuleBase::matrix wg;

	void set_isk()
	{
		this->isk.resize(nks, 0); // spin-unpolarized case, only 1 spin
	}

	void set_wk()
	{
		this->wk = dos_test_data::wk;
	}

	void set_istate_info()
	{
		this->ekb.create(nks, nbands);
		this->wg.create(nks, nbands);
		for (int ik = 0; ik < nks; ++ik)
		{
			for (int ib = 0; ib < nbands; ++ib)
			{
				this->ekb(ik, ib) = dos_test_data::ekb_ev[ik * nbands + ib];
				this->wg(ik, ib) = dos_test_data::wg[ik * nbands + ib];
			}
		}
		this->ekb *= 1.0 / ModuleBase::Ry_to_eV; // eV -> Ry
	}
};

#endif

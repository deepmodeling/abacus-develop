#include "write_dos_pw.h"

#include "cal_dos.h"
#include "module_base/parallel_reduce.h"
#include "module_io/input.h"

void ModuleIO::write_dos_pw(const ModuleBase::matrix &ekb,
	const ModuleBase::matrix &wg,
	const K_Vectors& kv,
	const double &dos_edelta_ev,
	const double &dos_scale,
	const double &bcoeff)
{
	ModuleBase::TITLE("ModuleIO","write_dos_pw");
	
	int nspin0=1;
	if(GlobalV::NSPIN==2) nspin0=2;


	//find energy range
	double emax = ekb(0, 0);
	double emin = ekb(0, 0);
	for(int ik=0; ik<kv.nks; ++ik)
	{
		for(int ib=0; ib<GlobalV::NBANDS; ++ib)
		{
			emax = std::max( emax, ekb(ik, ib) );
			emin = std::min( emin, ekb(ik, ib) );
		}
	}

#ifdef __MPI
	Parallel_Reduce::gather_max_double_all(emax);
	Parallel_Reduce::gather_min_double_all(emin);
#endif

	emax *= ModuleBase::Ry_to_eV;
	emin *= ModuleBase::Ry_to_eV;

	if(INPUT.dos_setemax)	emax = INPUT.dos_emax_ev;
	if(INPUT.dos_setemin)	emin = INPUT.dos_emin_ev;
	if(!INPUT.dos_setemax && !INPUT.dos_setemin)
	{
		//scale up a little bit so the end peaks are displaced better
		double delta=(emax-emin)*dos_scale;
		emax=emax+delta/2.0;
		emin=emin-delta/2.0;
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
	// 		atom_arrange::set_sr_NL();
	//		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21
	//determine #. energy points	
	const double de_ev = dos_edelta_ev;
	//std::cout << de_ev;

	const int npoints = static_cast<int>(std::floor ( ( emax - emin ) / de_ev ));
		const int np=npoints;

	for(int is=0; is<nspin0; ++is)
	{
		//DOS_ispin contains not smoothed dos
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "DOS" << is+1;
		std::stringstream ss1;
		ss1 << GlobalV::global_out_dir << "DOS" << is+1 << "_smearing.dat";
		ModuleIO::calculate_dos(is,
			ss.str(),
			ss1.str(),
			dos_edelta_ev,
			emax,
			emin,
			bcoeff,
			kv.nks,
			kv.nkstot,
			kv.wk,
			kv.isk,
			GlobalV::NBANDS,
			ekb,
			wg);
	}


}


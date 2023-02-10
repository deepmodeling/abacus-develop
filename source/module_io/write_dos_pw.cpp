#include "dos.h"
#include "write_dos_pw.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_elecstate/energy.h"
#include "src_parallel/parallel_reduce.h"
#include "module_elecstate/elecstate.h"

void ModuleIO::write_dos_pw(const elecstate::ElecState* pelec, 
		const int &out_dos, 
		const int &out_band, 
		const double &dos_edelta_ev,
		const double &dos_scale,
		const double &ef) 
{
	ModuleBase::TITLE("ModuleIO","write_dos_pw");

	if(out_dos !=0 || out_band !=0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }	
	
	int nspin0=1;
	if(GlobalV::NSPIN==2) nspin0=2;

	if(out_dos)
	{
//find energy range
		double emax = pelec->ekb(0, 0);
		double emin = pelec->ekb(0, 0);
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				emax = std::max( emax, pelec->ekb(ik, ib) );
				emin = std::min( emin, pelec->ekb(ik, ib) );
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

			 ModuleIO::calculate_dos(
					 is,
					 GlobalC::kv.isk,
					 ss.str(),
					 ss1.str(), 
					 dos_edelta_ev, 
					 emax, 
					 emin, 
					 GlobalC::kv.nks, GlobalC::kv.nkstot, GlobalC::kv.wk, pelec->wg, GlobalV::NBANDS, pelec->ekb );
	 	}

	}//out_dos=1
	if(out_band) //pengfei 2014-10-13
	{
		int nks=0;
		if(nspin0==1) 
		{
			nks = GlobalC::kv.nkstot;
		}
		else if(nspin0==2) 
		{
			nks = GlobalC::kv.nkstot/2;
		}

		for(int is=0; is<nspin0; is++)
		{
			std::stringstream ss2;
			ss2 << GlobalV::global_out_dir << "BANDS_" << is+1 << ".dat";
			GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
			ModuleIO::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, ef*0, pelec->ekb);
		}

	}
}


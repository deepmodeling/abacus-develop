#include "nscf_band.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"

void ModuleIO::nscf_band(
	const int &is,
	const std::string &out_band_dir, 
	const int &nks, 
	const int &nband,
	const double &fermie,
	const ModuleBase::matrix& ekb,
	const K_Vectors& kv,
	const Parallel_Kpoints* Pkpoints)
{
	ModuleBase::TITLE("ModuleIO","nscf_band");
	ModuleBase::timer::tick("ModuleIO", "nscf_band");

#ifdef __MPI
	if(GlobalV::MY_RANK==0)
	{
		std::ofstream ofs(out_band_dir.c_str());//make the file clear!!
		ofs.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<double> klength;
	klength.resize(nks);
	klength[0] = 0.0;
	for(int ik=0; ik<nks; ik++)
	{
		if (ik>0)
		{
			auto delta=kv.kvec_c[ik]-kv.kvec_c[ik-1];
			klength[ik] = klength[ik-1] + delta.norm();
		}
		if ( GlobalV::MY_POOL == Pkpoints->whichpool[ik] )
		{
			const int ik_now = ik - Pkpoints->startk_pool[GlobalV::MY_POOL];
			if( kv.isk[ik_now+is*nks] == is )
			{ 
				if ( GlobalV::RANK_IN_POOL == 0)
				{
					std::ofstream ofs(out_band_dir.c_str(),std::ios::app);
					ofs << std::setprecision(8);
					//start from 1
					ofs << ik+1;
					ofs << " " << klength[ik] << " ";
					for(int ib = 0; ib < nband; ib++)
					{
						ofs << " " << (ekb(ik_now+is*nks, ib)-fermie) * ModuleBase::Ry_to_eV;
					}
					ofs << std::endl;
					ofs.close();	
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// old version
	/*
	for(int ip=0;ip<GlobalV::KPAR;ip++)
	{
		if(GlobalV::MY_POOL == ip && GlobalV::RANK_IN_POOL == 0)
		{
			std::ofstream ofs(out_band_dir.c_str(),ios::app);
			for(int ik=0;ik<nks;ik++)
			{
				ofs<<std::setw(12)<<ik;
				for(int ib = 0; ib < nband; ib++)
				{
					ofs <<std::setw(12)<< ekb[ik][ib] * ModuleBase::Ry_to_eV;
				}
				ofs<<std::endl;
			}
			ofs.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	*/
#else
//	std::cout<<"\n nband = "<<nband<<std::endl;
//	std::cout<<out_band_dir<<std::endl;

	std::ofstream ofs(out_band_dir.c_str());
	for(int ik=0;ik<nks;ik++)
	{
		if( kv.isk[ik] == is)
		{
			ofs<<std::setw(12)<<ik + 1;
			for(int ibnd = 0; ibnd < nband; ibnd++)
			{
				ofs <<std::setw(15) << (ekb(ik, ibnd)-fermie) * ModuleBase::Ry_to_eV;
			}
			ofs<<std::endl;
		}
	}
	ofs.close();
#endif

	ModuleBase::timer::tick("ModuleIO", "nscf_band");
	return;
}

#include "nscf_fermi_surf.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

void ModuleIO::nscf_fermi_surface(const std::string &out_band_dir,
	const int &nband,
	const double &ef,
	const K_Vectors& kv,
	const Parallel_Kpoints& Pkpoints,
	const UnitCell& ucell,
	const ModuleBase::matrix &ekb)
{
	ModuleBase::TITLE("ModuleIO","nscf_fermi_surface");
	ModuleBase::timer::tick("ModuleIO", "nscf_fermi_surface");
#ifdef __MPI

	int start = 1;
	int end = GlobalV::NBANDS;

	std::ofstream ofs;
	if(GlobalV::MY_RANK==0)
	{
		ofs.open(out_band_dir.c_str());//make the file clear!!
		ofs << std::setprecision(6);
		ofs.close();	
	}

	for(int ik=0; ik<kv.nkstot; ik++)
	{
		if ( GlobalV::MY_POOL == Pkpoints.whichpool[ik] )
		{
			if( GlobalV::RANK_IN_POOL == 0)
			{
				std::ofstream ofs(out_band_dir.c_str(),std::ios::app);
				ofs << std::setprecision(8);

				if(ik==0)
				{
					ofs << " BEGIN_INFO" << std::endl;
					ofs << "   #" << std::endl;
					ofs << "   # this is a Band-XCRYSDEN-Structure-File" << std::endl;
					ofs << "   # aimed at Visualization of Fermi Surface" << std::endl;
					ofs << "   #" << std::endl;
					ofs << "   # Case: " << ucell.latName << std::endl;
					ofs << "   #" << std::endl;	
					ofs << " Fermi Energy: " << ef << std::endl;
					ofs << " END_INFO" << std::endl;
					ofs << " BEGIN_BLOCK_BANDGRID_3D" << std::endl;
					ofs << " band_energies" << std::endl;
					ofs << " BANDGRID_3D_BANDS" << std::endl;
					ofs << " " << end-start+1 << std::endl;
					ofs << " NKX NKY NKZ" << std::endl;
					ofs << " 0 0 0" << std::endl;
					ofs << " " << ucell.G.e11 << " " << ucell.G.e12 << " " << ucell.G.e13 << std::endl; 
					ofs << " " << ucell.G.e21 << " " << ucell.G.e22 << " " << ucell.G.e23 << std::endl; 
					ofs << " " << ucell.G.e31 << " " << ucell.G.e32 << " " << ucell.G.e33 << std::endl; 
				}

				const int ik_now = ik - Pkpoints.startk_pool[GlobalV::MY_POOL];
				ofs << "ik= " << ik << std::endl;
				ofs << kv.kvec_c[ik_now].x << " " << kv.kvec_c[ik_now].y << " " << kv.kvec_c[ik_now].z << std::endl;  

				for(int ib = 0; ib < nband; ib++)
				{
					ofs << " " << ekb(ik_now, ib) * ModuleBase::Ry_to_eV;
				}
				ofs << std::endl;

				// the last k point
				if(ik==kv.nkstot-1)
				{
					ofs << " END_BANDGRID_3D" << std::endl;
					ofs << " END_BLOCK_BANDGRID_3D" << std::endl;
				}
				ofs.close();

			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

#else


#endif
	ModuleBase::timer::tick("ModuleIO", "nscf_fermi_surface");
	return;
}

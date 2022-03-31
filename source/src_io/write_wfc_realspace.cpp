//======================
// AUTHOR : Peize Lin
// DATE :   2021-11-21
//======================

#include "write_wfc_realspace.h"
#include "src_pw/global.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include <fstream>
#include <stdexcept>
#include <cstdlib>

namespace Write_Wfc_Realspace
{
	// write ||wfc_r|| for all k-points and all bands
	// Input: wfc_g[ik](ib,ig)
	// loop order is for(z){for(y){for(x)}}
	void write_wfc_realspace_1(const ModuleBase::ComplexMatrix*const wfc_g, const std::string &folder_name, const bool& square)
	{
		ModuleBase::TITLE("Write_Wfc_Realspace", "write_wfc_realspace_1");
		ModuleBase::timer::tick("Write_Wfc_Realspace", "write_wfc_realspace_1");

		const string outdir = GlobalV::global_out_dir + folder_name + "/";
		const std::string command0 =  "test -d " + outdir + " || mkdir " + outdir;
    	if(GlobalV::MY_RANK==0)
        	system( command0.c_str() );

#ifdef __MPI		
		std::vector<MPI_Request> mpi_requests;
#endif
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			const int ik_out = (GlobalV::NSPIN!=2)
				? ik + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL]
				: ik - GlobalC::kv.nks/2*GlobalC::kv.isk[ik] + GlobalC::kv.nkstot/2*GlobalC::kv.isk[ik] + GlobalC::Pkpoints.startk_pool[GlobalV::MY_POOL];
			for(int ib=0; ib<wfc_g[ik].nr; ++ib)
			{
				const std::vector<std::complex<double>> wfc_r = cal_wfc_r(wfc_g[ik], ik, ib);

				std::vector<double> wfc_r2(wfc_r.size());
                if (square)
                    for (int ir = 0; ir < wfc_r2.size(); ++ir)
                        wfc_r2[ir] = std::norm(wfc_r[ir]);   // "std::norm(z)" returns |z|^2 
                else
                    for (int ir = 0; ir < wfc_r2.size(); ++ir)
                        wfc_r2[ir] = std::abs(wfc_r[ir]);

				const std::string file_name = outdir + "wfc_realspace_"
					+ ModuleBase::GlobalFunc::TO_STRING(ik_out)
					+ "_" + ModuleBase::GlobalFunc::TO_STRING(ib);
#ifdef __MPI
				mpi_requests.push_back({});
				write_charge_realspace_1(wfc_r2, file_name, mpi_requests.back());
#else
				write_charge_realspace_1(wfc_r2, file_name);
#endif
			}
		}
#ifdef __MPI
		MPI_Waitall( mpi_requests.size(), mpi_requests.data(), MPI_STATUSES_IGNORE );
#endif
		ModuleBase::timer::tick("Write_Wfc_Realspace", "write_wfc_realspace_1");
	}
	// processes output pipeline:
	//
	//           t0  t1  t2  t3  t4  t5  t6  t7
	//          -------------------------------->
	//  rank0    k0  k1  k2  k3  k4  k5
	//             \   \   \   \   \   \
	//  rank1        k0  k1  k2  k3  k4  k5
	//                 \   \   \   \   \   \
	//  rank2            k0  k1  k2  k3  k4  k5



	// Input: wfc_g(ib,ig)
	// Output: wfc_r[ir]
	std::vector<std::complex<double>> cal_wfc_r(const ModuleBase::ComplexMatrix &wfc_g, const int ik, const int ib)
	{
		ModuleBase::timer::tick("Write_Wfc_Realspace", "cal_wfc_r");
		ModuleBase::GlobalFunc::ZEROS(GlobalC::UFFT.porter, GlobalC::pw.nrxx);
		std::vector<std::complex<double>> wfc_r(GlobalC::pw.nrxx);
		for(int ig=0; ig<GlobalC::kv.ngk[ik]; ++ig)
			GlobalC::UFFT.porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik,ig)] ] = wfc_g(ib,ig);
		GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter,1);
		for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
			wfc_r[ir] = GlobalC::UFFT.porter[ir];
		ModuleBase::timer::tick("Write_Wfc_Realspace", "cal_wfc_r");
		return wfc_r;
	}


	
	// Input: chg_r[ir]
#ifdef  __MPI
	void write_charge_realspace_1(const std::vector<double> &chg_r, const std::string &file_name, MPI_Request &mpi_request)
#else
	void write_charge_realspace_1(const std::vector<double> &chg_r, const std::string &file_name)
#endif
	{
		ModuleBase::timer::tick("Write_Wfc_Realspace", "write_charge_realspace_1");
		std::ofstream ofs;

#ifdef  __MPI
		constexpr int mpi_tag=100;
		if(GlobalV::RANK_IN_POOL==0)
		{
#endif
			ofs.open(file_name);

			ofs<<"calculated by ABACUS"<<std::endl;
			ofs<<GlobalC::ucell.lat0_angstrom<<std::endl;
			ofs<<GlobalC::ucell.latvec.e11<<" "<<GlobalC::ucell.latvec.e12<<" "<<GlobalC::ucell.latvec.e13<<std::endl
			<<GlobalC::ucell.latvec.e21<<" "<<GlobalC::ucell.latvec.e22<<" "<<GlobalC::ucell.latvec.e23<<std::endl
			<<GlobalC::ucell.latvec.e31<<" "<<GlobalC::ucell.latvec.e32<<" "<<GlobalC::ucell.latvec.e33<<std::endl;

			for(int it=0; it<GlobalC::ucell.ntype; ++it)
				ofs<<GlobalC::ucell.atoms[it].label<<"\t";
			ofs<<std::endl;
			for(int it=0; it<GlobalC::ucell.ntype; ++it)
				ofs<<GlobalC::ucell.atoms[it].na<<"\t";
			ofs<<std::endl;

			ofs<<"Direct"<<std::endl;
			for(int it=0; it<GlobalC::ucell.ntype; ++it)
				for(int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
					ofs<<GlobalC::ucell.atoms[it].taud[ia].x<<" "<<GlobalC::ucell.atoms[it].taud[ia].y<<" "<<GlobalC::ucell.atoms[it].taud[ia].z<<std::endl;
			ofs<<std::endl;
			
			ofs<<GlobalC::pw.ncx<<" "<<GlobalC::pw.ncy<<" "<<GlobalC::pw.ncz<<std::endl;
#ifdef  __MPI
		}
		else
		{
			char recv_tmp;
			MPI_Recv( &recv_tmp, 1, MPI_CHAR, GlobalV::RANK_IN_POOL-1, mpi_tag, POOL_WORLD, MPI_STATUS_IGNORE);

			ofs.open(file_name, ofstream::app);
		}
#endif
		
		assert(GlobalC::pw.ncx * GlobalC::pw.ncy * GlobalC::pw.nczp == chg_r.size());
		for(int iz=0; iz<GlobalC::pw.nczp; ++iz)
		{
			for(int iy=0; iy<GlobalC::pw.ncy; ++iy)
			{
				for(int ix=0; ix<GlobalC::pw.ncx; ++ix)
				{
					const int ir = (ix*GlobalC::pw.ncy+iy)*GlobalC::pw.nczp+iz;
					ofs<<chg_r[ir]<<" ";
				}
				ofs<<"\n";
			}
		}
		ofs.close();

#ifdef  __MPI
		if(GlobalV::RANK_IN_POOL < GlobalV::NPROC_IN_POOL-1)
		{
			const char send_tmp = 'c';
			MPI_Isend( &send_tmp, 1, MPI_CHAR, GlobalV::RANK_IN_POOL+1, mpi_tag, POOL_WORLD, &mpi_request );
		}
		else
		{
			mpi_request = MPI_REQUEST_NULL;
		}
#endif
		ModuleBase::timer::tick("Write_Wfc_Realspace", "write_charge_realspace_1");
	}
};

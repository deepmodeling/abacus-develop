#include "module_io/write_occ.h"
#include "src_parallel/parallel_reduce.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void ModuleIO::print_occ(const elecstate::ElecState* pelec)
{
	
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "istate.info" ;
		if(GlobalV::MY_RANK==0)
		{
			std::ofstream ofsi( ss.str().c_str() ); // clear istate.info
			ofsi.close();
		}
#ifdef __MPI
		for(int ip=0; ip<GlobalV::KPAR; ip++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if( GlobalV::MY_POOL == ip )
			{
				if( GlobalV::RANK_IN_POOL != 0 || GlobalV::MY_STOGROUP != 0 ) continue;
#endif
				std::ofstream ofsi2( ss.str().c_str(), ios::app );
				if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
				{
					for (int ik = 0;ik < GlobalC::kv.nks;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Energy(ev)"
						<<std::setw(25)<<"Occupation"
#ifdef __MPI
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
#else
						<<std::setw(25)<<"Kpoint = "<<ik+1
#endif
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;
						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1<<std::setw(25)<<pelec->ekb(ik,ib)* ModuleBase::Ry_to_eV<<std::setw(25)<<pelec->wg(ik,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;
					}
				}
				else
				{
					for (int ik = 0;ik < GlobalC::kv.nks/2;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Spin up Energy(ev)"
						<<std::setw(25)<<"Occupation"
						<<std::setw(25)<<"Spin down Energy(ev)"
						<<std::setw(25)<<"Occupation"
#ifdef __MPI
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
#else
						<<std::setw(25)<<"Kpoint = "<<ik+1
#endif
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;

						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1
							<<std::setw(25)<<pelec->ekb(ik, ib)* ModuleBase::Ry_to_eV
							<<std::setw(25)<<pelec->wg(ik,ib)
							<<std::setw(25)<<pelec->ekb((ik+GlobalC::kv.nks/2), ib)* ModuleBase::Ry_to_eV
							<<std::setw(25)<<pelec->wg(ik+GlobalC::kv.nks/2,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;

					}
				}

				ofsi2.close();
#ifdef __MPI
			}
		}
#endif
}

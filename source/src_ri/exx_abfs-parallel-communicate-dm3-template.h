#include "exx_abfs-parallel-communicate-dm3.h"
#include "../src_pw/global.h"

#ifdef __MPI
const ModuleBase::matrix &Exx_Abfs::Parallel::Communicate::DM3::D_phase(
	const ModuleBase::matrix &DK, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	assert(box2 == Abfs::Vector3_Order<int>(0,0,0));
	return DK;
}
ModuleBase::matrix Exx_Abfs::Parallel::Communicate::DM3::D_phase(
	const ModuleBase::ComplexMatrix &DK, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	return (DK * exp( -ModuleBase::TWO_PI*ModuleBase::IMAG_UNIT * (GlobalC::kv.kvec_c[ik] * (box2*GlobalC::ucell.latvec)) )).real();
}

template<typename Tmatrix> std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
Exx_Abfs::Parallel::Communicate::DM3::K_to_R(const std::vector<Tmatrix> &DK_2D, const double threshold_D, const Parallel_Orbitals &pv) const
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM3::K_to_R");

	/*{
		static int istep=0;
		std::ofstream ofs("DK_2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<DK_2D<<std::endl;
	}*/
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DR_a2D(GlobalV::NSPIN);
	
	const std::map<int,int> nspin_2D = {{1,1}, {2,2}, {4,1}};
	const std::map<int,double> SPIN_multiple = {{1,0.5}, {2,1}, {4,1}};							// ???
	const Abfs::Vector3_Order<int> Born_von_Karman_period = ModuleBase::Vector3<int>{GlobalC::kv.nmp[0],GlobalC::kv.nmp[1],GlobalC::kv.nmp[2]};
	const std::vector<Abfs::Vector3_Order<int>> supercell_boxes = Abfs::get_Born_von_Karmen_boxes(Born_von_Karman_period);
	for(const Abfs::Vector3_Order<int> &box2 : supercell_boxes)
	{
		std::vector<ModuleBase::matrix> DR_2D( nspin_2D.at(GlobalV::NSPIN),
			{DK_2D[0].nr, DK_2D[0].nc} );
		for(int ik=0; ik!=DK_2D.size(); ++ik)
			DR_2D[GlobalC::kv.isk[ik]] += D_phase( DK_2D[ik], ik, box2);
		
		// C++: 0 1
		//      2 3
		std::vector<std::map<size_t,std::map<size_t,ModuleBase::matrix>>> DR_a2D_box2(GlobalV::NSPIN);
		for(int is_2D=0; is_2D!=nspin_2D.at(GlobalV::NSPIN); ++is_2D)
		{
			for(int iwt1_local=0; iwt1_local!=DR_2D[is_2D].nr; ++iwt1_local)
			{
				const int iwt1 = pv.MatrixInfo.col_set[iwt1_local];
				const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
				const int iw1 = GlobalC::ucell.iwt2iw[iwt1];
				for(int iwt2_local=0; iwt2_local!=DR_2D[is_2D].nc; ++iwt2_local)
				{
					const int iwt2 = pv.MatrixInfo.row_set[iwt2_local];
					const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
					const int iw2 = GlobalC::ucell.iwt2iw[iwt2];
					
					int is_R=is_2D, iw1_R=iw1, iw2_R=iw2;
					if(GlobalV::NSPIN==4)
					{
						iw1_R = iw1/2;
						iw2_R = iw2/2;
						is_R = iw1%2*2 + iw2%2;
					}
					ModuleBase::matrix &DR_a2D_tmp = DR_a2D_box2[is_R][iat1][iat2];
					if(!DR_a2D_tmp.c)
						DR_a2D_tmp.create(GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
					DR_a2D_tmp(iw1_R,iw2_R) = DR_2D[is_2D](iwt1_local,iwt2_local) * SPIN_multiple.at(GlobalV::NSPIN);
				}
			}
		}
		
		for(auto &DR_a2D_box2_is : DR_a2D_box2)
			Abfs::delete_threshold_ptrs(DR_a2D_box2_is,threshold_D);

		for(int is=0; is!=GlobalV::NSPIN; ++is)
		{
			for(auto & DR_a2D_box2_A : DR_a2D_box2[is])
			{
				const size_t iat1 = DR_a2D_box2_A.first;
				for(auto & DR_a2D_box2_B : DR_a2D_box2_A.second)
				{
					const size_t iat2 = DR_a2D_box2_B.first;
					DR_a2D[is][iat1][iat2][box2] = std::move(DR_a2D_box2_B.second);
				}
			}
		}
	}
	
	/*{
		static int istep=0;
		std::ofstream ofs("DR_a2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<DR_a2D<<std::endl;
	}*/
	return DR_a2D;
}
#endif

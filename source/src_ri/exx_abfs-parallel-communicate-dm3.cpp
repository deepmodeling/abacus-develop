#include "exx_abfs-parallel-communicate-dm3.h"
#include "exx_abfs-parallel-communicate-dm3-template.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"
#include "abfs-template.h"

#include "../src_external/src_test/test_function.h"
#include "../src_external/src_test/src_global/complexmatrix-test.h"

/*
template<> std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
Exx_Abfs::Parallel::Communicate::DM3::K_to_R(const std::vector<ModuleBase::matrix> &DK_2D, const double threshold_D) const
{

	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM3::K_to_R");
	assert(DK_2D.size()==GlobalV::NSPIN);
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DR_a2D(GlobalV::NSPIN);
	for(int is=0; is!=GlobalV::NSPIN; ++is)
	{
		for(int iwt1_local=0; iwt1_local!=DK_2D[is].nr; ++iwt1_local)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.col_set[iwt1_local];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
			const int iw1 = GlobalC::ucell.iwt2iw[iwt1];
			for(int iwt2_local=0; iwt2_local!=DK_2D[is].nc; ++iwt2_local)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.row_set[iwt2_local];
				const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
				const int iw2 = GlobalC::ucell.iwt2iw[iwt2];
				const double dm = DK_2D[is](iwt1_local,iwt2_local);
				if(abs(dm) > threshold_D)
				{
					ModuleBase::matrix &DR_a2D_box2 = DR_a2D[is][iat1][iat2][{0,0,0}];
					if(!DR_a2D_box2.c)
						DR_a2D_box2.create(GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
					DR_a2D_box2(iw1,iw2) = DK_2D[is](iwt1_local,iwt2_local) * SPIN_multiple;
				}
			}
		}
	}
	return DR_a2D;
}
*/

/*
template<> std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>>
Exx_Abfs::Parallel::Communicate::DM3::K_to_R(const std::vector<ModuleBase::ComplexMatrix> &DK_2D, const double threshold_D) const
{
	{
		static int istep=0;
		std::ofstream ofs("DK_2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<DK_2D<<std::endl;
	}
	
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM3::K_to_R");
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	const Abfs::Vector3_Order<int> Born_von_Karman_period = ModuleBase::Vector3<int>{GlobalC::kv.nmp[0],GlobalC::kv.nmp[1],GlobalC::kv.nmp[2]};
	const std::vector<Abfs::Vector3_Order<int>> supercell_boxes = Abfs::get_Born_von_Karmen_boxes(Born_von_Karman_period);
	
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DR_a2D(GlobalV::NSPIN);
	for(int ik=0; ik!=DK_2D.size(); ++ik)
	{
		for(int iwt1_local=0; iwt1_local!=DK_2D[ik].nr; ++iwt1_local)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.col_set[iwt1_local];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
			const int iw1 = GlobalC::ucell.iwt2iw[iwt1];
			for(int iwt2_local=0; iwt2_local!=DK_2D[ik].nc; ++iwt2_local)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.row_set[iwt2_local];
				const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
				const int iw2 = GlobalC::ucell.iwt2iw[iwt2];
				for(const Abfs::Vector3_Order<int> &box2 : supercell_boxes)
				{
					ModuleBase::matrix &DR_a2D_box2 = DR_a2D[GlobalC::kv.isk[ik]][iat1][iat2][box2];
					if(!DR_a2D_box2.c)
						DR_a2D_box2.create(GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
					DR_a2D_box2(iw1,iw2) += real( DK_2D[ik](iwt1_local,iwt2_local) * exp(-ModuleBase::TWO_PI*IMAG_UNIT*(GlobalC::kv.kvec_c[ik]*(box2*GlobalC::ucell.latvec))) ) * SPIN_multiple;
				}
			}
		}
	}
	for(auto &DR_a2D_is : DR_a2D)
		Abfs::delete_threshold_ptrs(DR_a2D_is,threshold_D);
	
	{
		static int istep=0;
		std::ofstream ofs("DR_a2D_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<DR_a2D<<std::endl;
	}	
	return DR_a2D;
}
*/

#ifdef __MPI
void Exx_Abfs::Parallel::Communicate::DM3::cal_DM(const double threshold_D,
    Local_Orbital_Charge &loc)
{
	ModuleBase::TITLE("Exx_Abfs::Parallel::Communicate::DM3::cal_DM");
	std::vector<std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,ModuleBase::matrix>>>> DR_a2D = GlobalV::GAMMA_ONLY_LOCAL
		? K_to_R(loc.dm_gamma, threshold_D, *loc.ParaV)
		: K_to_R(loc.dm_k, threshold_D, *loc.ParaV);
	DMr = allreduce.a2D_to_exx(DR_a2D);

	/*{
		static int istep=0;
		std::ofstream ofs("DMr_"+ModuleBase::GlobalFunc::TO_STRING(istep++)+"_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK));
		ofs<<DMr<<std::endl;
	}*/	
}
#endif

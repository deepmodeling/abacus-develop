#include "LCAO_deepks_test.h"

namespace Test_Deepks
{
	Grid_Driver GridD(GlobalV::test_deconstructor,GlobalV::test_grid_driver,GlobalV::test_grid);
}

test_deepks::test_deepks()
{}

test_deepks::~test_deepks()
{}


void test_deepks::check_dstable(void)
{
	OGT.talpha.print_Table_DSR(ORB);
	this->compare_with_ref("S_I_mu_alpha.dat","S_I_mu_alpha_ref.dat");
}

void test_deepks::check_psialpha(void)
{
	std::vector<int> na;
	na.resize(ucell.ntype);
	for(int it=0;it<ucell.ntype;it++)
	{
		na[it] = ucell.atoms[it].na;
	}
	ld.init(ORB,
		ucell.nat,
		ucell.ntype,
		na);

	ld.build_psialpha(GlobalV::CAL_FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);

	ld.check_psialpha(GlobalV::CAL_FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);
	this->compare_with_ref("psialpha.dat","psialpha_ref.dat");
	this->compare_with_ref("dpsialpha_x.dat","dpsialpha_x_ref.dat");
	this->compare_with_ref("dpsialpha_y.dat","dpsialpha_y_ref.dat");
	this->compare_with_ref("dpsialpha_z.dat","dpsialpha_z_ref.dat");	
}

void test_deepks::read_dm(void)
{
    std::ifstream ifs("dm");
    dm.resize(1);
    dm[0].create(GlobalV::NLOCAL, GlobalV::NLOCAL);

	for (int mu=0;mu<GlobalV::NLOCAL;mu++)
	{
		for (int nu=0;nu<GlobalV::NLOCAL;nu++)
		{
			double c;
			ifs >> c;
			dm[0](mu,nu)=c;
		}
	}
}

void test_deepks::read_dm_k(const int nks)
{
	dm_k.resize(nks);
	std::stringstream ss;
	for(int ik=0;ik<nks;ik++)
	{
        ss.str("");
        ss<<"dm_"<<ik;
        std::ifstream ifs(ss.str().c_str());
		dm_k[ik].create(GlobalV::NLOCAL,GlobalV::NLOCAL);

		for (int mu=0;mu<GlobalV::NLOCAL;mu++)
		{
			for (int nu=0;nu<GlobalV::NLOCAL;nu++)
			{
				std::complex<double> c;
				ifs >> c;
				dm_k[ik](mu,nu)=c;
			}
		}
	}
}

void test_deepks::check_pdm(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->read_dm();
		this->ld.cal_projected_DM(dm[0],
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col);
	}
	else
	{
		this->read_dm_k(kv.nkstot);
		this->ld.cal_projected_DM_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d);		
	}
	this->ld.check_projected_dm();
	this->compare_with_ref("pdm.dat","pdm_ref.dat");
}

void test_deepks::check_gdmx(void)
{
	this->ld.init_gdmx(ucell.nat);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_gdmx(dm[0],
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			0);
	}
	else
	{
		this->ld.cal_gdmx_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d,
			0);			
	}
	this->ld.check_gdmx(ucell.nat);

	for(int ia=0;ia<ucell.nat;ia++)
	{
		std::stringstream ss;
		std::stringstream ss1;
		ss.str("");
        ss<<"gdmx_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gdmx_"<<ia<<"_ref.dat";
		
		this->compare_with_ref(ss.str(),ss1.str());

        ss.str("");
        ss<<"gdmy_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gdmy_"<<ia<<"_ref.dat";
		this->compare_with_ref(ss.str(),ss1.str());

        ss.str("");
        ss<<"gdmz_"<<ia<<".dat";
        ss1.str("");
        ss1<<"gdmz_"<<ia<<"_ref.dat";
		this->compare_with_ref(ss.str(),ss1.str());
	}	
}

void test_deepks::check_descriptor(void)
{
	this->ld.cal_descriptor();
	this->ld.check_descriptor(ucell);
	this->compare_with_ref("descriptor.dat","descriptor_ref.dat");
}

void test_deepks::check_gvx(void)
{
	this->ld.cal_gvx(ucell.nat);
	this->ld.check_gvx(ucell.nat);

	for(int ia=0;ia<ucell.nat;ia++)
	{
		std::stringstream ss;
		std::stringstream ss1;
		ss.str("");
        ss<<"gvx_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gvx_"<<ia<<"_ref.dat";
		this->compare_with_ref(ss.str(),ss1.str());

        ss.str("");
        ss<<"gvy_"<<ia<<".dat";
		ss1.str("");
        ss1<<"gvy_"<<ia<<"_ref.dat";
		this->compare_with_ref(ss.str(),ss1.str());

        ss.str("");
        ss<<"gvz_"<<ia<<".dat";
        ss1.str("");
        ss1<<"gvz_"<<ia<<"_ref.dat";
		this->compare_with_ref(ss.str(),ss1.str());
	}
}

void test_deepks::check_edelta(void)
{
	this->ld.load_model("model.ptg");
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc);
	}
	else
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc,kv.nkstot);
	}
	this->ld.cal_gedm(ucell.nat);

	std::ofstream ofs("E_delta.dat");
	ofs << std::setprecision(10) << this->ld.E_delta << std::endl;
	ofs.close();
	this->compare_with_ref("E_delta.dat","E_delta_ref.dat");
	
	this->ld.check_gedm();
	this->compare_with_ref("gedm.dat","gedm_ref.dat");
}

void test_deepks::check_e_deltabands(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_e_delta_band(dm,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			ParaO.nrow);
	}
	else
	{
		this->folding_nnr(kv);
		this->ld.cal_e_delta_band_k(dm_k,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			ParaO.nrow,
			ParaO.ncol);
	}

	std::ofstream ofs("E_delta_bands.dat");
	ofs << std::setprecision(10) << this->ld.e_delta_band << std::endl;
	ofs.close();
	this->compare_with_ref("E_delta_bands.dat","E_delta_bands_ref.dat");
}

void test_deepks::check_v_delta()
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.add_v_delta(ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
            ParaO.trace_loc_col,
            ParaO.nrow,
            ParaO.ncol);
		this->ld.check_v_delta(ParaO.nrow,ParaO.ncol);
		this->compare_with_ref("H_V_delta.dat","H_V_delta_ref.dat");
	}
	else
	{
		this->cal_nnr();
		this->ld.allocate_V_deltaR(nnr);
		this->ld.add_v_delta_k(ucell,
        	ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			nnr);
		this->ld.check_v_delta_k(nnr);
		this->compare_with_ref("H_V_deltaR.dat","H_V_deltaR_ref.dat");
	}
}

void test_deepks::check_f_delta()
{
	ModuleBase::matrix svnl_dalpha;
	svnl_dalpha.create(3,3);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		ld.cal_f_delta_gamma(dm[0],
            ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
            1, svnl_dalpha);
	}
	else
	{
		ld.cal_f_delta_k(dm_k,
			ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d,
			1,svnl_dalpha);
	}
	ld.check_f_delta(ucell.nat, svnl_dalpha);

	this->compare_with_ref("F_delta.dat","F_delta_ref.dat");
}

void test_deepks::compare_with_ref(
	const std::string f1,
	const std::string f2)
{
	this->total_check+=1;
	std::ifstream file1(f1.c_str());
	std::ifstream file2(f2.c_str());
	double test_thr=1e-8;

	std::string word1;
	std::string word2;
	while(file1 >> word1)
	{
		file2 >> word2;
		if((word1[0]-'0'>=0 && word1[0]-'0'<10)||word1[0]=='-')
		{
			double num1 = std::stof(word1);
			double num2 = std::stof(word2);
			if(std::abs(num1-num2)>test_thr)
			{
				this->failed_check+=1;
				std::cout << "\e[1;31m [  FAILED  ] \e[0m" << f1.c_str() << " inconsistent!" << std::endl;
				return;
			}			
		}
		else
		{
			if(word1!=word2)
			{
				this->failed_check+=1;
				return;
			}
		}
	}
	return;
}
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include<streambuf>

/************************************************
 *  unit test of pseudo
 ***********************************************/

/**
 * - Tested Functions:
 *   - pseudo
 *   - ~pseudo
 *   - set_pseudo_h
 *   - set_pseudo_atom
 *   - set_pseudo_vl
 *   - set_pseudo_nc
 *   - print_pseudo_h
 *   - print_pseudo_atom
 *   - print_pseudo_vl
 *   - print_pseudo_nc
 */

#define private public
#include "module_cell/read_pp.h"
#include "module_cell/pseudo.h"

class NCPPTest : public testing::Test
{
protected:
	std::unique_ptr<Pseudopot_upf> upf{new Pseudopot_upf};
	std::unique_ptr<pseudo> ncpp{new pseudo};
};

TEST_F(NCPPTest, SetPseudoH)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs);
	//set_pseudo_h
	ncpp->set_pseudo_h(*upf);
	EXPECT_EQ(ncpp->nv,upf->nv);
	EXPECT_EQ(ncpp->psd,upf->psd);
	EXPECT_EQ(ncpp->pp_type,upf->pp_type);
	EXPECT_EQ(ncpp->tvanp,upf->tvanp);
	EXPECT_EQ(ncpp->nlcc,upf->nlcc);
	EXPECT_EQ(ncpp->xc_func,upf->xc_func);
	EXPECT_EQ(ncpp->zv,upf->zp);
	EXPECT_EQ(ncpp->etotps,upf->etotps);
	EXPECT_EQ(ncpp->ecutwfc,upf->ecutwfc);
	EXPECT_EQ(ncpp->ecutrho,upf->ecutrho);
	EXPECT_EQ(ncpp->lmax,upf->lmax);
	EXPECT_EQ(ncpp->mesh,upf->mesh);
	EXPECT_EQ(ncpp->nchi,upf->nwfc);
	EXPECT_EQ(ncpp->nbeta,upf->nbeta);
	for (int i=0;i<ncpp->nchi;i++)
	{
		EXPECT_EQ(ncpp->els[i],upf->els[i]);
		EXPECT_EQ(ncpp->lchi[i],upf->lchi[i]);
		EXPECT_EQ(ncpp->oc[i],upf->oc[i]);
	}
	EXPECT_EQ(ncpp->has_so,upf->has_so);
	if(ncpp->has_so)
	{
		for (int i=0;i<ncpp->nchi;i++)
		{
			EXPECT_EQ(ncpp->nn[i],upf->nn[i]);
			EXPECT_EQ(ncpp->jchi[i],upf->jchi[i]);
		}
		for (int i=0;i<ncpp->nbeta;i++)
		{
			EXPECT_EQ(ncpp->jjj[i],upf->jjj[i]);
		}
	}
	else
	{
		for (int i=0;i<ncpp->nchi;i++)
		{
			EXPECT_EQ(ncpp->nn[i],0);
			EXPECT_EQ(ncpp->jchi[i],0);
		}
		for (int i=0;i<ncpp->nbeta;i++)
		{
			EXPECT_EQ(ncpp->jjj[i],0);
		}
	}
	ifs.close();
}

TEST_F(NCPPTest, SetPseudoAtom)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs);
	//set_pseudo_atom
	ncpp->set_pseudo_atom(*upf);
	EXPECT_EQ(ncpp->rcut,GlobalV::PSEUDORCUT);
	for(int i=0;i<ncpp->nchi;i++)
	{
		for(int j=0;j<ncpp->mesh;j++)
		{
			EXPECT_EQ(ncpp->chi(i,j),upf->chi(i,j));
		}
	}
	for(int i=0;i<ncpp->mesh;i++)
	{
		EXPECT_EQ(ncpp->r[i],upf->r[i]);
		EXPECT_EQ(ncpp->rab[i],upf->rab[i]);
		EXPECT_EQ(ncpp->rho_at[i],upf->rho_at[i]);
	}
	if(ncpp->nlcc)
	{
		for(int i=0;i<ncpp->mesh;i++)
		{
			EXPECT_EQ(ncpp->rho_atc[i],upf->rho_atc[i]);
		}
	}
	else
	{
		for(int i=0;i<ncpp->mesh;i++)
		{
			EXPECT_EQ(ncpp->rho_atc[i],0.0);
		}
	}
	EXPECT_EQ(ncpp->msh,ncpp->mesh);
	ifs.close();
}

TEST_F(NCPPTest, SetPseudoNC)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs);
	//set_pseudo_nc
	ncpp->set_pseudo_nc(*upf);
	for(int i=0;i<ncpp->nbeta;i++)
	{
		EXPECT_EQ(ncpp->lll[i],upf->lll[i]);
	}
	EXPECT_EQ(ncpp->nh,14);
	EXPECT_EQ(ncpp->kkbeta,132);
	ifs.close();
}

TEST_F(NCPPTest, PrintNC)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs);
	ncpp->set_pseudo_nc(*upf);
	ifs.close();
	//print
	std::ofstream ofs;
	ofs.open("./tmp_log");
	ncpp->print_pseudo_nc(ofs);
	ofs.close();
	ifs.open("./tmp_log");
	std::string str((std::istreambuf_iterator<char>(ifs)),std::istreambuf_iterator<char>());
	EXPECT_THAT(str,testing::HasSubstr("psd  C"));
	EXPECT_THAT(str,testing::HasSubstr("pp_type  NC"));
	EXPECT_THAT(str,testing::HasSubstr("dft  PBE"));
	EXPECT_THAT(str,testing::HasSubstr("zv       4"));
	EXPECT_THAT(str,testing::HasSubstr("nchi 3"));
	EXPECT_THAT(str,testing::HasSubstr("nbeta    6"));
	EXPECT_THAT(str,testing::HasSubstr("dion :  nr=6 nc=6"));
	EXPECT_THAT(str,testing::HasSubstr("msh\t1247"));
	ifs.close();
	remove("./tmp_log");
}
#undef private

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
 *   - set_pseudo
 *   - print_pseudo_h
 *   - print_pseudo_atom
 *   - print_pseudo_vl
 *   - print_pseudo
 */

#define private public
#include "module_cell/read_pp.h"
#include "module_cell/atom_pseudo.h"

class NCPPTest : public testing::Test
{
protected:
	std::unique_ptr<Pseudopot_upf> upf{new Pseudopot_upf};
	std::unique_ptr<Atom_pseudo> ncpp{new Atom_pseudo};
};

TEST_F(NCPPTest, SetPseudoH)
{
	std::ifstream ifs;
	//set
	ifs.open("./support/C.upf");
	GlobalV::PSEUDORCUT = 15.0;
	upf->read_pseudo_upf201(ifs, *ncpp);
	//set_pseudo_h
	ncpp->set_pseudo_h();

	if(!ncpp->has_so)
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
	upf->read_pseudo_upf201(ifs, *ncpp);
	//set_pseudo_atom
	ncpp->set_pseudo_h();
	ncpp->set_pseudo_atom();
	EXPECT_EQ(ncpp->rcut,GlobalV::PSEUDORCUT);

	if(!ncpp->nlcc)
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
	// set pseudo nbeta = 0
	upf->read_pseudo_upf201(ifs, *ncpp);
	ncpp->nbeta = 0;
	ncpp->set_pseudo();
	EXPECT_EQ(ncpp->nh,0);
    // set_pseudo nbeta > 0
	upf->read_pseudo_upf201(ifs, *ncpp);
    ncpp->set_pseudo();
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
	upf->read_pseudo_upf201(ifs, *ncpp);
    ncpp->set_pseudo();
    ifs.close();
	//print
	std::ofstream ofs;
	ofs.open("./tmp_log");
	ncpp->print_pseudo(ofs);
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

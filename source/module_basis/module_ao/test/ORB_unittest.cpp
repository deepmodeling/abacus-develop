#include "ORB_unittest.h"

void test_orb::SetUp()
{
	//test constructor
	/*Center2_Orb::Orb11 testcto = Center2_Orb::Orb11(
		ORB.Phi[0].PhiLN(0, 0),
		ORB.Phi[0].PhiLN(0, 0),
		OGT.MOT, Center2_MGT);*/
	// 1. setup orbitals
	this->ofs_running.open("log.txt");
	this->count_ntype();
    this->set_files();
    this->set_ekcut();


	//2. setup 2-center-integral tables by basic methods 
	// not including center2orb, it will be set up when needed
	// in some test cases.
	this->set_orbs();
	//this->set_center2orbs();

}

void test_orb::TearDown()
{
	int* nproj = new int[ORB.get_ntype()];
	for (int i = 0;i < ORB.get_ntype();++i)
		nproj[i] = 0;
	ooo.clear_after_ions(OGT, ORB, 0, nproj);
	delete[] nproj;
    delete[] orbital_fn;
	return;
}

void test_orb::set_ekcut()
{
	std::cout << "set lcao_ecut from LCAO files" << std::endl;
	//set as max of ekcut from every element

	lcao_ecut=0.0;
	std::ifstream in_ao;

	for(int it=0;it<ntype_read;it++)
	{
		double ek_current;

		in_ao.open((this->case_dir+ORB.orbital_file[it].c_str()));
		if(!in_ao)
		{
			std::cout << "error : cannot find LCAO file : " << ORB.orbital_file[it] << std::endl;
		}
		ORB.orbital_file[it] = this->case_dir + ORB.orbital_file[it].c_str();
		std::string word;
		while (in_ao.good())
		{
			in_ao >> word;
			if(word == "Cutoff(Ry)") break;
		}
		in_ao >> ek_current;
		lcao_ecut = std::max(lcao_ecut,ek_current);

		in_ao.close();
	}

	ORB.ecutwfc=lcao_ecut;
	std::cout << "lcao_ecut : " << lcao_ecut << std::endl;
	return;
}

void test_orb::set_orbs()
{

	ooo.read_orb_first(
		ofs_running,
		ORB,
		ntype_read,
        "./",
        orbital_fn,
        descriptor_file,
		lmax,
		lcao_ecut,
		lcao_dk,
		lcao_dr,
		lcao_rmax,
		0,
		0,
		1,//force
		0);//myrank
	
	int* nproj = new int[ORB.get_ntype()];
	for (int i = 0;i < ORB.get_ntype();++i)
		nproj[i] = 0;
	const Numerical_Nonlocal beta_[ORB.get_ntype()];

	ooo.set_orb_tables(
		ofs_running,
		OGT,
		ORB,
		lat0,
		0,  //no out_descriptor
		lmax,
		0,  //no nproj
		nproj,
		beta_);
	
	delete[] nproj;
	return;
}

void test_orb::set_files()
{
	std::cout << "read names of atomic basis set files" << std::endl;
	std::ifstream ifs((this->case_dir + "STRU"),std::ios::in);

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"NUMERICAL_ORBITAL");

    orbital_fn = new std::string [ntype_read];

	for(int it=0;it<ntype_read;it++)
	{
		ifs >> orbital_fn[it];
		ORB.orbital_file.push_back(orbital_fn[it]);

		std::cout << "Numerical orbital file : " << orbital_fn[it] << std::endl;
	}
	
	return;
}

void test_orb::count_ntype()
{
	std::cout << "count number of atom types" << std::endl;
		std::cout << this->case_dir +"STRU" << std::endl;
	std::ifstream ifs( (this->case_dir+ "STRU"), std::ios::in);

	if (!ifs)
	{
		std::cout << "ERROR : file STRU does not exist" <<std::endl;
		exit(1);
	}

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"ATOMIC_SPECIES");	
	
	ntype_read = 0;

	std::string x;
	ifs.rdstate();
	while ( ifs.good() )
	{
//read a line
		std::getline(ifs,x);

//trim white space
		const char* typeOfWhitespaces = " \t\n\r\f\v";
		x.erase(x.find_last_not_of(typeOfWhitespaces) + 1);
		x.erase(0,x.find_first_not_of(typeOfWhitespaces));

		if(x=="LATTICE_CONSTANT" || x=="NUMERICAL_ORBITAL" || x=="LATTICE_VECTORS" || x=="ATOMIC_POSITIONS") break;

		std::string tmpid=x.substr(0,1);
		if(!x.empty() && tmpid!="#") ntype_read++;
	}
	std::cout <<"ntype="<< ntype_read << std::endl;
	ifs.close();

	return;
}


void test_orb::set_center2orbs()
{
	//1. setup Gaunt coeffs
    Center2_MGT.init_Gaunt_CH( lmax );
	Center2_MGT.init_Gaunt(lmax);
	//2. setup tables
	
	for (int TA = 0; TA < ORB.get_ntype(); TA++)
    {
        for (int TB = 0;  TB < ORB.get_ntype(); TB++)
        {
            for (int LA=0; LA <= ORB.Phi[TA].getLmax() ; LA++)
            {
                for (int NA = 0; NA < ORB.Phi[TA].getNchi(LA); ++NA)
                {
                    for (int LB = 0; LB <= ORB.Phi[TB].getLmax(); ++LB)
                    {
                        for (int NB = 0; NB < ORB.Phi[TB].getNchi(LB); ++NB)
						{
							this->set_single_c2o<Center2_Orb::Orb11>(TA, TB, LA, NA, LB, NB);
							// test_center2_orb11[TA][TB][LA][NA][LB].insert(
							// 	make_pair(NB, MockCenter2Orb11(ORB.Phi[TA].PhiLN(LA, NA),
							// 		ORB.Phi[TB].PhiLN(LB, NB), OGT.MOT, Center2_MGT)));
						}
                    }
                }
            }
        }
	}
	
	for (auto& co1 : this->test_center2_orb11)
        for( auto &co2 : co1.second )
            for( auto &co3 : co2.second )
                for( auto &co4 : co3.second )
                    for( auto &co5 : co4.second )
                        for( auto &co6 : co5.second )
                            co6.second->init_radial_table();
}
template <class c2o>
void test_orb::set_single_c2o(int TA, int TB, int LA, int NA, int LB, int NB)
{
	this->test_center2_orb11[TA][TB][LA][NA][LB].insert(
	std::make_pair(NB, std::make_unique<c2o>(ORB.Phi[TA].PhiLN(LA, NA),
		ORB.Phi[TB].PhiLN(LB, NB), OGT.MOT, Center2_MGT)));
}
double test_orb::randr(double Rmax)
{
    return double(rand()) / double(RAND_MAX) * Rmax;
}

/*
void test_orb::test() {
	ModuleBase::Vector3<double> R1(0, 0, 0);
	ModuleBase::Vector3<double> R2(randr(50), randr(50), randr(50));
	std::cout << "random R2=(" << R2.x << "," << R2.y << "," << R2.z << ")" << std::endl;
	ModuleBase::Vector3<double> dR = ModuleBase::Vector3<double>(0.001, 0.001, 0.001);
	std::cout << this->test_center2_orb11[0][0][0][0][0][0].cal_overlap(R1, R2, 0, 0);
}*/
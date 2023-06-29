#include "LCAO_deepks_test.h"

void test_deepks::preparation()
{
    this->count_ntype();
	this->set_parameters();

    this->setup_cell();

	this->setup_kpt();

	this->set_ekcut();
	this->set_orbs(ucell.lat0);
	this->prep_neighbour();

	this->ParaO.set_trace();
}

void test_deepks::set_parameters()
{
	GlobalV::BASIS_TYPE = "lcao";
	// GlobalV::global_pseudo_type= "auto";
	GlobalV::PSEUDORCUT = 15.0;
	GlobalV::global_out_dir="./";
	GlobalV::ofs_warning.open("warning.log");
	GlobalV::ofs_running.open("running.log");
	GlobalV::deepks_setorb=1;
	GlobalV::CAL_FORCE=1;
	
	std::ifstream ifs("INPUT");
	char word[80];
	ifs >> word;
	ifs >> GlobalV::GAMMA_ONLY_LOCAL;
	ifs.close();

	ucell.latName = "none";
	ucell.ntype = ntype;
	return;
}

void test_deepks::count_ntype()
{
	GlobalV::ofs_running << "count number of atom types" << std::endl;
	std::ifstream ifs("STRU",std::ios::in);

	if (!ifs)
	{
		GlobalV::ofs_running << "ERROR : file STRU does not exist" <<std::endl;
		exit(1);
	}

	ModuleBase::GlobalFunc::SCAN_BEGIN(ifs,"ATOMIC_SPECIES");	
	
	ntype = 0;

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

		if(x=="LATTICE_CONSTANT" || x=="NUMERICAL_ORBITAL" 
            || x=="LATTICE_VECTORS" || x=="ATOMIC_POSITIONS"
            || x=="NUMERICAL_DESCRIPTOR") break;

		std::string tmpid=x.substr(0,1);
		if(!x.empty() && tmpid!="#") ntype++;
	}

	GlobalV::ofs_running << "ntype : "<< ntype << std::endl;
	ifs.close();

	return;
}

void test_deepks::set_ekcut()
{
	GlobalV::ofs_running << "set lcao_ecut from LCAO files" << std::endl;
	//set as max of ekcut from every element

	lcao_ecut=0.0;
	std::ifstream in_ao;

	for(int it=0;it<ntype;it++)
	{
		double ek_current;

		in_ao.open(ucell.orbital_fn[it].c_str());
		if(!in_ao)
		{
			GlobalV::ofs_running << "error : cannot find LCAO file : " << ucell.orbital_fn[it] << std::endl;
		}

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
	GlobalV::ofs_running << "lcao_ecut : " << lcao_ecut << std::endl;
	
	return;
}

void test_deepks::setup_cell()
{
	ucell.setup_cell("STRU", GlobalV::ofs_running);
    ucell.read_pseudo(GlobalV::ofs_running);

	return;
}

void test_deepks::prep_neighbour()
{
    double search_radius = atom_arrange::set_sr_NL(
        GlobalV::ofs_running,
        GlobalV::OUT_LEVEL,
        ORB.get_rcutmax_Phi(),
		ucell.infoNL.get_rcutmax_Beta(),
        GlobalV::GAMMA_ONLY_LOCAL);

    atom_arrange::search(
        GlobalV::SEARCH_PBC,
        GlobalV::ofs_running,
        Test_Deepks::GridD,
        ucell,
        search_radius,
        GlobalV::test_atom_input);
}

void test_deepks::set_orbs(const double &lat0_in)
{
	for(int it=0;it<ntype;it++)
	{
		ooo.read_orb_first(
			GlobalV::ofs_running,
			ORB,
			ucell.ntype,
            GlobalV::global_orbital_dir,
            ucell.orbital_fn,
            ucell.descriptor_file,
			ucell.lmax,
			lcao_ecut,
			lcao_dk,
			lcao_dr,
			lcao_rmax,
			GlobalV::deepks_setorb,
			out_mat_r,
			GlobalV::CAL_FORCE,
			my_rank);
		
		ucell.infoNL.setupNonlocal(
			ucell.ntype,
			ucell.atoms,
			GlobalV::ofs_running,
			ORB);

		ooo.set_orb_tables(
			GlobalV::ofs_running,
			OGT,
			ORB,
			ucell.lat0,
			GlobalV::deepks_setorb,
			lmax,
			ucell.infoNL.nprojmax,
			ucell.infoNL.nproj,
			ucell.infoNL.Beta);
        GlobalV::ofs_running << "read and set from orbital_file : " << ORB.orbital_file[it] << std::endl;
        GlobalV::ofs_running << "ucell.ntype, ucell.lmax : " << ucell.ntype << " " << ucell.lmax << std::endl;
	}
	return;
}

void test_deepks::setup_kpt()
{
	this->kv.set(
		"KPT",
		GlobalV::NSPIN,
		ucell.G,
		ucell.latvec,
		GlobalV::GAMMA_ONLY_LOCAL,
		GlobalV::ofs_running,
		GlobalV::ofs_warning);
}
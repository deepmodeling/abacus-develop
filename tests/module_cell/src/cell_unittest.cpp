#include "cell_unittest.h"
#ifdef __LCAO
namespace Test_Cell
{
	Grid_Driver GridD(GlobalV::test_deconstructor, GlobalV::test_grid_driver,GlobalV::test_grid);
}

test_cell_orb::test_cell_orb()
{}
test_cell_orb::~test_cell_orb()
{}
#else
test_cell::test_cell()
{}
test_cell::~test_cell()
{}
#endif

#ifdef __LCAO
void test_cell_orb::set_parameters()
#else
void test_cell::set_parameters()
#endif
{
	GlobalV::BASIS_TYPE = "lcao";
	GlobalV::global_pseudo_type= "auto";
	GlobalV::PSEUDORCUT = 15.0;
	GlobalV::global_out_dir="./";
	GlobalV::ofs_warning.open("warning.log");
	ofs_running.open("log.txt");

	ucell.latName = "test";
	ucell.ntype = ntype;
	return;
}

#ifdef __LCAO
void test_cell_orb::setup_cell()
#else
void test_cell::setup_cell()
#endif
{
	ucell.setup_cell(
#ifdef __LCAO
	ORB,
#endif
	"./",
	"STRU", 
	ofs_running);

	return;
}

#ifndef __LCAO
void test_cell::count_ntype()
{
	std::cout << "count number of atom types" << std::endl;
	std::ifstream ifs("STRU",std::ios::in);

	if (!ifs)
	{
		std::cout << "ERROR : file STRU does not exist" <<std::endl;
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

		if(x=="LATTICE_CONSTANT" || x=="NUMERICAL_ORBITAL" || x=="LATTICE_VECTORS" || x=="ATOMIC_POSITIONS") break;

		std::string tmpid=x.substr(0,1);
		if(!x.empty() && tmpid!="#") ntype++;
	}

	std::cout << "ntype : "<< ntype << std::endl;
	ifs.close();

	return;
}
#endif

#ifdef __LCAO
void test_cell_orb::setup_kpt()
{
	this->kv.set(
		"KPT",
		GlobalV::NSPIN,
		ucell.G,
		ucell.latvec,
		GAMMA_ONLY_LOCAL,
		ofs_running,
		ofs_warning);
}

void test_cell_orb::alloc_ST()
{
	if(this->GAMMA_ONLY_LOCAL)//gamma only
	{
		Sgamma = new double*[GlobalV::NLOCAL];
		Tgamma = new double*[GlobalV::NLOCAL];
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			Sgamma[i] = new double[GlobalV::NLOCAL];
			Tgamma[i] = new double[GlobalV::NLOCAL];
			for(int j=0; j<GlobalV::NLOCAL; j++)
			{
				Sgamma[i][j] = 0.0;
				Tgamma[i][j] = 0.0;
			}
		}
	}//k-point sampling
	else
	{
		Sk = new std::complex<double>**[this->kv.nkstot];
		Tk = new std::complex<double>**[this->kv.nkstot];
		for(int ik=0; ik<this->kv.nkstot; ik++)
		{
			Sk[ik] = new std::complex<double>*[GlobalV::NLOCAL];
			Tk[ik] = new std::complex<double>*[GlobalV::NLOCAL];
			for(int i=0; i<GlobalV::NLOCAL; i++)
			{
				Sk[ik][i] = new std::complex<double>[GlobalV::NLOCAL];
				Tk[ik][i] = new std::complex<double>[GlobalV::NLOCAL];
				for(int j=0; j<GlobalV::NLOCAL; j++)
				{
					Sk[ik][i][j] = std::complex<double>(0.0,0.0);
					Tk[ik][i][j] = std::complex<double>(0.0,0.0);
				}
			}
		}
	}
}

void test_cell_orb::print_ST()
{

	if(this->GAMMA_ONLY_LOCAL)
	{
		std::ofstream ofs_s("mat_S.dat");
		std::ofstream ofs_t("mat_T.dat");
		for(int i=0; i<GlobalV::NLOCAL; i++)
		{
			for(int j=0; j<GlobalV::NLOCAL; j++)
			{
				int index = i * GlobalV::NLOCAL+j;
				ofs_s << i << " " << j << " " << Sgamma[i][j] << endl;
				ofs_t << i << " " << j << " " << Tgamma[i][j] << endl;
			}
		}
	}
	else
	{
		for(int ik=0;ik<this->kv.nkstot;ik++)
		{
			std::stringstream sss;
			std::stringstream sst;
			sss << "mat_S_k" << ik <<".dat";
			sst << "mat_T_k" << ik <<".dat";
			ofstream ofs_s(sss.str());
			ofstream ofs_t(sst.str());
			for(int i=0; i<GlobalV::NLOCAL; i++)
			{
				for(int j=0; j<GlobalV::NLOCAL; j++)
				{
					int index = i * GlobalV::NLOCAL+j;
					ofs_s << i << " " << j << " " << Sk[ik][i][j] << endl;
					ofs_t << i << " " << j << " " << Tk[ik][i][j] << endl;
				}
			}
			ofs_s.close();
			ofs_t.close();
		}
	}

}

void test_cell_orb::prep_neighbour()
{
    GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
        ofs_running,
        GlobalV::OUT_LEVEL,
        ORB.get_rcutmax_Phi(),
        ORB.get_rcutmax_Beta(),
        GAMMA_ONLY_LOCAL);

    atom_arrange::search(
        GlobalV::SEARCH_PBC,
        ofs_running,
        Test_Cell::GridD,
        ucell,
        GlobalV::SEARCH_RADIUS,
        GlobalV::test_atom_input);
}

///build_ST_new : adapted from build_ST_new in src_lcao
///generates overlap (S) and kinetic (T) matrices
///S_ij = <phi_i | phi_j>, T_ij = <phi_i | T |phi_j>
void test_cell_orb::build_ST_new(const char& dtype)
{
    //array to store data
    double olm[3]={0.0,0.0,0.0};

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>
	ModuleBase::Vector3<double> tau1, tau2, dtau;
	ModuleBase::Vector3<double> dtau1, dtau2, tau0;
    for (int T1=0; T1<ucell.ntype; ++T1)
    {
		Atom* atom1 = &ucell.atoms[T1];
        for (int I1=0; I1<atom1->na; ++I1)
        {
			tau1 = atom1->tau[I1];
            //Test_Cell::GridD.Find_atom(tau1);
            Test_Cell::GridD.Find_atom(ucell, tau1, T1, I1);
            for (int ad = 0; ad < Test_Cell::GridD.getAdjacentNum()+1; ++ad)
            {
                const int T2 = Test_Cell::GridD.getType(ad);
				const int I2 = Test_Cell::GridD.getNatom(ad);
				Atom* atom2 = &ucell.atoms[T2];
				tau2 = Test_Cell::GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();
				if(distance < rcut)
				{
					int iw1_all = ucell.itiaiw2iwt( T1, I1, 0) ; //iw1_all = combined index (it, ia, iw)
					
					for(int jj=0; jj<atom1->nw*GlobalV::NPOL; ++jj)
					{
						const int jj0 = jj/GlobalV::NPOL;
						const int L1 = atom1->iw2l[jj0];
						const int N1 = atom1->iw2n[jj0];
						const int m1 = atom1->iw2m[jj0];

						int iw2_all = ucell.itiaiw2iwt( T2, I2, 0);//zhengdy-soc
						for(int kk=0; kk<atom2->nw*GlobalV::NPOL; ++kk)
						{
							const int kk0 = kk/GlobalV::NPOL;
							const int L2 = atom2->iw2l[kk0];
							const int N2 = atom2->iw2n[kk0];
							const int m2 = atom2->iw2m[kk0];

							olm[0] = olm[1] = olm[2] = 0.0;

							std::complex<double> olm1[4]={ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO, ModuleBase::ZERO};
							std::complex<double> *olm2 = &olm1[0];
							// PLEASE use UOT as an input parameter of this subroutine
							// mohan add 2021-03-30
							OGT.snap_psipsi( ORB, olm, 0, dtype,
									tau1, T1, L1, m1, N1, 
									tau2, T2, L2, m2, N2,
									GlobalV::NSPIN,	olm2//for soc
									);
							if(GAMMA_ONLY_LOCAL)
							{
								if(dtype=='S') Sgamma[iw1_all][iw2_all]+=olm[0];
								if(dtype=='T') Tgamma[iw1_all][iw2_all]+=olm[0];
							}
							else
							{
								for(int ik=0;ik<this->kv.nkstot;ik++)
								{
									ModuleBase::Vector3<double> dR(Test_Cell::GridD.getBox(ad).x, Test_Cell::GridD.getBox(ad).y, Test_Cell::GridD.getBox(ad).z);
									const double arg = ( this->kv.kvec_d[ik] * dR ) * ModuleBase::TWO_PI;
									const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
									if(dtype=='S') Sk[ik][iw1_all][iw2_all]+=olm[0]*kphase;
									if(dtype=='T') Tk[ik][iw1_all][iw2_all]+=olm[0]*kphase;
								}
							}
							++iw2_all;
						}// nw2 
						++iw1_all;
					}// nw1
				}// distance
			}// ad
		}// I1
	}// T1

    return;
}
#endif

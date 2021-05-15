#include "print_info.h"
#include "../src_global/global_variable.h"
#include "../src_pw/global.h"

Print_Info::Print_Info(){}

Print_Info::~Print_Info(){}


void Print_Info::setup_parameters(void)
{
	TITLE("Print_Info","setup_parameters");
	
    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax" || CALCULATION=="nscf"
	        || CALCULATION=="istate" || CALCULATION=="ienvelope" || CALCULATION=="md")
	{
		cout << " ---------------------------------------------------------" << endl;
		if(CALCULATION=="scf")
		{
			cout << " Self-consistent calculations for electrons" << endl;
		}
		else if(CALCULATION=="test")
		{
			cout << " Test run" << endl;
		}
		if(CALCULATION=="relax")
		{
            cout << " Ion relaxation calculations" << endl;
		}
        if(CALCULATION=="cell-relax")
        {
            cout << " Cell relaxation calculations" << endl;
        }
		if(CALCULATION=="md")
		{
			cout << " Molecular Dynamics simulations" << endl;

			cout << " ---------------------------------------------------------" << endl;

			if(INPUT.mdp.mdtype ==1 || INPUT.mdp.mdtype==2)
			{
				cout << " ENSEMBLE                 : " << "NVT" << endl;
				cout << " Qmass for NVT(a.u.)      : " << INPUT.mdp.Qmass/6.02/9.109*1e5 << endl;
			}
			else if(INPUT.mdp.mdtype==0)
			{
				cout << " ENSEMBLE                 : " << "NVE" << endl;
			}
			
			cout << " Time interval(fs)        : " << INPUT.mdp.dt << endl;
		}
		cout << " ---------------------------------------------------------" << endl;


		cout << " " << setw(8) << "SPIN"
		     << setw(16) << "KPOINTS"
		     << setw(12) << "PROCESSORS";

		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
		{
			cout << setw(12) << "NBASE";
		}

		cout << endl;



		cout << " " << setw(8) << NSPIN;

		if(GAMMA_ONLY_LOCAL)
		{
			if(COLOUR && MY_RANK==0)
			{
				// red
				//printf( "\e[31m%-16s\e[0m", "Gamma");
				printf( "[31m%-16s[0m", "Gamma");
			}
			else
			{
				cout << setw(16) << "Gamma";
			}
		}
		else
		{
			if(COLOUR && MY_RANK==0)
			{
				// zi
				//printf( "\e[35m%-16d\e[0m", kv.nkstot);
				printf( "[35m%-16d[0m", kv.nkstot);
			}
			else
			{
				cout << setw(16) << kv.nkstot;
			}
		}

		cout << setw(12) << NPROC;

		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
		{
			cout << setw(12) << NLOCAL;
		}

		cout << endl;




		cout << " ---------------------------------------------------------" << endl;
		if(BASIS_TYPE=="lcao") 
		{
			if(COLOUR && MY_RANK==0)
			{
				string a = "Use Systematically Improvable Atomic bases";
				//printf( " \e[36m%-45s\e[0m\n", a.c_str());
				printf( " [36m%-45s[0m\n", a.c_str());
			}
			else
			{
				cout << " Use Systematically Improvable Atomic bases" << endl;
			}
		}
		else if(BASIS_TYPE=="lcao_in_pw")
		{
			cout << " Expand Atomic bases into plane waves" << endl;
		}
		else if(BASIS_TYPE=="pw")
		{
			cout << " Use plane wave basis" << endl;
		}
		cout << " ---------------------------------------------------------" << endl;



		//----------------------------------
		// second part
		//----------------------------------

		cout << " " << setw(8) << "ELEMENT";

		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
		{
			cout << setw(16) << "ORBITALS";
			cout << setw(12) << "NBASE";
		}
		cout << setw(12) << "NATOM";

		cout << setw(12) << "XC";
		cout << endl;



		for(int it=0; it<ucell.ntype; ++it)
		{
			if(COLOUR && MY_RANK==0)
			{
				printf( " [36m%-8s[0m", ucell.atoms[it].label.c_str());
			}
			else
			{
				cout << " " << setw(8) << ucell.atoms[it].label;
			}

			if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw")
			{
				stringstream orb;

				int norb = 0;
				/*for(int L=0; L<=ORB.Phi[it].getLmax(); ++L)
				{
					norb += (2*L+1)*ORB.Phi[it].getNchi(L);
					orb << ORB.Phi[it].getNchi(L);
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << ORB.Phi[it].getRcut() << "au";*/


				for(int L=0; L<=ucell.atoms[it].nwl; ++L)        // pengfei Li 16-2-29
				{
					norb += (2*L+1)* ucell.atoms[it].l_nchi[L];
					orb << ucell.atoms[it].l_nchi[L];
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << ucell.atoms[it].Rcut << "au";

				if(COLOUR && MY_RANK==0)
				{
					//printf( "\e[36m%-16s\e[0m", orb.str().c_str());
					//printf( "\e[36m%-12d\e[0m", norb);
					printf( "[36m%-16s[0m", orb.str().c_str());
					printf( "[36m%-12d[0m", norb);
				}
				else
				{
					cout << setw(16) << orb.str();
					cout << setw(12) << norb;
				}
			}


			cout << setw(12) << ucell.atoms[it].na;

//				if(ucell.atoms[it].dft[1]=="PZ")    // pengfei Li added 2015-1-31 cancelled by zws
//				{
//					//cout << setw(12) << "PZ-LDA";
//
//				}
//				else
//				{
//					//cout << setw(12) << ucell.atoms[it].dft[0];
//                                        cout << setw(12) << "PBE";
//				}
			xcf.ostreamdft(cout); // zws add 20150108
			//cout << " ( "  << setw(3) << xcf.iexch << setw(3) << xcf.icorr << setw(3) << xcf.igcx << setw(3) << xcf.igcc << ")";
			cout << endl;
		}

		cout << " ---------------------------------------------------------" << endl;
		cout << " Initial plane wave basis and FFT box" << endl;
		cout << " ---------------------------------------------------------" << endl;

	}

	return;
}

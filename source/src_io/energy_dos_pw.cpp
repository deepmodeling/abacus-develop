#include "dos.h"
#include "../src_pw/global.h"
#include "../src_pw/energy.h"

void energy::perform_dos_pw(void)
{
	TITLE("energy","perform_dos_pw");

	if(out_dos !=0 || out_band !=0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        GlobalV::ofs_running << " |                                                                    |" << endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << endl;
        GlobalV::ofs_running << " | done here.                                                         |" << endl;
        GlobalV::ofs_running << " |                                                                    |" << endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }

	//qianrui modify 2020-10-18
	if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax")
	{
		stringstream ss;
		ss << GlobalV::global_out_dir << "istate.info" ;
		if(GlobalV::MY_RANK==0)
		{
			ofstream ofsi( ss.str().c_str() ); // clear istate.info
			ofsi.close();
		}
#ifdef __MPI
		for(int ip=0; ip<GlobalV::NPOOL; ip++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if( GlobalV::MY_POOL == ip )
			{
				if( GlobalV::RANK_IN_POOL != 0 ) continue;
#endif
				ofstream ofsi2( ss.str().c_str(), ios::app );
				if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
				{
					for (int ik = 0;ik < GlobalC::kv.nks;ik++)
					{
						ofsi2<<"BAND"
						<<setw(25)<<"Energy(ev)"
						<<setw(25)<<"Occupation"
#ifdef __MPI
						<<setw(25)<<"Kpoint = "<<Pkpoints.startk_pool[ip]+ik+1
#else
						<<setw(25)<<"Kpoint = "<<ik+1
#endif
						<<setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<endl;
						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<setw(6)<<ib+1<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV<<setw(25)<<wf.wg(ik,ib)<<endl;
						}
						ofsi2 <<endl;
						ofsi2 <<endl;
					}
				}
				else
				{
					for (int ik = 0;ik < GlobalC::kv.nks/2;ik++)
					{
						ofsi2<<"BAND"
						<<setw(25)<<"Spin up Energy(ev)"
						<<setw(25)<<"Occupation"
						<<setw(25)<<"Spin down Energy(ev)"
						<<setw(25)<<"Occupation"
#ifdef __MPI
						<<setw(25)<<"Kpoint = "<<Pkpoints.startk_pool[ip]+ik+1
#else
						<<setw(25)<<"Kpoint = "<<ik+1
#endif
						<<setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<endl;

						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<setw(6)<<ib+1
							<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV
							<<setw(25)<<wf.wg(ik,ib)
							<<setw(25)<<wf.ekb[(ik+GlobalC::kv.nks/2)][ib]* Ry_to_eV
							<<setw(25)<<wf.wg(ik+GlobalC::kv.nks/2,ib)<<endl;
						}
						ofsi2 <<endl;
						ofsi2 <<endl;

					}
				}

				ofsi2.close();
#ifdef __MPI
			}
		}
#endif
	}
	
	int nspin0=1;
	if(GlobalV::NSPIN==2) nspin0=2;

	if(this->out_dos)
	{
//find energy range
		double emax = wf.ekb[0][0];
		double emin = wf.ekb[0][0];
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				emax = std::max( emax, wf.ekb[ik][ib] );
				emin = std::min( emin, wf.ekb[ik][ib] );
			}
		}

#ifdef __MPI
		Parallel_Reduce::gather_max_double_all(emax);
		Parallel_Reduce::gather_min_double_all(emin);
#endif

		emax *= Ry_to_eV;
		emin *= Ry_to_eV;

//scale up a little bit so the end peaks are displaced better
		const double scl=this->dos_scale;
		double delta=(emax-emin)*scl;
		cout << scl;
		emax=emax+delta/2.0;
		emin=emin-delta/2.0;

				OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
				OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
		// 		atom_arrange::set_sr_NL();
		//		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21
		
//determine #. energy points	
		const double de_ev = this->dos_edelta_ev;
		cout << de_ev;

		const int npoints = static_cast<int>(std::floor ( ( emax - emin ) / de_ev ));
		const int np=npoints;

	 	for(int is=0; is<nspin0; ++is)
	 	{
//DOS_ispin contains not smoothed dos
			 stringstream ss;
			 ss << GlobalV::global_out_dir << "DOS" << is+1;

			 Dos::calculate_dos(
					 is,
					 GlobalC::kv.isk,
					 ss.str(), 
					 this->dos_edelta_ev, 
					 emax, 
					 emin, 
					 GlobalC::kv.nks, GlobalC::kv.nkstot, GlobalC::kv.wk, wf.wg, GlobalV::NBANDS, wf.ekb );
			 ifstream in(ss.str().c_str());
			 if(!in)
			 {
				       //cout<<"\n Can't find file : "<< name << endl;
				       //return 0;
			 }

			 //----------------------------------------------------------
			 // FOUND LOCAL VARIABLES :
			 // NAME : number(number of DOS points)
			 // NAME : nk(number of k point used)
			 // NAME : energy(energy range,from emin_ev to emax_ev)
			 // NAME : dos(old,count k points in the energy range)
			 // NAME : dos2(new,count k points in the energy range)
			 //----------------------------------------------------------
			 int number=0;
			 int nk=0;
			 in >> number;
			 in >> nk;
			 double *energy = new double[number];
			 double *dos = new double[number];
			 double *dos2 = new double[number];
			 for(int i=0 ;i<number; i++)
			 {
				 energy[i] = 0.0;
				 dos[i] = 0.0;
				 dos2[i] =0.0;
			 }

			 for(int i=0;i<number;i++)
			 {
				 in >> energy[i] >> dos[i];
			 }
			 if(!in.eof())
			 {
				 //cout<<"\n Read Over!"<<endl;
			 }
			 in.close();

//now use Gaussian smearing to smooth the dos and write to DOS_is_smearing

			 //----------------------------------------------------------
			 // EXPLAIN : b is an empirical value.
			 // please DIY b!!
			 //----------------------------------------------------------

			 //double b = INPUT.b_coef;
			 double b = bcoeff;
			 for(int i=0;i<number;i++)
			 {
				 double Gauss=0.0;
	
				 for(int j=0;j<number;j++)
				 {
					 double de = energy[j] - energy[i];
					 double de2 = de * de;
					 //----------------------------------------------------------
					 // EXPLAIN : if en
					 //----------------------------------------------------------
					 Gauss = exp(-de2/b/b)/sqrt(3.1415926)/b;
					 dos2[j] += dos[i]*Gauss;
				 }
			 }

		 //----------------------------------------------------------
		 // EXPLAIN : output DOS2.txt
		 //----------------------------------------------------------
		 stringstream sss;
		 sss << GlobalV::global_out_dir << "DOS" << is+1 << "_smearing" << ".dat" ;
		 ofstream out(sss.str().c_str());
		 double sum2=0.0;
		 for(int i=0;i<number;i++)
		 {
			 sum2 += dos2[i];
			 //            if(dos2[i]<1e-5)
			 //            {
			 //                    dos2[i] = 0.00;
			 //            }
			 out <<setw(20)<<energy[i]
				 <<setw(20)<<dos2[i]
				 <<setw(20)<<sum2<<"\n";
		 }
		 out.close();

		 //----------------------------------------------------------
		 // DELETE
		 //----------------------------------------------------------
		 delete[] dos;
		 delete[] dos2;
		 delete[] energy;

		 //cout<<" broden spectrum over, success : ) "<<endl;

	 }

	}//out_dos=1
	if(this->out_band) //pengfei 2014-10-13
	{
		int nks=0;
		if(nspin0==1) 
		{
			nks = GlobalC::kv.nkstot;
		}
		else if(nspin0==2) 
		{
			nks = GlobalC::kv.nkstot/2;
		}

		for(int is=0; is<nspin0; is++)
		{
			stringstream ss2;
			ss2 << GlobalV::global_out_dir << "BANDS_" << is+1 << ".dat";
			GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << endl;
			Dos::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, this->ef*0, wf.ekb);
		}

	}
}


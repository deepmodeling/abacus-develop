#include "dos.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"
#include "../src_pw/energy.h"
#include "../src_pw/wavefunc.h"
#ifdef __LCAO
#include "mulliken_charge.h"
#include "../src_lcao/LCAO_nnr.h"
#include "../src_lcao/LCAO_gen_fixedH.h"    
#include "../src_lcao/local_orbital_charge.h"
#include "../src_lcao/LCAO_matrix.h"
#include "../src_lcao/global_fp.h"
#include "../src_lcao/wfc_dm_2d.h"
#include "../module_neighbor/sltk_atom_arrange.h"//qifeng-2019-01-21
#endif
#include "../module_base/lapack_connector.h"
#include "../module_base/scalapack_connector.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include <vector>
#ifdef __MPI
#include<mpi.h>
#endif
#include <sys/time.h>

void energy::perform_dos(void)
{
	ModuleBase::TITLE("energy","perform_dos");


    if(out_dos !=0 || out_band !=0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }


/*	if(GlobalV::MY_RANK==0)
	{
		if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax")
		{
			std::stringstream ss;
			ss << GlobalV::global_out_dir << "istate.info" ;
			std::ofstream ofsi( ss.str().c_str() );
			*for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				ofsi << "0 " << ib+1;
				for(int is=0; is<GlobalV::NSPIN; ++is)
				{
					ofsi << " " << GlobalC::wf.ekb[is][ib];
					}
					ofsi << std::endl;
					}
			// pengfei 2015-4-1
			if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
			{ 
				for (int ik = 0;ik < GlobalC::kv.nks;ik++)
				{                       
					ofsi<<"BAND"
					<<std::setw(25)<<"Energy(ev)"
					<<std::setw(25)<<"Occupation"
					<<std::setw(25)<<"Kpoint = "<<ik+1
					<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;
					for(int ib=0;ib<GlobalV::NBANDS;ib++)
					{
						ofsi<<ib+1<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* Ry_to_eV<<std::setw(25)<<GlobalC::wf.wg(ik,ib)<<std::endl;
					}
					ofsi <<std::endl; 
					ofsi <<std::endl;                              
				}
			}
			else
			{
				for (int ik = 0;ik < GlobalC::kv.nks/2;ik++)
				{
					ofsi<<"BAND"<<std::setw(25)<<"Spin up Energy(ev)"
					<<std::setw(25)<<"Occupation"
					<<std::setw(25)<<"Spin down Energy(ev)"
					<<std::setw(25)<<"Occupation"
					<<std::setw(25)<<"Kpoint = "<<ik+1
					<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;
					for(int ib=0;ib<GlobalV::NBANDS;ib++)
					{
						ofsi<<ib+1<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* Ry_to_eV
						<<std::setw(25)<<GlobalC::wf.wg(ik,ib)
						<<std::setw(25)<<GlobalC::wf.ekb[(ik+GlobalC::kv.nks/2)][ib]* Ry_to_eV
						<<std::setw(25)<<GlobalC::wf.wg(ik+GlobalC::kv.nks/2,ib)<<std::endl;
					}
					ofsi <<std::endl;
					ofsi <<std::endl;

				}
			}

			ofsi.close();
		}
	}
*/

	//qianrui modify 2020-10-18
	if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax")
	{
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "istate.info" ;
		if(GlobalV::MY_RANK==0)
		{
			std::ofstream ofsi( ss.str().c_str() ); // clear istate.info
			ofsi.close();
		}
		for(int ip=0; ip<GlobalV::NPOOL; ip++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if( GlobalV::MY_POOL == ip )
			{
				if( GlobalV::RANK_IN_POOL != 0 ) continue;
				std::ofstream ofsi2( ss.str().c_str(), ios::app );
				if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
				{
					for (int ik = 0;ik < GlobalC::kv.nks;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Energy(ev)"
						<<std::setw(25)<<"Occupation"
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;
						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* ModuleBase::Ry_to_eV<<std::setw(25)<<GlobalC::wf.wg(ik,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;
					}
				}
				else
				{
					for (int ik = 0;ik < GlobalC::kv.nks/2;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Spin up Energy(ev)"
						<<std::setw(25)<<"Occupation"
						<<std::setw(25)<<"Spin down Energy(ev)"
						<<std::setw(25)<<"Occupation"
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;

						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1
							<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* ModuleBase::Ry_to_eV
							<<std::setw(25)<<GlobalC::wf.wg(ik,ib)
							<<std::setw(25)<<GlobalC::wf.ekb[(ik+GlobalC::kv.nks/2)][ib]* ModuleBase::Ry_to_eV
							<<std::setw(25)<<GlobalC::wf.wg(ik+GlobalC::kv.nks/2,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;

					}
				}

				ofsi2.close();
			}
		}
	}


	// GlobalV::mulliken charge analysis
#ifdef __LCAO
	if(GlobalV::mulliken == 1)
	{
		Mulliken_Charge   MC;
		MC.stdout_mulliken();			
	}//qifeng add 2019/9/10
#endif

	int nspin0=1;
	if(GlobalV::NSPIN==2) nspin0=2;

	if(this->out_dos)
	{
		// find the maximal and minimal band energy.
		double emax = GlobalC::wf.ekb[0][0];
		double emin = GlobalC::wf.ekb[0][0];
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				emax = std::max( emax, GlobalC::wf.ekb[ik][ib] );
				emin = std::min( emin, GlobalC::wf.ekb[ik][ib] );
			}
		}

#ifdef __MPI
		Parallel_Reduce::gather_max_double_all(emax);
		Parallel_Reduce::gather_min_double_all(emin);
#endif

		emax *= ModuleBase::Ry_to_eV;
		emin *= ModuleBase::Ry_to_eV;

//scale up a little bit so the end peaks are displaced better
		double delta=(emax-emin)*dos_scale;
		std::cout << dos_scale;
		emax=emax+delta/2.0;
		emin=emin-delta/2.0;

			//	OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
			//	OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
		//  output the PDOS file.////qifeng-2019-01-21
		// 		atom_arrange::set_sr_NL();
		//		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21
		const double de_ev = this->dos_edelta_ev;


		const int npoints = static_cast<int>(std::floor ( ( emax - emin ) / de_ev ));

		int NUM=GlobalV::NLOCAL*npoints;
		Wfc_Dm_2d D;
		D.init();
		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			for(int in=0;in<GlobalV::NSPIN;in++)
			{

				D.wfc_gamma[in]=GlobalC::LOC.wfc_dm_2d.wfc_gamma[in];
			}

		}
		else 
		{

			for(int in=0;in<GlobalC::kv.nks;in++)
			{

				D.wfc_k[in] = GlobalC::LOC.wfc_dm_2d.wfc_k[in];
			}
		}


		const int np=npoints;
		ModuleBase::matrix*  pdosk = new ModuleBase::matrix[nspin0];

		for(int is=0; is<nspin0; ++is)
		{


			pdosk[is].create(GlobalV::NLOCAL,np,true);



		}
		ModuleBase::matrix*  pdos = new ModuleBase::matrix[nspin0];
		for(int is=0; is<nspin0; ++is)
		{
			pdos[is].create(GlobalV::NLOCAL,np,true);

		}

		double a = bcoeff;
		double c=2*3.141592653;
		double b = sqrt(c)*a;                                         

		std::complex<double>*waveg = new std::complex<double>[GlobalV::NLOCAL];

		double*Gauss = new double[np];

		for(int is=0; is<nspin0; ++is)
		{
			if(GlobalV::GAMMA_ONLY_LOCAL)
			{
				std::vector<ModuleBase::matrix>   Mulk;
				Mulk.resize(1);
				Mulk[0].create(GlobalC::ParaO.ncol,GlobalC::ParaO.nrow);


				ModuleBase::matrix Dwf = D.wfc_gamma[is];
				for (int i=0; i<GlobalV::NBANDS; ++i)		  
				{     
					ModuleBase::GlobalFunc::ZEROS(waveg, GlobalV::NLOCAL);

					ModuleBase::GlobalFunc::ZEROS(Gauss,np);
					for (int n=0; n<npoints; ++n)		  
					{  
						double en=emin+n * de_ev;
						double en0=GlobalC::wf.ekb[0][i]*ModuleBase::Ry_to_eV;
						double de = en-en0;
						double de2 = 0.5*de * de;
						Gauss[n] = GlobalC::kv.wk[0]*exp(-de2/a/a)/b;
					}

					const int NB= i+1;

					const double one_float=1.0, zero_float=0.0;
					const int one_int=1;


					const char T_char='T';		
					pdgemv_(
							&T_char,
							&GlobalV::NLOCAL,&GlobalV::NLOCAL,
							&one_float,
							GlobalC::LM.Sloc, &one_int, &one_int, GlobalC::ParaO.desc,
							Dwf.c, &one_int, &NB, GlobalC::ParaO.desc, &one_int,
							&zero_float,
							Mulk[0].c, &one_int, &NB, GlobalC::ParaO.desc,
							&one_int);

					for (int j=0; j<GlobalV::NLOCAL; ++j)
					{

						if ( GlobalC::ParaO.in_this_processor(j,i) )
						{

							const int ir = GlobalC::ParaO.trace_loc_row[j];
							const int ic = GlobalC::ParaO.trace_loc_col[i];
							waveg[j] = Mulk[0](ic,ir)*D.wfc_gamma[is](ic,ir);
							const double x = waveg[j].real();
							LapackConnector::axpy(np , x,Gauss, 1,pdosk[is].c+j*pdosk[is].nc,1);
						}
					} 
				}//ib
			}//if
			else
			{
				GlobalV::SEARCH_RADIUS = atom_arrange::set_sr_NL(
					GlobalV::ofs_running,
					GlobalV::OUT_LEVEL,
					GlobalC::ORB.get_rcutmax_Phi(), 
					GlobalC::ucell.infoNL.get_rcutmax_Beta(), 
					GlobalV::GAMMA_ONLY_LOCAL);

				atom_arrange::search(
					GlobalV::SEARCH_PBC,
					GlobalV::ofs_running,
					GlobalC::GridD, 
					GlobalC::ucell, 
					GlobalV::SEARCH_RADIUS, 
					GlobalV::test_atom_input);//qifeng-2019-01-21

				// mohan update 2021-04-16
				GlobalC::LOWF.orb_con.set_orb_tables(
						GlobalV::ofs_running,
						GlobalC::UOT, 
						GlobalC::ORB,
						GlobalC::ucell.ntype,
						GlobalC::ucell.lmax,
						INPUT.lcao_ecut,
						INPUT.lcao_dk,
						INPUT.lcao_dr,
						INPUT.lcao_rmax, 
						GlobalC::ucell.lat0, 
						INPUT.out_descriptor,
						INPUT.out_r_matrix,
						Exx_Abfs::Lmax,
						GlobalV::FORCE,
						GlobalV::MY_RANK,
						GlobalC::ucell.infoNL.nprojmax,
						GlobalC::ucell.infoNL.nproj,
						GlobalC::ucell.infoNL.Beta);

				GlobalC::LM.allocate_HS_R(GlobalC::LNNR.nnr);
				GlobalC::LM.zeros_HSR('S', GlobalC::LNNR.nnr);
				GlobalC::UHM.genH.calculate_S_no();
				GlobalC::UHM.genH.build_ST_new('S', false, GlobalC::ucell);
				std::vector<ModuleBase::ComplexMatrix> Mulk;
				Mulk.resize(1);
				Mulk[0].create(GlobalC::ParaO.ncol,GlobalC::ParaO.nrow);


				for(int ik=0;ik<GlobalC::kv.nks;ik++)
				{

					if(is == GlobalC::kv.isk[ik])
					{
						GlobalC::LM.allocate_HS_k(GlobalC::ParaO.nloc);
						GlobalC::LM.zeros_HSk('S');
						GlobalC::LNNR.folding_fixedH(ik);


						ModuleBase::ComplexMatrix Dwfc = conj(D.wfc_k[ik]);

						for (int i=0; i<GlobalV::NBANDS; ++i)		  
						{     

							ModuleBase::GlobalFunc::ZEROS(waveg, GlobalV::NLOCAL);


							ModuleBase::GlobalFunc::ZEROS(Gauss,np);
							for (int n=0; n<npoints; ++n)		  
							{  
								double en=emin+n * de_ev;
								double en0=GlobalC::wf.ekb[ik][i]*ModuleBase::Ry_to_eV;
								double de = en-en0;
								double de2 = 0.5*de * de;
								Gauss[n] = GlobalC::kv.wk[ik]*exp(-de2/a/a)/b;
							}

							const int NB= i+1;

							const double one_float=1.0, zero_float=0.0;
							const int one_int=1;
							//   const int two_int=2;
							const char T_char='T';		// N_char='N',U_char='U'

							pzgemv_(
									&T_char,
									&GlobalV::NLOCAL,&GlobalV::NLOCAL,
									&one_float,
									GlobalC::LM.Sloc2, &one_int, &one_int, GlobalC::ParaO.desc,
									Dwfc.c, &one_int, &NB, GlobalC::ParaO.desc, &one_int,
									&zero_float,
									Mulk[0].c, &one_int, &NB, GlobalC::ParaO.desc,
									&one_int);



							for (int j=0; j<GlobalV::NLOCAL; ++j)
							{

								if ( GlobalC::ParaO.in_this_processor(j,i) )
								{

									const int ir = GlobalC::ParaO.trace_loc_row[j];
									const int ic = GlobalC::ParaO.trace_loc_col[i];

									waveg[j] = Mulk[0](ic,ir)*D.wfc_k[ik](ic,ir);
									const double x = waveg[j].real();
									LapackConnector::axpy(np , x,Gauss, 1,pdosk[is].c+j*pdosk[is].nc,1);

								}
							}                             


						}//ib

					}//if                       
				}//ik
#ifdef __MPI
				atom_arrange::delete_vector(
					GlobalV::ofs_running,
					GlobalV::SEARCH_PBC, 
					GlobalC::GridD, 
					GlobalC::ucell, 
					GlobalV::SEARCH_RADIUS, 
					GlobalV::test_atom_input);
#endif
				// mohan update 2021-02-10
				GlobalC::LOWF.orb_con.clear_after_ions(GlobalC::UOT, GlobalC::ORB, INPUT.out_descriptor, GlobalC::ucell.infoNL.nproj);
			}//else

		 MPI_Reduce(pdosk[is].c, pdos[is].c , NUM , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);
	 }//is                                              
	 delete[] pdosk;                                               
	 delete[] waveg;
	 if(GlobalV::MY_RANK == 0)
	 {
		 {  std::stringstream ps;
			 ps << GlobalV::global_out_dir << "TDOS";
			 std::ofstream out(ps.str().c_str());
			 if (GlobalV::NSPIN==1)
			 {

				 for (int n=0; n<npoints; ++n)
				 { double y=0.0;
					 double en=emin + n * de_ev;
					 for (int i=0; i<GlobalV::NLOCAL; i++)
					 {
						 y +=  pdos[0](i,n);
					 }  

					 out <<std::setw(20)<< en <<std::setw(30)<< y << std::endl;
				 }
			 }
			 else if (GlobalV::NSPIN==2)
			 {
				 for (int n=0; n<npoints; ++n)
				 { double y=0.0;
					 double z=0.0;
					 double en=emin + n * de_ev;
					 for (int i=0; i<GlobalV::NLOCAL; i++)
					 {
						 y +=  pdos[0](i,n);
						 z +=  pdos[1](i,n);

					 }  

					 out <<std::setw(20)<< en <<std::setw(30)<< y << std::setw(30)<< z<< std::endl;
				 }
			 }
			 out.close();
		 }

		 std::string Name_Angular[5][11];
		 /* decomposed Mulliken charge */

		 Name_Angular[0][0] = "s          ";
		 Name_Angular[1][0] = "px         ";
		 Name_Angular[1][1] = "py         ";
		 Name_Angular[1][2] = "pz         ";
		 Name_Angular[2][0] = "d3z^2-r^2  ";
		 Name_Angular[2][1] = "dxy        ";
		 Name_Angular[2][2] = "dxz        ";
		 Name_Angular[2][3] = "dx^2-y^2   ";
		 Name_Angular[2][4] = "dyz        ";
		 Name_Angular[3][0] = "f5z^2-3r^2 ";
		 Name_Angular[3][1] = "f5xz^2-xr^2";
		 Name_Angular[3][2] = "f5yz^2-yr^2";
		 Name_Angular[3][3] = "fzx^2-zy^2 ";
		 Name_Angular[3][4] = "fxyz       ";
		 Name_Angular[3][5] = "fx^3-3*xy^2";
		 Name_Angular[3][6] = "f3yx^2-y^3 ";
		 Name_Angular[4][0] = "g1         ";
		 Name_Angular[4][1] = "g2         ";
		 Name_Angular[4][2] = "g3         ";
		 Name_Angular[4][3] = "g4         ";
		 Name_Angular[4][4] = "g5         ";
		 Name_Angular[4][5] = "g6         ";
		 Name_Angular[4][6] = "g7         ";
		 Name_Angular[4][7] = "g8         ";
		 Name_Angular[4][8] = "g9         ";

		 {std::stringstream as;
			 as << GlobalV::global_out_dir << "PDOS";
			 std::ofstream out(as.str().c_str());

			 out << "<"<<"pdos"<<">" <<std::endl;
			 out << "<"<<"nspin"<<">" << GlobalV::NSPIN<< "<"<<"/"<<"nspin"<<">"<< std::endl;
			 out << "<"<<"norbitals"<<">" <<std::setw(2) <<GlobalV::NLOCAL<< "<"<<"/"<<"norbitals"<<">"<< std::endl;
			 out << "<"<<"energy"<<"_"<<"values units"<<"="<<"\""<<"eV"<<"\""<<">"<<std::endl;

			 for (int n=0; n<npoints; ++n)
			 { double y=0.0;
				 double en=emin + n * de_ev;
				 out <<std::setw(20)<< en << std::endl;
			 }
			 out << "<"<<"/"<<"energy"<<"_"<<"values"<<">" <<std::endl;
			 for (int i=0; i<GlobalC::ucell.nat; i++)
			 {   
				 int a = GlobalC::ucell.iat2ia[i];
				 int t = GlobalC::ucell.iat2it[i];
				 Atom* atom1 = &GlobalC::ucell.atoms[t];
				 for(int j=0; j<atom1->nw; ++j)
				 {
					 const int L1 = atom1->iw2l[j];
					 const int N1 = atom1->iw2n[j];
					 const int m1 = atom1->iw2m[j];
					 const int w = GlobalC::ucell.itiaiw2iwt(t, a, j);

					 //out << "<"<<"/"<<"energy"<<"_"<<"values"<<">" <<std::endl;
					 out << "<"<<"orbital" <<std::endl;
					 out <<std::setw(6)<< "index"<<"="<<"\""<<std::setw(40) <<w+1<<"\""<<std::endl;
					 out <<std::setw(5)<< "atom"<<"_"<<"index"<<"="<<"\""<<std::setw(40) <<i+1<<"\""<<std::endl;
					 out <<std::setw(8)<< "species"<<"="<<"\""<<GlobalC::ucell.atoms[t].label<<"\""<<std::endl;
					 out<<std::setw(2)<< "l"<<"="<<"\""<<std::setw(40)<<L1<<"\""<<std::endl;
					 out <<std::setw(2)<< "m"<<"="<<"\""<<std::setw(40)<<m1<<"\""<<std::endl;
					 out <<std::setw(2)<< "z"<<"="<<"\""<<std::setw(40)<<N1+1<<"\""<<std::endl;
					 out << ">" <<std::endl;
					 out << "<"<<"data"<<">" <<std::endl;
					 if (GlobalV::NSPIN==1)
					 {
						 for (int n=0; n<npoints; ++n)
						 {


							 out <<std::setw(13)<< pdos[0](w,n)<<std::endl;
						 }//n
					 }
					 else if (GlobalV::NSPIN==2)
					 {
						 for (int n=0; n<npoints; ++n)
						 {
							 out <<std::setw(20)<< pdos[0](w,n)<< std::setw(30)<< pdos[1](w,n)<<std::endl;
						 }//n
					 }

					 out << "<"<<"/"<<"data"<<">" <<std::endl;
					 out << "<"<<"/"<<"orbital"<<">" <<std::endl;
				 }//j
			 }//i

			 out << "<"<<"/"<<"pdos"<<">" <<std::endl;
			 out.close();}
		 {  std::stringstream os;
			 os<<GlobalV::global_out_dir<<"Orbital";
			 std::ofstream out(os.str().c_str());
			 out<< std::setw(5)<<"io"<< std::setw(8) <<"spec" <<std::setw(5)<<"l"<<std::setw(5)<<"m"<<std::setw(5)<<"z"<<std::setw(5)<<"sym"<<std::endl;


			 for (int i=0; i<GlobalC::ucell.nat; i++)
			 {
				 int   t = GlobalC::ucell.iat2it[i];
				 Atom* atom1 = &GlobalC::ucell.atoms[t];  
				 for(int j=0; j<atom1->nw; ++j)
				 {
					 const int L1 = atom1->iw2l[j];
					 const int N1 = atom1->iw2n[j];
					 const int m1 = atom1->iw2m[j];
					 out <<std::setw(5) << i << std::setw(8) 
						<< GlobalC::ucell.atoms[t].label <<std::setw(5)
							<<L1<<std::setw(5) <<m1<<std::setw(5)<<N1+1<<std::setw(15)<< Name_Angular[L1][m1] << std::endl;
				 }
			 }
			 out <<std::endl<<std::endl;
			 out <<std::setw(5)<< "io"<<std::setw(2)<<"="<<std::setw(2)<<"Orbital index in supercell"<<std::endl;
			 out <<std::setw(5)<< "spec"<<std::setw(2)<<"="<<std::setw(2)<<"Atomic species label"<<std::endl;
			 out <<std::setw(5)<< "l"<<std::setw(2)<<"="<<std::setw(2)<<"Angular mumentum quantum number"<<std::endl;
			 out <<std::setw(5)<< "m"<<std::setw(2)<<"="<<std::setw(2)<<"Magnetic quantum number"<<std::endl;
			 out <<std::setw(5)<< "z"<<std::setw(2)<<"="<<std::setw(2)<<"Zeta index of orbital"<<std::endl;
			 out <<std::setw(5)<< "sym"<<std::setw(2)<<"="<<std::setw(2)<<"Symmetry name of real orbital"<<std::endl;
			 out.close();}

	 }       
	 delete[] pdos;

	 // output the DOS file.
	 for(int is=0; is<nspin0; ++is)
	 {
		 std::stringstream ss;
		 ss << GlobalV::global_out_dir << "DOS" << is+1;

		 Dos::calculate_dos(
				 is,
				 GlobalC::kv.isk,
				 ss.str(), 
				 this->dos_edelta_ev, 
				 emax, 
				 emin, 
				 GlobalC::kv.nks, GlobalC::kv.nkstot, GlobalC::kv.wk, GlobalC::wf.wg, GlobalV::NBANDS, GlobalC::wf.ekb );
		 std::ifstream in(ss.str().c_str());
		 if(!in)
		 {
			 //      std::cout<<"\n Can't find file : "<< name << std::endl;
			 //      return 0;
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
			 //std::cout<<"\n Read Over!"<<std::endl;
		 }
		 in.close();

		 //----------------------------------------------------------
		 // EXPLAIN : b is an empirical value.
		 // please DIY b!!
		 //----------------------------------------------------------

		 //double b = INPUT.b_coef;
		 double b = sqrt(2.0)*bcoeff;
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
		 std::stringstream sss;
		 sss << GlobalV::global_out_dir << "DOS" << is+1 << "_smearing" << ".dat" ;
		 std::ofstream out(sss.str().c_str());
		 double sum2=0.0;
		 for(int i=0;i<number;i++)
		 {
			 sum2 += dos2[i];
			 //            if(dos2[i]<1e-5)
			 //            {
			 //                    dos2[i] = 0.00;
			 //            }
			 out <<std::setw(20)<<energy[i]
				 <<std::setw(20)<<dos2[i]
				 <<std::setw(20)<<sum2<<"\n";
		 }
		 out.close();

		 //----------------------------------------------------------
		 // DELETE
		 //----------------------------------------------------------
		 delete[] dos;
		 delete[] dos2;
		 delete[] energy;

		 //std::cout<<" broden spectrum over, success : ) "<<std::endl;

	 }


		// GlobalV::mulliken charge analysis
		if(out_dos == 2)
		{
			std::stringstream sp;
			sp << GlobalV::global_out_dir << "Mulliken.dat";
			Dos::calculate_Mulliken(sp.str());
		}
	
		if(nspin0==1)
		{
			GlobalV::ofs_running << " Fermi energy is " << this->ef << " Rydberg" << std::endl;
		}
		else if(nspin0==2)
		{
			GlobalV::ofs_running << " Fermi energy (spin = 1) is " << this->ef_up << " Rydberg" << std::endl;
			GlobalV::ofs_running << " Fermi energy (spin = 2) is " << this->ef_dw << " Rydberg" << std::endl;
		}

		//int nks;
		//if(nspin0==1) nks = GlobalC::kv.nkstot;
		//else if(nspin0==2) nks = GlobalC::kv.nkstot/2;



		/*for(int is=0; is<GlobalV::NSPIN; is++)
		{
			std::stringstream ss2;
			ss2 << GlobalV::global_out_dir << "BANDS_" << is+1 << ".dat";
			GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
			Dos::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, this->ef, GlobalC::wf.ekb);
		}*/
		
		if(out_dos==3)
		{
			for(int i=0; i<nspin0; i++)
			{
				std::stringstream ss3;
				ss3 << GlobalV::global_out_dir << "Fermi_Surface_" << i << ".bxsf";
				Dos::nscf_fermi_surface(ss3.str(),GlobalC::kv.nks,GlobalV::NBANDS,GlobalC::wf.ekb);
			}
		}
	}
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
			std::stringstream ss2;
			ss2 << GlobalV::global_out_dir << "BANDS_" << is+1 << ".dat";
			GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
			Dos::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, this->ef*0, GlobalC::wf.ekb);
		}

	}
	return;
}

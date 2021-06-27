#include "src_pw/tools.h"
#include "src_pw/global.h"
#include "src_pw/energy.h"
#include "src_io/mulliken_charge.h"
#include "src_lcao/LCAO_nnr.h"
#include "src_lcao/LCAO_gen_fixedH.h"    
#include "module_neighbor/sltk_atom_arrange.h"//qifeng-2019-01-21
#include "src_lcao/local_orbital_charge.h"
#include "src_pw/global.h"
#include "src_pw/wavefunc.h"
#include "src_lcao/LCAO_matrix.h"
#include "src_lcao/global_fp.h"
#include "src_lcao/wfc_dm_2d.h"
#include "module_base/lapack_connector.h"
#include "module_base/scalapack_connector.h"
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "src_io/dos.h"
#include <vector>
#include <mpi.h>
#include <sys/time.h>


void energy::perform_dos(void)
{
	TITLE("energy","perform_dos");


    if(out_dos !=0 || out_band !=0)
    {
        ofs_running << "\n\n\n\n";
        ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
        ofs_running << " |                                                                    |" << endl;
        ofs_running << " | Post-processing of data:                                           |" << endl;
        ofs_running << " | DOS (density of states) and bands will be output here.             |" << endl;
        ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << endl;
        ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << endl;
        ofs_running << " | done here.                                                         |" << endl;
        ofs_running << " |                                                                    |" << endl;
        ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        ofs_running << "\n\n\n\n";
    }


	if(MY_RANK==0)
	{
		if(CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax")
		{
			stringstream ss;
			ss << global_out_dir << "istate.info" ;
			ofstream ofsi( ss.str().c_str() );
			/*for(int ib=0; ib<NBANDS; ++ib)
			{
				ofsi << "0 " << ib+1;
				for(int is=0; is<NSPIN; ++is)
				{
					ofsi << " " << wf.ekb[is][ib];
					}
					ofsi << endl;
					}*/
			// pengfei 2015-4-1
			if(NSPIN == 1||NSPIN == 4)
			{ 
				for (int ik = 0;ik < kv.nks;ik++)
				{                       
					ofsi<<"BAND"
					<<setw(25)<<"Energy(ev)"
					<<setw(25)<<"Occupation"
					<<setw(25)<<"Kpoint = "<<ik+1
					<<setw(25)<<"("<<kv.kvec_d[ik].x<<" "<<kv.kvec_d[ik].y<<" "<<kv.kvec_d[ik].z<<")"<<endl;
					for(int ib=0;ib<NBANDS;ib++)
					{
						ofsi<<ib+1<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV<<setw(25)<<wf.wg(ik,ib)<<endl;
					}
					ofsi <<endl; 
					ofsi <<endl;                              
				}
			}
			else
			{
				for (int ik = 0;ik < kv.nks/2;ik++)
				{
					ofsi<<"BAND"<<setw(25)<<"Spin up Energy(ev)"
					<<setw(25)<<"Occupation"
					<<setw(25)<<"Spin down Energy(ev)"
					<<setw(25)<<"Occupation"
					<<setw(25)<<"Kpoint = "<<ik+1
					<<setw(25)<<"("<<kv.kvec_d[ik].x<<" "<<kv.kvec_d[ik].y<<" "<<kv.kvec_d[ik].z<<")"<<endl;
					for(int ib=0;ib<NBANDS;ib++)
					{
						ofsi<<ib+1<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV
						<<setw(25)<<wf.wg(ik,ib)
						<<setw(25)<<wf.ekb[(ik+kv.nks/2)][ib]* Ry_to_eV
						<<setw(25)<<wf.wg(ik+kv.nks/2,ib)<<endl;
					}
					ofsi <<endl;
					ofsi <<endl;

				}
			}

			ofsi.close();
		}
	}


	//qianrui modify 2020-10-18
	if(CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax")
	{
		stringstream ss;
		ss << global_out_dir << "istate.info" ;
		if(MY_RANK==0)
		{
			ofstream ofsi( ss.str().c_str() ); // clear istate.info
			ofsi.close();
		}
		for(int ip=0; ip<NPOOL; ip++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if( MY_POOL == ip )
			{
				if( RANK_IN_POOL != 0 ) continue;
				ofstream ofsi2( ss.str().c_str(), ios::app );
				if(NSPIN == 1||NSPIN == 4)
				{
					for (int ik = 0;ik < kv.nks;ik++)
					{
						ofsi2<<"BAND"
						<<setw(25)<<"Energy(ev)"
						<<setw(25)<<"Occupation"
						<<setw(25)<<"Kpoint = "<<Pkpoints.startk_pool[ip]+ik+1
						<<setw(25)<<"("<<kv.kvec_d[ik].x<<" "<<kv.kvec_d[ik].y<<" "<<kv.kvec_d[ik].z<<")"<<endl;
						for(int ib=0;ib<NBANDS;ib++)
						{
							ofsi2<<setw(6)<<ib+1<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV<<setw(25)<<wf.wg(ik,ib)<<endl;
						}
						ofsi2 <<endl;
						ofsi2 <<endl;
					}
				}
				else
				{
					for (int ik = 0;ik < kv.nks/2;ik++)
					{
						ofsi2<<"BAND"
						<<setw(25)<<"Spin up Energy(ev)"
						<<setw(25)<<"Occupation"
						<<setw(25)<<"Spin down Energy(ev)"
						<<setw(25)<<"Occupation"
						<<setw(25)<<"Kpoint = "<<Pkpoints.startk_pool[ip]+ik+1
						<<setw(25)<<"("<<kv.kvec_d[ik].x<<" "<<kv.kvec_d[ik].y<<" "<<kv.kvec_d[ik].z<<")"<<endl;

						for(int ib=0;ib<NBANDS;ib++)
						{
							ofsi2<<setw(6)<<ib+1
							<<setw(25)<<wf.ekb[ik][ib]* Ry_to_eV
							<<setw(25)<<wf.wg(ik,ib)
							<<setw(25)<<wf.ekb[(ik+kv.nks/2)][ib]* Ry_to_eV
							<<setw(25)<<wf.wg(ik+kv.nks/2,ib)<<endl;
						}
						ofsi2 <<endl;
						ofsi2 <<endl;

					}
				}

				ofsi2.close();
			}
		}
	}


	// mulliken charge analysis
	if(mulliken == 1)
	{
		Mulliken_Charge   MC;
		MC.stdout_mulliken();			
	}//qifeng add 2019/9/10

	int nspin0=1;
	if(NSPIN==2) nspin0=2;

	if(this->out_dos)
	{
		// find the maximal and minimal band energy.
		double emax = wf.ekb[0][0];
		double emin = wf.ekb[0][0];
		for(int ik=0; ik<kv.nks; ++ik)
		{
			for(int ib=0; ib<NBANDS; ++ib)
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


		//		OUT(ofs_running,"minimal energy is (eV)", emin);
		//		OUT(ofs_running,"maximal energy is (eV)", emax);
		//  output the PDOS file.////qifeng-2019-01-21
		// 		atom_arrange::set_sr_NL();
		//		atom_arrange::search( SEARCH_RADIUS );//qifeng-2019-01-21
		const double de_ev = this->dos_edelta_ev;


		const int npoints = static_cast<int>(std::floor ( ( emax - emin ) / de_ev ));

		int NUM=NLOCAL*npoints;

		Wfc_Dm_2d D;
		D.init();
		if(GAMMA_ONLY_LOCAL)
		{
			for(int in=0;in<NSPIN;in++)
			{

				D.wfc_gamma[in]=LOC.wfc_dm_2d.wfc_gamma[in];
			}

		}
		else 
		{

			for(int in=0;in<kv.nks;in++)
			{

				D.wfc_k[in] = LOC.wfc_dm_2d.wfc_k[in];
			}
		}


		const int np=npoints;
		matrix*  pdosk = new matrix[nspin0];

		for(int is=0; is<nspin0; ++is)
		{


			pdosk[is].create(NLOCAL,np,true);



		}
		matrix*  pdos = new matrix[nspin0];
		for(int is=0; is<nspin0; ++is)
		{
			pdos[is].create(NLOCAL,np,true);

		}

		double a = bcoeff;
		double c=2*3.141592653;
		double b = sqrt(c)*a;                                         

		complex<double>*waveg = new complex<double>[NLOCAL];

		double*Gauss = new double[np];

		for(int is=0; is<nspin0; ++is)
		{
			if(GAMMA_ONLY_LOCAL)
			{
				std::vector<matrix>   Mulk;
				Mulk.resize(1);
				Mulk[0].create(ParaO.ncol,ParaO.nrow);


				matrix Dwf = D.wfc_gamma[is];
				for (int i=0; i<NBANDS; ++i)		  
				{     
					ZEROS(waveg, NLOCAL);

					ZEROS(Gauss,np);
					for (int n=0; n<npoints; ++n)		  
					{  
						double en=emin+n * de_ev;
						double en0=wf.ekb[0][i]*Ry_to_eV;
						double de = en-en0;
						double de2 = 0.5*de * de;
						Gauss[n] = kv.wk[0]*exp(-de2/a/a)/b;
					}

					const int NB= i+1;

					const double one_float=1.0, zero_float=0.0;
					const int one_int=1;


					const char T_char='T';		
					pdgemv_(
							&T_char,
							&NLOCAL,&NLOCAL,
							&one_float,
							LM.Sloc, &one_int, &one_int, ParaO.desc,
							Dwf.c, &one_int, &NB, ParaO.desc, &one_int,
							&zero_float,
							Mulk[0].c, &one_int, &NB, ParaO.desc,
							&one_int);

					for (int j=0; j<NLOCAL; ++j)
					{

						if ( ParaO.in_this_processor(j,i) )
						{

							const int ir = ParaO.trace_loc_row[j];
							const int ic = ParaO.trace_loc_col[i];
							waveg[j] = Mulk[0](ic,ir)*D.wfc_gamma[is](ic,ir);
							const double x = waveg[j].real();
							LapackConnector::axpy(np , x,Gauss, 1,pdosk[is].c+j*pdosk[is].nc,1);
						}
					} 
				}//ib
			}//if
			else
			{
				SEARCH_RADIUS = atom_arrange::set_sr_NL(
					ofs_running,
					OUT_LEVEL,
					ORB.get_rcutmax_Phi(), 
					ORB.get_rcutmax_Beta(), 
					GAMMA_ONLY_LOCAL);

				atom_arrange::search(
					SEARCH_PBC,
					ofs_running,
					GridD, 
					ucell, 
					SEARCH_RADIUS, 
					test_atom_input);//qifeng-2019-01-21

				// mohan update 2021-04-16
				LOWF.orb_con.set_orb_tables(
						ofs_running,
						UOT, 
						ORB,
						ucell.ntype,
						ucell.lmax,
						INPUT.lcao_ecut,
						INPUT.lcao_dk,
						INPUT.lcao_dr,
						INPUT.lcao_rmax, 
						ucell.lat0, 
						INPUT.out_descriptor,
						INPUT.out_r_matrix,
						Exx_Abfs::Lmax,
						FORCE,
						MY_RANK);

				LM.allocate_HS_R(LNNR.nnr);
				LM.zeros_HSR('S', LNNR.nnr);
				UHM.genH.calculate_S_no();
				UHM.genH.build_ST_new('S', false);
				std::vector<ComplexMatrix> Mulk;
				Mulk.resize(1);
				Mulk[0].create(ParaO.ncol,ParaO.nrow);


				for(int ik=0;ik<kv.nks;ik++)
				{

					if(is == kv.isk[ik])
					{
						LM.allocate_HS_k(ParaO.nloc);
						LM.zeros_HSk('S');
						LNNR.folding_fixedH(ik);


						ComplexMatrix Dwfc = conj(D.wfc_k[ik]);

						for (int i=0; i<NBANDS; ++i)		  
						{     

							ZEROS(waveg, NLOCAL);


							ZEROS(Gauss,np);
							for (int n=0; n<npoints; ++n)		  
							{  
								double en=emin+n * de_ev;
								double en0=wf.ekb[ik][i]*Ry_to_eV;
								double de = en-en0;
								double de2 = 0.5*de * de;
								Gauss[n] = kv.wk[ik]*exp(-de2/a/a)/b;
							}

							const int NB= i+1;

							const double one_float=1.0, zero_float=0.0;
							const int one_int=1;
							//   const int two_int=2;
							const char T_char='T';		// N_char='N',U_char='U'

							pzgemv_(
									&T_char,
									&NLOCAL,&NLOCAL,
									&one_float,
									LM.Sloc2, &one_int, &one_int, ParaO.desc,
									Dwfc.c, &one_int, &NB, ParaO.desc, &one_int,
									&zero_float,
									Mulk[0].c, &one_int, &NB, ParaO.desc,
									&one_int);



							for (int j=0; j<NLOCAL; ++j)
							{

								if ( ParaO.in_this_processor(j,i) )
								{

									const int ir = ParaO.trace_loc_row[j];
									const int ic = ParaO.trace_loc_col[i];

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
					ofs_running,
					SEARCH_PBC, 
					GridD, 
					ucell, 
					SEARCH_RADIUS, 
					test_atom_input);
#endif
				// mohan update 2021-02-10
				LOWF.orb_con.clear_after_ions(UOT, ORB, INPUT.out_descriptor);
			}//else

		 MPI_Reduce(pdosk[is].c, pdos[is].c , NUM , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);
	 }//is                                              
	 delete[] pdosk;                                               
	 delete[] waveg;
	 if(MY_RANK == 0)
	 {
		 {  stringstream ps;
			 ps << global_out_dir << "TDOS";
			 ofstream out(ps.str().c_str());
			 if (NSPIN==1)
			 {

				 for (int n=0; n<npoints; ++n)
				 { double y=0.0;
					 double en=emin + n * de_ev;
					 for (int i=0; i<NLOCAL; i++)
					 {
						 y +=  pdos[0](i,n);
					 }  

					 out <<setw(20)<< en <<setw(30)<< y << endl;
				 }
			 }
			 else if (NSPIN==2)
			 {
				 for (int n=0; n<npoints; ++n)
				 { double y=0.0;
					 double z=0.0;
					 double en=emin + n * de_ev;
					 for (int i=0; i<NLOCAL; i++)
					 {
						 y +=  pdos[0](i,n);
						 z +=  pdos[1](i,n);

					 }  

					 out <<setw(20)<< en <<setw(30)<< y << setw(30)<< z<< endl;
				 }
			 }
			 out.close();
		 }

		 string Name_Angular[5][11];
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

		 {stringstream as;
			 as << global_out_dir << "PDOS";
			 ofstream out(as.str().c_str());

			 out << "<"<<"pdos"<<">" <<endl;
			 out << "<"<<"nspin"<<">" << NSPIN<< "<"<<"/"<<"nspin"<<">"<< endl;
			 out << "<"<<"norbitals"<<">" <<setw(2) <<NLOCAL<< "<"<<"/"<<"norbitals"<<">"<< endl;
			 out << "<"<<"energy"<<"_"<<"values units"<<"="<<"\""<<"eV"<<"\""<<">"<<endl;

			 for (int n=0; n<npoints; ++n)
			 { double y=0.0;
				 double en=emin + n * de_ev;
				 out <<setw(20)<< en << endl;
			 }
			 out << "<"<<"/"<<"energy"<<"_"<<"values"<<">" <<endl;
			 for (int i=0; i<ucell.nat; i++)
			 {   
				 int a = ucell.iat2ia[i];
				 int t = ucell.iat2it[i];
				 Atom* atom1 = &ucell.atoms[t];
				 for(int j=0; j<atom1->nw; ++j)
				 {
					 const int L1 = atom1->iw2l[j];
					 const int N1 = atom1->iw2n[j];
					 const int m1 = atom1->iw2m[j];
					 const int w = ucell.itiaiw2iwt(t, a, j);

					 //out << "<"<<"/"<<"energy"<<"_"<<"values"<<">" <<endl;
					 out << "<"<<"orbital" <<endl;
					 out <<setw(6)<< "index"<<"="<<"\""<<setw(40) <<w+1<<"\""<<endl;
					 out <<setw(5)<< "atom"<<"_"<<"index"<<"="<<"\""<<setw(40) <<i+1<<"\""<<endl;
					 out <<setw(8)<< "species"<<"="<<"\""<<ucell.atoms[t].label<<"\""<<endl;
					 out<<setw(2)<< "l"<<"="<<"\""<<setw(40)<<L1<<"\""<<endl;
					 out <<setw(2)<< "m"<<"="<<"\""<<setw(40)<<m1<<"\""<<endl;
					 out <<setw(2)<< "z"<<"="<<"\""<<setw(40)<<N1+1<<"\""<<endl;
					 out << ">" <<endl;
					 out << "<"<<"data"<<">" <<endl;
					 if (NSPIN==1)
					 {
						 for (int n=0; n<npoints; ++n)
						 {


							 out <<setw(13)<< pdos[0](w,n)<<endl;
						 }//n
					 }
					 else if (NSPIN==2)
					 {
						 for (int n=0; n<npoints; ++n)
						 {
							 out <<setw(20)<< pdos[0](w,n)<< setw(30)<< pdos[1](w,n)<<endl;
						 }//n
					 }

					 out << "<"<<"/"<<"data"<<">" <<endl;

				 }//j
			 }//i
			 out << "<"<<"/"<<"orbital"<<">" <<endl;
			 out << "<"<<"/"<<"pdos"<<">" <<endl;
			 out.close();}
		 {  stringstream os;
			 os<<global_out_dir<<"Orbital";
			 ofstream out(os.str().c_str());
			 out<< setw(5)<<"io"<< setw(8) <<"spec" <<setw(5)<<"l"<<setw(5)<<"m"<<setw(5)<<"z"<<setw(5)<<"sym"<<endl;


			 for (int i=0; i<ucell.nat; i++)
			 {
				 int   t = ucell.iat2it[i];
				 Atom* atom1 = &ucell.atoms[t];  
				 for(int j=0; j<atom1->nw; ++j)
				 {
					 const int L1 = atom1->iw2l[j];
					 const int N1 = atom1->iw2n[j];
					 const int m1 = atom1->iw2m[j];
					 out <<setw(5) << i << setw(8) 
						<< ucell.atoms[t].label <<setw(5)
							<<L1<<setw(5) <<m1<<setw(5)<<N1+1<<setw(15)<< Name_Angular[L1][m1] << endl;
				 }
			 }
			 out <<endl<<endl;
			 out <<setw(5)<< "io"<<setw(2)<<"="<<setw(2)<<"Orbital index in supercell"<<endl;
			 out <<setw(5)<< "spec"<<setw(2)<<"="<<setw(2)<<"Atomic species label"<<endl;
			 out <<setw(5)<< "l"<<setw(2)<<"="<<setw(2)<<"Angular mumentum quantum number"<<endl;
			 out <<setw(5)<< "m"<<setw(2)<<"="<<setw(2)<<"Magnetic quantum number"<<endl;
			 out <<setw(5)<< "z"<<setw(2)<<"="<<setw(2)<<"Zeta index of orbital"<<endl;
			 out <<setw(5)<< "sym"<<setw(2)<<"="<<setw(2)<<"Symmetry name of real orbital"<<endl;
			 out.close();}

	 }       
	 delete[] pdos;

	 // output the DOS file.
	 for(int is=0; is<nspin0; ++is)
	 {
		 stringstream ss;
		 ss << global_out_dir << "DOS" << is+1;

		 Dos::calculate_dos(
				 is,
				 kv.isk,
				 ss.str(), 
				 this->dos_edelta_ev, 
				 emax, 
				 emin, 
				 kv.nks, kv.nkstot, kv.wk, wf.wg, NBANDS, wf.ekb );
		 ifstream in(ss.str().c_str());
		 if(!in)
		 {
			 //      cout<<"\n Can't find file : "<< name << endl;
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
			 //cout<<"\n Read Over!"<<endl;
		 }
		 in.close();

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
		 sss << global_out_dir << "DOS" << is+1 << "_smearing" << ".dat" ;
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


		// mulliken charge analysis
		if(out_dos == 2)
		{
			stringstream sp;
			sp << global_out_dir << "Mulliken.dat";
			Dos::calculate_Mulliken(sp.str());
		}
	
		if(nspin0==1)
		{
			ofs_running << " Fermi energy is " << this->ef << " Rydberg" << endl;
		}
		else if(nspin0==2)
		{
			ofs_running << " Fermi energy (spin = 1) is " << this->ef_up << " Rydberg" << endl;
			ofs_running << " Fermi energy (spin = 2) is " << this->ef_dw << " Rydberg" << endl;
		}

		//int nks;
		//if(nspin0==1) nks = kv.nkstot;
		//else if(nspin0==2) nks = kv.nkstot/2;



		/*for(int is=0; is<NSPIN; is++)
		{
			stringstream ss2;
			ss2 << global_out_dir << "BANDS_" << is+1 << ".dat";
			ofs_running << "\n Output bands in file: " << ss2.str() << endl;
			Dos::nscf_band(is, ss2.str(), nks, NBANDS, this->ef, wf.ekb);
		}*/
		
		if(out_dos==3)
		{
			for(int i=0; i<nspin0; i++)
			{
				stringstream ss3;
				ss3 << global_out_dir << "Fermi_Surface_" << i << ".bxsf";
				Dos::nscf_fermi_surface(ss3.str(),kv.nks,NBANDS,wf.ekb);
			}
		}
	}
	if(this->out_band) //pengfei 2014-10-13
	{
		int nks=0;
		if(nspin0==1) 
		{
			nks = kv.nkstot;
		}
		else if(nspin0==2) 
		{
			nks = kv.nkstot/2;
		}

		for(int is=0; is<nspin0; is++)
		{
			stringstream ss2;
			ss2 << global_out_dir << "BANDS_" << is+1 << ".dat";
			ofs_running << "\n Output bands in file: " << ss2.str() << endl;
			Dos::nscf_band(is, ss2.str(), nks, NBANDS, this->ef*0, wf.ekb);
		}

	}
	return;
}

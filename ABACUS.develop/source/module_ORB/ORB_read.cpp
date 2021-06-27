#include "ORB_read.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02
#include <cassert>
#include "../src_global/math_integral.h"
#include "../src_global/tool_check.h"
#include <algorithm>
//#include "src_pw/soc.h"
using namespace std;

//==============================
// Define an object here! 
//==============================
/// PLEASE avoid using 'ORB' as global variable 
// mohan note 2021-03-23
LCAO_Orbitals ORB;

LCAO_Orbitals::LCAO_Orbitals()
{
	this->nchimax = 0;// this initialzied must specified
	this->Phi = new Numerical_Orbital[1];	
	this->Beta = new Numerical_Nonlocal[1];
	this->Alpha = new Numerical_Orbital[1];

	this->nproj = new int[1];
	this->nprojmax = 0;
	this->read_in_flag = false;	

	this->dr_uniform = 0.001;

}

LCAO_Orbitals::~LCAO_Orbitals()
{
	delete[] Phi;
	delete[] Beta;
	delete[] Alpha;
	delete[] nproj;
}

#ifdef __MPI
// be called in unitcell_pseudo.
void LCAO_Orbitals::bcast_files(
	const int &ntype_in, 
	const int &my_rank)
{
	TITLE("LCAO_Orbitals","bcast_files");

	// 'read_in_flag' is true when there is a
	// block "NUMERICAL_ORBITAL" in structure
	// file.
	Parallel_Common::bcast_bool(read_in_flag);
	if(!read_in_flag)
	{
		return;
	}

	assert(ntype_in > 0 );

	ofs_running << "\n READING ORBITAL FILE NAMES FOR LCAO" << endl;
	for(int it=0; it<ntype_in; it++)
	{
		string ofile;
		string nfile;

		if(my_rank==0)
		{
			ofile = orbital_file[it];
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
//			nfile = nonlocal_file[it];
		}

// PLEASE avoid using 'bcast_string' as global variable 
// mohan note 2021-03-23
		Parallel_Common::bcast_string(ofile);
		//-----------------------------------
		// Turn off the read in NONLOCAL file
		// function since 2013-08-02 by mohan
		//-----------------------------------
//		Parallel_Common::bcast_string(nfile);

		if(my_rank!=0)
		{
			orbital_file.push_back( ofile );
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
//			nonlocal_file.push_back ( nfile );
		}

		ofs_running << " orbital file: " << orbital_file[it] << endl;
//		ofs_running << " nonlocal file: " << nonlocal_file[it] << endl;
	}
	return;
}
#endif


void LCAO_Orbitals::Read_Orbitals(
	ofstream &ofs_in,
	const int &ntype_in, 
	const int &lmax_in,
	const int &out_descriptor,
	const int &out_r_matrix,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank) // mohan add 2021-04-26
{
	TITLE("LCAO_Orbitals", "Read_Orbitals");
	timer::tick("LCAO_Orbitals","Read_Orbitals",'C');

	ofs_in << "\n\n\n\n";
	ofs_in << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_in << " |                                                                    |" << endl;
	ofs_in << " | Setup numerical orbitals:                                          |" << endl;
	ofs_in << " | This part setup: numerical atomic orbitals, non-local projectors   |" << endl;
	ofs_in << " | and neutral potential (1D). The atomic orbitals information        |" << endl;
	ofs_in << " | including the radius, angular momentum and zeta number.            |" << endl;
	ofs_in << " | The neutral potential is the sum of local part of pseudopotential  |" << endl;
	ofs_in << " | and potential given by atomic charge, they will cancel out beyond  |" << endl;
	ofs_in << " | a certain radius cutoff, because the Z/r character.                |" << endl;
	ofs_in << " |                                                                    |" << endl;
	ofs_in << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_in << "\n\n\n\n";	

	//////////////////////
	/// (1) check dk, dR, Rmax.
	//////////////////////

	ofs_in << "\n SETUP ONE DIMENSIONAL ORBITALS/POTENTIAL" << endl;

	if(!read_in_flag)
	{
		WARNING_QUIT("LCAO_Orbitals::Read_Orbitals","Set the NUMERICAL_ORBITAL block in structure file.");
	}


	//OUT(ofs_in,"ecutwfc for kmesh",ecutwfc);
	OUT(ofs_in,"delta k  (1/Bohr)",dk);
	OUT(ofs_in,"delta r    (Bohr)",dR);
	OUT(ofs_in,"dr_uniform (Bohr)",dr_uniform);
	OUT(ofs_in,"rmax       (Bohr)",Rmax);

	// check the read in data.
    assert(dk > 0.0);
    assert(ecutwfc > 0.0);
    assert(dR > 0.0);
    assert(Rmax > 0.0);

	/// ntype: number of atom species
	this->ntype = ntype_in; 
	assert(ntype>0);

	/// lmax: lmax used in local orbitals as basis sets
	assert(lmax_in>=0); // mohan add 2021-04-16
	this->lmax = lmax_in;

	//////////////////////////////////////////////////////////
	/// (2) set the kmesh according to ecutwfc and dk. 
	//////////////////////////////////////////////////////////

	//-----------------------------------------------------------------
	/// calculate number of k mesh according to energy cutoff.
	/// Mohan choose ecutwfc according to interpolation requirement.
	//	cout << " ecutwfc=" << ecutwfc << endl;
	//LiuXh modified 2016-01-25, 2016-07-20
	if(ecutwfc< 20)
	{
		this->kmesh = static_cast<int>( 2 * sqrt(ecutwfc) / dk )  + 4;
	}
	else
	{
		this->kmesh = static_cast<int>( sqrt(ecutwfc) / dk )  + 4;
	}

	// jingan add for calculate r(R) matrix
	if(out_r_matrix) 
	{
		kmesh = kmesh * 4;
	}

	//	this->kmesh = static_cast<int> (PI / 0.01 / 4 / this->dk);
	if(kmesh%2==0) kmesh++;
	OUT(ofs_in,"kmesh",kmesh);
	//-----------------------------------------------------------------


	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//~~~~~~~~~~~~~~~~~~~~~~   1    ~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read in numerical atomic orbitals for each atom type.
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	delete[] this->Phi;

	this->Phi = new Numerical_Orbital[ntype];
	for(int it=0; it<ntype; it++)
	{
		this->Read_PAO(ofs_in, it, force_flag, my_rank);
		//caoyu add 2021-05-24	to reconstruct atom_arrange::set_sr_NL
		this->rcutmax_Phi = std::max(this->rcutmax_Phi, this->Phi[it].getRcut());
	}

	
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	//~~~~~~~~~~~~~~~~~~~~~~   2    ~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read in non-local projector for each atom type.
	// In fact this should be improved,
	// the non-local projector should be transferred
	// from .UPF file directly.
	// mohan note 2011-03-04
	//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	delete[] this->Beta;
	this->Beta = new Numerical_Nonlocal[ntype];

	delete[] nproj;
	this->nproj = new int[ntype];
	ZEROS(nproj, ntype);
	
	this->nprojmax = 0;
	

	// if true: read in the nonlocal file from file.
	// if false: get nonlocal information from .upf or .vwr directly
	bool readin_nonlocal = false;

	for(int it=0; it<ntype; it++)
	{
#ifdef __NORMAL
		// need to be added 2021-04-26 mohan
#else
		if(readin_nonlocal)
		{
			this->Read_NonLocal(it, this->nproj[it], my_rank);	
		}
		else
		{
			this->Set_NonLocal(it, this->nproj[it]);
		}
#endif
		this->nprojmax = std::max(this->nprojmax, this->nproj[it]);
		//caoyu add 2021-05-24 to reconstruct atom_arrange::set_sr_NL
		this->rcutmax_Beta = std::max(this->rcutmax_Beta, this->Beta[it].get_rcut_max());
	}

	ofs_in << " max number of nonlocal projetors among all species is " << nprojmax << endl; 

	//caoyu add 2021-3-16
	///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	///~~~~~~~~~~~~~~~~~~~~~~   3    ~~~~~~~~~~~~~~~~~~~~~~~~~
	/// Read in numerical basis for descriptor.
	///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


	if (out_descriptor>0)	//condition: descriptor in lcao line
	{
		
		delete[] this->Alpha;
		this->Alpha = new Numerical_Orbital[1];	//not related to atom type -- remain to be discussed

		this->Read_Descriptor(ofs_in, force_flag, my_rank);

	}

	timer::tick("LCAO_Orbitals","Read_Orbitals",'C');
	return;
}


#ifdef __NORMAL

#else

#include "../src_pw/global.h"
// mohan add 2013-08-02
// In order to get rid of the read in file .NONLOCAL.
void LCAO_Orbitals::Set_NonLocal(const int &it, int &n_projectors)
{
	TITLE("LCAO_Orbitals","Set_NonLocal");

	// set a pointer
	Atom* atom = &ucell.atoms[it];

	// get the number of non-local projectors
	n_projectors = atom->nbeta;

	const int nh = atom->nh;//zhengdy-soc

	// set the nonlocal projector objects
	Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];

	ComplexMatrix coefficient_D_nc_in(nh*2, nh*2);//zhengdy-soc

	if(!atom->has_so)
	{
		for(int p1 = 0; p1<n_projectors; p1++)
		{
			const int lnow = atom->lll[p1];

			// only keep the nonzero part.
			int cut_mesh = atom->mesh; 
			for(int ir=atom->mesh-1; ir>=0; --ir)
			{
				if( abs( atom->betar(p1,ir) ) > 1.0e-10 )
				{
					cut_mesh = ir; 
					break;
				}
			}
			if(cut_mesh %2 == 0) ++cut_mesh;

//		cout << " cut_mesh=" << cut_mesh << endl;
			double* beta_r = new double[cut_mesh];
			ZEROS(beta_r, cut_mesh);
			for(int ir=0; ir<cut_mesh; ++ir)
			{
				beta_r[ir] = atom->betar(p1,ir);
			}

			tmpBeta_lm[p1].set_NL_proj(
					atom->label,
					it, //type
					lnow, // angular momentum L
					cut_mesh, // number of radial mesh
					atom->rab,
					atom->r, // radial mesh value (a.u.)
					beta_r,
					this->kmesh,
					this->dk,
					dr_uniform); // delta k mesh in reciprocal space

			tmpBeta_lm[p1].plot(MY_RANK);

			delete[] beta_r;
				
		}
		
		// mohan comment out 2021-04-26
		//WARNING("LCAO_Orbitals::Set_NonLocal","bug in line "+TO_STRING(__LINE__)+", matrix ic>=nc");		


		// Peize Lin add 2019-01-23
		this->Beta[it].set_type_info(
			it, 
			atom->label, 
			atom->pp_type, 
			atom->lmax, 
			n_projectors, 
			tmpBeta_lm);//LiuXh 2016-01-14, 2016-07-19

		// mohan add 2021-05-07
		atom->set_d_so(coefficient_D_nc_in,n_projectors,0,0);
	}
	else//added by zhengdy-soc
	{
		int lmaxkb = - 1;
		for (int ibeta = 0; ibeta < atom->nbeta; ibeta++)
		{
			lmaxkb = max( lmaxkb, atom->lll[ibeta]);
		}
		Soc soc;
		soc.rot_ylm(lmaxkb);
		soc.fcoef.create(ucell.ntype, atom->nh, atom->nh);

		int ip1=0;
		for(int p1 = 0; p1<n_projectors; p1++)//nbeta
		{
			const int lnow = atom->lll[p1];
			
			const int l1 = atom->lll[p1];
			const double j1 = atom->jjj[p1];
			for(int m1=0; m1<2*l1+1; m1++)
			{
				int ip2=0;
				for(int p2 = 0; p2<n_projectors; p2++)
				{
					const int l2 = atom->lll[p2];
					const double j2 = atom->jjj[p2];
					for(int m2 = 0;m2<2*l2+1;m2++)
					{
						if(l1==l2 && fabs(j1-j2)<1e-7)
						{
							for(int is1=0;is1<2;is1++)
							{
								for(int is2=0;is2<2;is2++)
								{
									soc.set_fcoef(l1, l2,
										is1, is2,
										m1, m2,
										j1, j2,
										it, ip1, ip2);
									coefficient_D_nc_in(ip1 + nh*is1, ip2 + nh*is2) = atom->dion(p1,p2) 
									* soc.fcoef(it, is1, is2, ip1, ip2);
									if(p1 != p2) 
									{
										soc.fcoef(it, is1, is2, ip1, ip2) = complex<double>(0.0,0.0);
									}
								}// end is2
							}// end is1
						}// end l1==l2
						ip2++;
					}// end m2
				}// end p2
				assert(ip2==nh);
				ip1++;
			}// end m1

			// only keep the nonzero part.
			int cut_mesh = atom->mesh; 
			for(int ir=atom->mesh-1; ir>=0; --ir)
			{
				if( abs( atom->betar(p1,ir) ) > 1.0e-10 )
				{
					cut_mesh = ir; 
					break;
				}
			}
			if(cut_mesh %2 == 0) 
			{
				++cut_mesh;
			}

			double* beta_r = new double[cut_mesh];
			ZEROS(beta_r, cut_mesh);
			for(int ir=0; ir<cut_mesh; ++ir)
			{
				beta_r[ir] = atom->betar(p1,ir);
			}

			tmpBeta_lm[p1].set_NL_proj(
					atom->label,
					it, //type
					lnow, // angular momentum L
					cut_mesh, // number of radial mesh
					atom->rab,
					atom->r, // radial mesh value (a.u.)
					beta_r,
					this->kmesh,
					this->dk,
					dr_uniform); // delta k mesh in reciprocal space

			tmpBeta_lm[p1].plot(MY_RANK);

			delete[] beta_r;
		}

		assert(ip1==nh);

		this->Beta[it].set_type_info(
			it, 
			atom->label, 
			atom->pp_type, 
			atom->lmax, 
			n_projectors, 
			tmpBeta_lm);//zhengdy-soc 2018-09-10

		// mohan add 2021-05-07
		atom->set_d_so(coefficient_D_nc_in,n_projectors,nh,1);

	}//end if

	delete[] tmpBeta_lm;

	cout << " SET NONLOCAL PSEUDOPOTENTIAL PROJECTORS" << endl;
	return;
}


void LCAO_Orbitals::Read_NonLocal(
	const int &it, 
	int &n_projectors,
	const int &my_rank)
{
	TITLE("LCAO_Orbitals","Read_NonLocal");

	ifstream ifs;

	// mohan add 2010-09-08.
	// check if the non-local pseudopotential file exist.
	bool open = false;	
	if(my_rank==0)
	{
		ifs.open( this->nonlocal_file[it].c_str() );
		if(ifs)
		{
			open = true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(open);
#endif
	if(!open)
	{
		cout << " Non-local File : " << nonlocal_file[it] << endl;
		WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","Can not find the NONLOCAL file.");
	}
	else
	{
//		ofs_running << " open nonLocal pseudopotential file: " << nonlocal_file[it] << endl;
	}


	string label;
	string ps_type;

	// maximal lmax allowed in this calculation
	int nlmax = 0;

	if(my_rank==0)
	{
		if(SCAN_BEGIN(ifs, "<HEADER>"))
		{
			READ_VALUE(ifs, label);
			READ_VALUE(ifs, ps_type);
			if(ps_type != "NC")
			{
				WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","Only available for NC nonlocal pseudopotential");
			}
			READ_VALUE(ifs, nlmax);
//			cout << " " << label << " " << ps_type << " " << nlmax << endl; 
			assert(nlmax >= -1);
			SCAN_END(ifs,"</HEADER>");
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_string(label);
	Parallel_Common::bcast_string(ps_type);
	Parallel_Common::bcast_int(nlmax);
#endif

	//mohan add 2012-06-09
	if( nlmax != -1 )
	{
		bool find_lmax = false;
		for(int ic=0; ic<ucell.atoms[it].nbeta; ic++)
		{
			if( nlmax == ucell.atoms[it].lll[ic] )
			{
				//			cout << " nlmax = " << nlmax << endl;
				//			cout << " lchi = " << ucell.atoms[it].lll[ic] << endl;
				find_lmax = true;
				break;
			}
		}

		if( !find_lmax )
		{
			cout << " For element " << label << endl;
			cout << " Max L Read in from NONLOCAL = " << nlmax << endl;
			for(int ib=0; ib<ucell.atoms[it].nbeta; ++ib)
			{
				cout << " Max L Read in from pseudopotential file = " << ucell.atoms[it].lll[ib] << endl;
			}
			WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","nlmax != ucell.atoms[it].lll");
		}
	}


//	OUT(ofs_running,"Type",it);
	OUT(ofs_running,"label",label);
//	OUT(ofs_running,"ps_type",ps_type);
	OUT(ofs_running,"nlmax",nlmax);

	//-------------------------------------------
	// if each L has projectors more than once,
	// this needed to be modified.	
	//-------------------------------------------
	int nproj_allowed = nlmax+1;
	matrix coefficient_D_in(nproj_allowed, nproj_allowed);
	ComplexMatrix coefficient_D_nc_in(nproj_allowed*2, nproj_allowed*2);

//	OUT(ofs_running,"nproj_allowed",nproj_allowed);

	if(my_rank==0)
	{
		if(SCAN_BEGIN(ifs, "<DIJ>"))
		{
			//--------------------------------------
			// this parameter is very important!!!
			//--------------------------------------
			READ_VALUE(ifs, n_projectors);
			OUT(ofs_running,"n_projectors",n_projectors);
			
			for (int p1 = 0; p1 < n_projectors; p1++)
        	{
            	for (int p2 = 0; p2 < n_projectors; p2++)
            	{
                	int L1_read, L2_read;
                	
					ifs >> L1_read >> L2_read;
					
					assert(L1_read <= nlmax);
					assert(L2_read <= nlmax);
                	
					ifs >> coefficient_D_in(L1_read, L2_read);
					
//					ofs_running << " L1=" << L1_read << " L2=" << L2_read << " Coef=" << coefficient_D_in(L1_read,L2_read) << endl;
            	}
        	}
			SCAN_END(ifs,"</DIJ>");
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_int(n_projectors); // mohan add 2010-12-20
//	Parallel_Common::bcast_double(coefficient_D_in.c, coefficient_D_in.nr * coefficient_D_in.nc);
#endif

	Numerical_Nonlocal_Lm* tmpBeta_lm = new Numerical_Nonlocal_Lm[n_projectors];
	int* LfromBeta = new int[n_projectors];
	ZEROS(LfromBeta, n_projectors);

	for(int p1 = 0; p1<n_projectors; p1++)
	{
		int meshr_ps = 0;
		if(my_rank==0)
		{
			if(SCAN_BEGIN(ifs, "<PP_BETA>", 0))
			{
				int iproj;
				READ_VALUE(ifs, iproj);
				if(iproj!=p1)
				{
					cout << " iproj=" << iproj << " p1=" << p1 << endl;
					WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","Check non-local projector index.");
				}
				
				READ_VALUE(ifs, LfromBeta[p1]);
				assert( LfromBeta[p1] >= 0 );
				assert( LfromBeta[p1] <= nlmax );

				READ_VALUE(ifs, meshr_ps);
				if(meshr_ps%2==0)
				{
					cout << " meshr_ps = " << meshr_ps << endl;
					WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","meshr_ps must be odd!");
				}
			}
			else
			{
				WARNING_QUIT("LCAO_Orbitals::Read_NonLocal","<PP_BETA> doesn't match!");
			}
		}// end my_rank==0

//		OUT(ofs_running,"meshr_ps",meshr_ps);

#ifdef __MPI
		Parallel_Common::bcast_int(meshr_ps);
		Parallel_Common::bcast_int(LfromBeta[p1]);
#endif		

		double* radial_ps = new double[meshr_ps];
		double* rab_ps = new double[meshr_ps];
		double* beta_r = new double[meshr_ps];
		ZEROS(radial_ps, meshr_ps);
		ZEROS(rab_ps, meshr_ps);
		ZEROS(beta_r, meshr_ps);

		if(my_rank==0)
		{
            for (int ir = 0; ir < meshr_ps; ir++)
            {
               	ifs >> radial_ps[ir];
               	ifs >> beta_r[ir];
               	ifs >> rab_ps[ir];
            }
		}

#ifdef __MPI
		Parallel_Common::bcast_double(radial_ps, meshr_ps);
		Parallel_Common::bcast_double(beta_r, meshr_ps);
		Parallel_Common::bcast_double(rab_ps, meshr_ps);
#endif
			
//		OUT(ofs_running,"radial_ps max",radial_ps[meshr_ps-1]);

//		cout << this->kmesh << endl;
        tmpBeta_lm[p1].set_NL_proj(
        		label,
                it, //type
                LfromBeta[p1], //angular momentum L
                meshr_ps, // number of radial mesh
                rab_ps,
                radial_ps,// radial mesh value(a.u.)
                beta_r,
                this->kmesh,
                this->dk,
				dr_uniform); // delta k mesh in reciprocal space

		tmpBeta_lm[p1].plot(my_rank);

		delete[] radial_ps;
		delete[] rab_ps;
		delete[] beta_r;
		
		if(my_rank==0)
		{
			SCAN_END(ifs,"</PP_BETA>");
		}
	}// end projectors.
	
	this->Beta[it].set_type_info(
		it, 
		label, 
		ps_type, 
		nlmax, 
		n_projectors, 
		tmpBeta_lm);
		
	ifs.close();

	delete[] LfromBeta;
	delete[] tmpBeta_lm;

	return;
}

#endif // end for Read and Set Nonlocal, mohan add 2021-04-26


//-------------------------------------------------------
// mohan note 2021-04-26
// to_caoyu: 
// 1. read in lmaxt and nchi directly from orbital files
// 2. pass nchi to phi via this->Phi[it].set_orbital_info 
// be careful! nchi[l] may be different for differnt phi
//-------------------------------------------------------
void LCAO_Orbitals::Read_PAO(
	ofstream &ofs_in,
	const int& it, 
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank) // mohan add 2021-04-26
{
	TITLE("LCAO_Orbitals","Read_PAO");

	ifstream in_ao;
	bool open=false;
	if(my_rank==0)
	{
		in_ao.open(this->orbital_file[it].c_str());
		if(in_ao)
		{
			open=true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool( open );
#endif
	if(!open)
	{
		cout << " Orbital file : " << this->orbital_file[it] << endl;
		WARNING_QUIT("LCAO_Orbitals::Read_PAO","Couldn't find orbital files");
	}

	ofs_in << " " << setw(12) << "ORBITAL" << setw(3) << "L" 
	<< setw(3) << "N" << setw(8) << "nr" << setw(8) << "dr"
	<< setw(8) << "RCUT" << setw(12) << "CHECK_UNIT"
		<< setw(12) << "NEW_UNIT" << endl;
	
	//lmax and nchimax for type it
	int lmaxt=0;
	int nchimaxt=0;

	this->read_orb_file(ofs_in, in_ao, it, lmaxt, nchimaxt, this->Phi, force_flag, my_rank);

	//lmax and nchimax for all types
	this->lmax = std::max(this->lmax, lmaxt);
	this->nchimax = std::max(this->nchimax, nchimaxt);
	
	in_ao.close();
	return;
}


//caoyu add 2021-3-16
void LCAO_Orbitals::Read_Descriptor(
	ofstream &ofs_in,
	const bool &force_flag, // mohan add 2021-05-07
	const int &my_rank)	//read descriptor basis
{
	TITLE("LCAO_Orbitals", "Read_Descriptor");

	ifstream in_de;
	ofs_in << " " << setw(12) << "DESCRIPTOR" << setw(3) << "L"
		<< setw(3) << "N" << setw(8) << "nr" << setw(8) << "dr"
		<< setw(8) << "RCUT" << setw(12) << "CHECK_UNIT"
		<< setw(12) << "NEW_UNIT" << endl;

	// check if the descriptor file exists.
	bool open = false;
	if (my_rank == 0)
	{
		in_de.open(this->descriptor_file.c_str());
		if (in_de)
		{
			open = true;
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(open);
#endif
	if (!open)
	{
		cout << " Orbital file : " << this->descriptor_file << endl;
		WARNING_QUIT("LCAO_Orbitals::Read_Descriptor", "Couldn't find orbital files for descriptor");
	}

	this->lmax_d = 0;
	this->nchimax_d = 0;

	this->read_orb_file(ofs_in, in_de, 0, this->lmax_d, this->nchimax_d, this->Alpha, force_flag, my_rank);

	in_de.close();

	return;
}


void LCAO_Orbitals::read_orb_file(
	ofstream &ofs_in, // ofs_running
	ifstream &ifs,
	const int &it, 
	int &lmax, 
	int &nchimax, 
	Numerical_Orbital* ao,
	const bool &force_flag,
	const int &my_rank)
{
	TITLE("LCAO_Orbitals","read_orb_file");
	char word[80];
	string orb_label;
	if (my_rank == 0)
	{
		while (ifs.good())
		{
			ifs >> word;
			if (std::strcmp(word, "Element") == 0)
			{
				ifs >> orb_label;
				continue;
			}
			if (std::strcmp(word, "Lmax") == 0)
			{
				ifs >> lmax;
				break;
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_int(lmax);
#endif	
	
	int* nchi = new int[lmax];		// allocate space: number of chi for each L.
	
	if (my_rank == 0)
	{	
		for (int l = 0; l <= lmax; l++)
		{
			ifs >> word >> word >> word >> nchi[l];
			nchimax = std::max(nchimax, nchi[l]);
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_int(nchimax);
	Parallel_Common::bcast_int(nchi, lmax + 1);
#endif		

	// calculate total number of chi
	int total_nchi = 0;
	for (int l = 0; l <= lmax; l++)
	{
		total_nchi += nchi[l];
	}

	//OUT(ofs_running,"Total number of chi(l,n)",total_nchi);
	delete[] ao[it].phiLN;
	ao[it].phiLN = new Numerical_Orbital_Lm[total_nchi];

	int meshr=0; // number of mesh points
	int meshr_read=0;
	double dr=0.0; 

	if (my_rank == 0)
	{
		while (ifs.good())
		{
			ifs >> word;
			if (std::strcmp(word, "END") == 0)		// Peize Lin fix bug about strcmp 2016-08-02
			{
				break;
			}
		}
		CHECK_NAME(ifs, "Mesh");
		ifs >> meshr;
		meshr_read = meshr;
		if (meshr % 2 == 0)
		{
			++meshr;
		}
		CHECK_NAME(ifs, "dr");
		ifs >> dr;
	}

#ifdef __MPI
	Parallel_Common::bcast_int(meshr);
	Parallel_Common::bcast_int(meshr_read);
	Parallel_Common::bcast_double(dr);
#endif		

	int count = 0;
	string name1;
	string name2;
	string name3;
	int tmp_it=0;
	int tmp_l=0;
	int tmp_n=0;

	for (int L = 0; L <= lmax; L++)
	{
		for (int N = 0; N < nchi[L]; N++)
		{
			ofs_in << " " << setw(12) << count + 1 << setw(3) << L << setw(3) << N;

			double* radial; // radial mesh
			double* psi; // radial local orbital
			double* psir;// psi * r
			double* rab;// dr

			// set the number of mesh and the interval distance.
			ofs_in << setw(8) << meshr << setw(8) << dr;

			radial = new double[meshr];
			psi = new double[meshr];
			psir = new double[meshr];
			rab = new double[meshr];

			for(int im=0; im<meshr; ++im)
			{
				radial[im]=0.0;
				psi[im]=0.0;
				psir[im]=0.0;
				rab[im]=0.0;
			}

			for (int ir = 0; ir < meshr; ir++)
			{
				rab[ir] = dr;
				// plus one because we can't read in r = 0 term now.
				// change ir+1 to ir, because we need psi(r==0) information.
				radial[ir] = ir * dr; //mohan 2010-04-19
			}

			// set the length of orbital
			ofs_in << setw(8) << radial[meshr - 1];

			// mohan update 2010-09-07
			bool find = false;
			if (my_rank == 0)
			{
				while (!find)
				{
					if (ifs.eof())
					{
						cout << " Can't find l="
							<< L << " n=" << N << " orbital." << endl;
						break;
					}

					ifs >> name1 >> name2>> name3;
					ifs >> tmp_it >> tmp_l >> tmp_n;
					assert( name1 == "Type" );
					if (L == tmp_l && N == tmp_n)
					{
						// meshr_read is different from meshr if meshr is even number.
						for (int ir = 0; ir < meshr_read; ir++)
						{
							ifs >> psi[ir];
							psir[ir] = psi[ir] * radial[ir];
						}
						find = true;
					}
					else
					{
						double no_use;
						for (int ir = 0; ir < meshr_read; ir++)
						{
							ifs >> no_use;
						}
					}
				}//end find
			}

#ifdef __MPI
			Parallel_Common::bcast_bool(find);
#endif
			if (!find)
			{
				WARNING_QUIT("LCAO_Orbitals::read_orb_file", "Can't find orbitals.");
			}

#ifdef __MPI
			Parallel_Common::bcast_double(psi, meshr_read);
			Parallel_Common::bcast_double(psir, meshr_read);
#endif

			// renormalize radial wave functions
			double* inner = new double[meshr]();
			for (int ir = 0; ir < meshr; ir++)
			{
				inner[ir] = psir[ir] * psir[ir];
			}
			double unit = 0.0;

			Integral::Simpson_Integral(meshr, inner, rab, unit);

			assert(unit>0.0);

			// check unit: \sum ( psi[r] * r )^2 = 1
			ofs_in << setprecision(3) << setw(12) << unit;

			for (int ir = 0; ir < meshr; ir++)
			{
				psi[ir] /= sqrt(unit);
				psir[ir] /= sqrt(unit);
			}

			for (int ir = 0; ir < meshr; ir++)
			{
				inner[ir] = psir[ir] * psir[ir];
			}
			Integral::Simpson_Integral(meshr, inner, rab, unit);
			delete[] inner;
			ofs_in << setw(12) << unit << endl;

			ao[it].phiLN[count].set_orbital_info(
                orb_label,
                it, //type
				L, //angular momentum L
				N, // number of orbitals of this L
				meshr, // number of radial mesh
				rab,
				radial,// radial mesh value(a.u.)
				Numerical_Orbital_Lm::Psi_Type::Psi,// psi type next
				psi, // radial wave function
				this->kmesh,
				this->dk,
				this->dr_uniform,
				true,
				true,
				force_flag); // delta k mesh in reciprocal space

			delete[] radial;
			delete[] rab;
			delete[] psir;
			delete[] psi;

			++count;
		}
	}
	ao[it].set_orbital_info(
        it, // type
        orb_label, // label	
		lmax,
		nchi,
		total_nchi); //copy twice !
	
	delete[] nchi;
	return;
}

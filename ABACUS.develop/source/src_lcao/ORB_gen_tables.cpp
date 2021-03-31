#include "src_pw/global.h"
#include "ORB_read.h"
#include "ORB_gen_tables.h"
#include "src_global/ylm.h"

// here is a member of ORB_gen_tables class
ORB_gen_tables UOT;

ORB_gen_tables::ORB_gen_tables(){}
ORB_gen_tables::~ORB_gen_tables(){}

// call in hamilt_linear::init_before_ions.
void ORB_gen_tables::gen_tables( const int &job0 )
{
	TITLE("ORB_gen_tables","gen_tables");
	timer::tick("ORB_gen_tables","gen_tables",'C');

	ofs_running << "\n SETUP THE TWO-CENTER INTEGRATION TABLES" << endl;
	
	//=========================================
	// (1) MOT: make overlap table.
	//=========================================
	MOT.allocate(
		ORB.get_ntype(),// number of atom types
        ORB.get_lmax(),// max L used to calculate overlap
        ORB.get_kmesh(), // kpoints, for integration in k space
        ORB.get_Rmax(),// max value of radial table
        ORB.get_dR(),// delta R, for making radial table
        ORB.get_dk() ); // delta k, for integration in k space

	tbeta.allocate(
		ORB.get_ntype(),// number of atom types
        ORB.get_lmax(),// max L used to calculate overlap
        ORB.get_kmesh(), // kpoints, for integration in k space
        ORB.get_Rmax(),// max value of radial table
        ORB.get_dR(),// delta R, for making radial table
        ORB.get_dk() ); // delta k, for integration in k space

	//caoyu add 2021-03-18
	if (INPUT.out_descriptor && BASIS_TYPE == "lcao") {
		talpha.allocate(
			ORB.get_ntype(),// number of atom types
			ORB.get_lmax(),// max L used to calculate overlap
			ORB.get_kmesh(), // kpoints, for integration in k space
			ORB.get_Rmax(),// max value of radial table
			ORB.get_dR(),// delta R, for making radial table
			ORB.get_dk()); // delta k, for integration in k space
	}

	// OV: overlap
	MOT.init_OV_Tpair();
	MOT.init_OV_Opair();

	// NL: nonlocal
	tbeta.init_NL_Tpair();
	tbeta.init_NL_Opair(); // add 2009-5-8

	//caoyu add 2021-03-18
	// DS: Descriptor
	if (INPUT.out_descriptor && BASIS_TYPE == "lcao") {
		talpha.init_DS_Opair();
		talpha.init_DS_2Lplus1();
	}

	//=========================================
	// (2) init Ylm Coef
	//=========================================
	//liaochen add 2010/4/29
	Ylm::set_coefficients ();

	// Peize Lin update 2016-01-26
	int Lmax_used, Lmax;
	MOT.init_Table_Spherical_Bessel (2,1, Lmax_used, Lmax);
	
	//calculate S(R) for interpolation
	MOT.init_Table(job0);
	tbeta.init_Table_Beta( MOT.pSB );// add 2009-5-8

	//caoyu add 2021-03-18
	if (INPUT.out_descriptor && BASIS_TYPE == "lcao") {
		talpha.init_Table_Alpha(MOT.pSB);
		talpha.print_Table_DSR();	
	}

	//=========================================
	// (3) make Gaunt coefficients table
	//=========================================

	const int lmax = (Lmax_used-1) / 2 ;
	//MGT.init_Ylm_Gaunt(ORB.get_lmax()+1, 0.0,PI,0.0,TWO_PI);
	MGT.init_Gaunt_CH( lmax );
	//MGT.init_Gaunt(ORB.get_lmax()+1);
	MGT.init_Gaunt( lmax );



	timer::tick("ORB_gen_tables","gen_tables",'C');
	return;
}

void ORB_gen_tables::snap_psibeta(
	double nlm[],
	const int& job,
	const Vector3<double> &R1,
	const int &T1,
	const int &L1,
	const int &m1,
	const int &N1,
	const Vector3<double> &R2,
	const int &T2,
	const int &L2,
	const int &m2,
	const int &N2,
	const Vector3<double> &R0,// The projector.
	const int &T0,
	complex<double> *nlm1,
	const int is) const
{
	//TITLE ("ORB_gen_tables","snap_psibeta");

	//optimized by zhengdy-soc
	if(NSPIN==4 && ORB.Beta[T0].get_count_soc(is)==0) 
	{
		return;
	}

	timer::tick ("ORB_gen_tables","snap_psibeta",'X');

	bool has_so = 0;
	if(ORB.Beta[T0].get_count_soc(0)>0 ) has_so = 1;

	const int nproj = ORB.nproj[T0];
	bool *calproj = new bool[nproj];
	int* rmesh1 = new int[nproj];
	int* rmesh2 = new int[nproj];

	//rcut of orbtials and projectors
	const double Rcut1 = ORB.Phi[T1].getRcut();
	const double Rcut2 = ORB.Phi[T2].getRcut();
	
	//in our calculation, we always put orbital phi at the left side of <phi|beta>
	//because <phi|beta> = <beta|phi>
	const Vector3<double> dRa = (R0-R1)*this->lat0 ; 
	const Vector3<double> dRb = (R0-R2)*this->lat0 ;
	
	double distance10 = dRa.norm();
	double distance20 = dRb.norm();

	// mohan add 2011-03-10
	// because the table length is different accordint to each length
	// of projector, so sometimes some shorter projectors need not be 
	// calculated.
	bool all_out = true;
	for(int ip=0; ip<nproj; ip++)
	{
		const double Rcut0 = ORB.Beta[T0].Proj[ip].getRcut();
		if( distance10 > (Rcut1 + Rcut0) || distance20 > (Rcut2 + Rcut0) )  
		{
			calproj[ip] = false;
		}
		else
		{
			all_out = false;
			calproj[ip] = true;
			//length of table for interpolation
			rmesh1[ip] = tbeta.get_rmesh(Rcut1, Rcut0);
			rmesh2[ip] = tbeta.get_rmesh(Rcut2, Rcut0);
		}
	}

	if(all_out)
	{
		delete[] calproj;
		delete[] rmesh1;
		delete[] rmesh2;
		timer::tick ("ORB_gen_tables","snap_psibeta",'X');
		return;
	}


	//FOR INTERPOLATION
	double* curr; //current pointer
	int iqa, iqb;
	double psa, psb;
	double x0a,x1a,x2a,x3a,x123a,x120a,x032a,x031a;
	double x0b,x1b,x2b,x3b,x123b,x120b,x032b,x031b;
	
	psa = distance10 / tbeta.dr;
	iqa = static_cast<int>(psa);
   	x0a = psa - static_cast<double>(iqa);
  	x1a = 1.0 - x0a;
   	x2a = 2.0 - x0a;
    x3a = 3.0 - x0a;
	x123a = x1a*x2a*x3a/6.0;
	x120a = x1a*x2a*x0a/6.0;
	x032a = x0a*x3a*x2a/2.0;
	x031a = x0a*x3a*x1a/2.0;
	
	psb = distance20 / tbeta.dr;
	iqb = (int) psb;
   	x0b = psb - (double)iqb ;
  	x1b = 1.0 - x0b;
   	x2b = 2.0 - x0b;
    x3b = 3.0 - x0b;
	x123b = x1b*x2b*x3b/6.0;
	x120b = x1b*x2b*x0b/6.0;
	x032b = x0b*x3b*x2b/2.0;
	x031b = x0b*x3b*x1b/2.0;
	
	//UNIT VECTOR
			
	//double unit_vec_dRa[3];
	//unit_vec_dRa[0] = dRa.x;
	//unit_vec_dRa[1] = dRa.y;
	//unit_vec_dRa[2] = dRa.z;
	
	double unit_vec_dRb[3];
	unit_vec_dRb[0] = dRb.x;
	unit_vec_dRb[1] = dRb.y;
	unit_vec_dRb[2] = dRb.z;
	
	//special case for R = 0;
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;

	if(distance10 < tiny1) distance10 += tiny1;
	if(distance20 < tiny1) distance20 += tiny1;
	

	// Find three dimension of 'Table_NR' '
	// Notice!!! T1 must be orbital, 
	// T0 must be nonlocal orbital
	// usage : pairs_nonlocal_type(T1 : orbital, T0 : projector);
	const int Tpair1 = tbeta.NL_Tpair(T1, T0);
	const int Tpair2 = tbeta.NL_Tpair(T2, T0);
	const int T1_2Lplus1 = tbeta.NL_L2plus1(T1, T0);
	const int T2_2Lplus1 = tbeta.NL_L2plus1(T2, T0);

	//gaunt index
	const int gindex1 = L1*L1+m1;
	const int gindex2 = L2*L2+m2;

	// Peize Lin change rlya, rlyb, grlyb 2016-08-26
	vector<double> rlya;
	vector<double> rlyb;
	vector<vector<double>> grlyb;
	
	Ylm::rl_sph_harm (T1_2Lplus1-1, dRa.x, dRa.y, dRa.z, rlya);
	if (job == 0) 
	{
		Ylm::rl_sph_harm (T2_2Lplus1-1, dRb.x, dRb.y, dRb.z, rlyb);
	}
	else 
	{
		Ylm::grad_rl_sph_harm (T2_2Lplus1-1, dRb.x, dRb.y, dRb.z, rlyb, grlyb);
	}
	//==============================================================================
	// Formula :                         T1       T0          T0        T2
	// sum_{L0}sum_{m0}
	// 			D_{L0,L0} <psi1_{L1,N1}|Beta_{L0,m0}><Beta_{L0,m0}|psi2_{L2,N2}>
	//==============================================================================
	//double v = 0.0;

	// mohan update 2011-03-07
	int n_projection =1;
	if(has_so) 
	{
		n_projection = ORB.Beta[T0].get_nproj_soc();
	}

	vector<complex<double>> term_a_nc(n_projection,{0,0});		// Peize Lin change ptr to vector at 2020.01.31
	vector<complex<double>> term_b_nc(n_projection,{0,0});		// Peize Lin change ptr to vector at 2020.01.31
	int ip = -1;

	for(int nb=0; nb<nproj; nb++)
	{
		if( !calproj[nb] ) continue;

		const int L0 = ORB.Beta[T0].getL_Beta(nb);
		//const int next_ip = 2* L0 +1;
	

//-------------------------------------------------------------------
// move iterations for psi1 and psi2 from cal_fvnl_dbeta 
// to here --- 2021/03/20 mohan chen
//-------------------------------------------------------------------


		// <psi1 | Beta>
		const int Opair1 = tbeta.NL_Opair(Tpair1, L1, N1, nb); 
		// <psi2 | Beta>
		const int Opair2 = tbeta.NL_Opair(Tpair2, L2, N2, nb); 
		
			
		for(int m0=0; m0<2*L0+1; m0++)
		{
			++ip;
			int gindex0 = L0*L0+m0;
			
			//loop of {lmn}
			double term_a = 0.0;
			double term_b = 0.0;
			double term_c[3] = {0,0,0};	
			
			//=============
			// FIRST PART	
			//=============
			for(int L=0; L<T1_2Lplus1; L++)
			{
				//triangle rule for gaunt coefficients
				int AL = L1 + L0;
				int SL = abs (L1 - L0);
				if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) continue;
			
				//prefac = (i)^{lphi - lbeta - l}
				//R0-R1 ==> <phi|beta>
				double i_exp = pow(-1.0, (L1-L0-L)/2);
				double rl1 = pow(distance10, L);			
				double Interp_Vnla = 0.0;
				if (distance10 > tiny2)
				{	
					curr = tbeta.Table_NR[0][Tpair1][Opair1][L];
					if( iqa >= rmesh1[nb]-4)
					{
						Interp_Vnla = 0.0;
					}
					else
					{
						Interp_Vnla = i_exp * (x123a*curr[iqa]+x120a*curr[iqa+3]+x032a*curr[iqa+1]-x031a*curr[iqa+2]);
					}
					Interp_Vnla /= rl1;
				}
				else 
				{
					Interp_Vnla = i_exp * tbeta.Table_NR[0][Tpair1][Opair1][L][0];
				}
	
				//------------------------------------------
				//  Overlap value = S_from_table * G * Ylm				
				//------------------------------------------
				for(int m=0; m<2*L+1; m++)
				{
					int gindexa = L*L+m;
					//double tmpGaunt = this->MGT.Get_Gaunt_SH(L1, m1, L0, m0, L, m); 
					double tmpGaunt = this->MGT.Gaunt_Coefficients (gindex1, gindex0, gindexa);
					term_a += tmpGaunt * Interp_Vnla * rlya[ MGT.get_lm_index(L, m) ];
				}
			} //end L

			//=============
			// SECOND PART	
			//=============
			for(int L=0; L<T2_2Lplus1; L++)
			{
				//triangle rule for gaunt coefficients
				int AL = L2 + L0;
				int SL = abs (L2 - L0);
				if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) continue;

				double Interp_Vnlb = 0.0;
				double Interp_Vnlc = 0.0;
				
				//prefac
				double i_exp = pow(-1.0, (L2-L0-L)/2);
				double rl2 = pow (distance20, L);	
				if (distance20 > tiny2)
				{
					curr = tbeta.Table_NR[0][Tpair2][Opair2][L];
   					
					if( iqb >= rmesh2[nb]-4) 
					{
						Interp_Vnlb = 0.0;
					}
					else 
					{
						Interp_Vnlb = i_exp * (x123b*curr[iqb]+x120b*curr[iqb+3]+x032b*curr[iqb+1]-curr[iqb+2]*x031b);
					}
					
					Interp_Vnlb /= rl2;
				}
				else 
				{
					Interp_Vnlb = i_exp * tbeta.Table_NR[0][Tpair2][Opair2][L][0];
				}

				
				if (job == 1) // 1 means calculate the derivative part.
				{
					if (distance20 > tiny2)
					{
						curr = tbeta.Table_NR[1][Tpair2][Opair2][L];
   					
						if( iqb >= rmesh2[nb]-4) 
						{
							Interp_Vnlc = 0.0;
						}
						else 
						{
							Interp_Vnlc = i_exp * (x123b*curr[iqb]+x120b*curr[iqb+3]+x032b*curr[iqb+1]-curr[iqb+2]*x031b);
						}
						Interp_Vnlc = Interp_Vnlc / pow(distance20, L) - Interp_Vnlb * L / distance20;
					}
					else 
					{
						Interp_Vnlc = 0.0;
					}
				}
				
				// sum up the second part.	
				for(int m=0; m<2*L+1; m++)
				{
					int gindexb = L*L+m;
					//double tmpGaunt = this->MGT.Get_Gaunt_SH(L0, m0, L2, m2, L, m);
					double tmpGaunt = this->MGT.Gaunt_Coefficients (gindex0, gindex2, gindexb);
					const int lm = MGT.get_lm_index(L, m);
					
					switch (job)
					{
						case 0:// calculate the overlap part.
						{
							term_b += tmpGaunt * Interp_Vnlb * rlyb[lm];
							break;
						}
						case 1: // calculate the derivative part.
						{
							double tt1 = tmpGaunt * Interp_Vnlc * rlyb[lm] / distance20;
							double tt2 = tmpGaunt * Interp_Vnlb;
										
							for(int ir = 0; ir < 3; ir++)
							{
								term_c[ir] += tt1 * unit_vec_dRb[ir] 
											+ tt2 * grlyb[lm][ir];
							}

							break;
						}
						default: break;
					}
				}// end m of SECOND PART
			}// end L of SECOND PART
		
		
			//added by zhengdy-soc, store them for soc case
			if(has_so)
			{
				term_a_nc[ip] = term_a;
				term_b_nc[ip] = term_b;
			}
		
			//===============================================
			// THIRD PART: SUM THE VALUE FROM ALL PROJECTS.
			//===============================================
			switch (job)
			{
				case 0://calculate the overlap part.
				{
					//nlm[0] += term_a * term_b * ORB.Beta[T0].getCoefficient_D(L0, L0);//LiuXh 2016-01-14
					if(!has_so) 
					{
						nlm[0] += term_a * term_b * ORB.Beta[T0].getCoefficient_D(nb, nb);//LiuXh 2016-01-14
					}
					break;
				}
				case 1: //calculate the derivative part.
				{
					for(int jr = 0; jr < 3; jr++) 
					{
						//nlm[jr] += term_c[jr] * term_a * ORB.Beta[T0].getCoefficient_D(L0, L0);//LiuXh 2016-01-14
						if(!has_so) 
						{
							nlm[jr] += term_c[jr] * term_a * ORB.Beta[T0].getCoefficient_D(nb, nb);//LiuXh 2016-01-14
						}
					}
					break;
				}
				default: break;
			}
		}//!m0
	}//!L0

	//zhengdy-soc, calculate non-local term
	if(has_so)
	{
		switch (job)
		{
			case 0://overlap part
				for(int no=0;no<ORB.Beta[T0].get_count_soc(is);no++)
				{
					const int p1 = ORB.Beta[T0].get_index1_soc(is, no);
					const int p2 = ORB.Beta[T0].get_index2_soc(is, no);
					if(NSPIN==4 && nlm1!=NULL)
					{
						nlm1[is] += term_a_nc[p1] * term_b_nc[p2] * ORB.Beta[T0].getCoefficient_D_so(is, p2, p1);
					}
					else if(NSPIN!=4)
					{
						nlm[0] += (term_a_nc[p1] * term_b_nc[p2] * ORB.Beta[T0].getCoefficient_D_so(0, p2, p1)).real();
					}
					else
					{
						WARNING_QUIT("ORB_gen_tables::snap_psibeta","Conflict! Didn't count non-local part");
					}
				}
				break;
			case 1://need to be added later
			{break;}
			default: break;
		}
	}

	delete[] calproj;
	delete[] rmesh1;
	delete[] rmesh2;

	timer::tick ("ORB_gen_tables","snap_psibeta",'X');
	return;
}

void ORB_gen_tables::snap_psipsi(
	double olm[],
	const int &job, //0, 1
	const char &dtype, // derivative type: S or T
	const Vector3<double> &R1,
    const int &T1,
    const int &L1,
    const int &m1,
    const int &N1,
    const Vector3<double> &R2,
    const int &T2,
    const int &L2,
    const int &m2,
    const int &N2,
	complex<double> *olm1)const
{
	//TITLE("ORB_gen_tables","snap_psipsi");
	//timer::tick ("ORB_gen_tables", "snap_psipsi");
	if(job != 0 && job != 1)
	{
		WARNING_QUIT("ORB_gen_tables::snap_psipsi","job must be equal to 0 or 1!");
	}
	
	Numerical_Orbital::set_position(R1, R2);
	assert(this->lat0>0.0);

	// (1) get distance between R1 and R2 (a.u.)
	// judge if there exist overlap
	double distance = Numerical_Orbital::get_distance()*this->lat0;
	
	const double Rcut1 = ORB.Phi[T1].getRcut();
	const double Rcut2 = ORB.Phi[T2].getRcut();

	if(job == 0) ZEROS(olm, 1);
	else if(job == 1) ZEROS(olm, 3);
	
	if( distance > (Rcut1 + Rcut2) ) return;
	
	//if distance == 0
	//\int psi(r) psi(r-R) dr independent of R if R == 0
	//distance += tiny1 avoid overflow during calculation
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;
	if(distance < tiny1) distance += tiny1;
	
	// (2) if there exist overlap, calculate the mesh number
	// between two atoms
	const int rmesh = this->MOT.get_rmesh(Rcut1, Rcut2);
	
	// (3) Find three dimension of 'Table_S' or 'Table_T'
	// dim1 : type pairs,
	// dim2 : radial orbital pairs,
	// dim3 : find lmax between T1 and T2, and get lmax*2+1
	const int dim1 = this->MOT.OV_Tpair(T1, T2);
	const int dim3 = this->MOT.OV_L2plus1(T1, T2); //2*lmax+1
	
	int dim2;
	if (T1 <= T2) dim2 = this->MOT.OV_Opair(dim1, L1, L2, N1, N2); 
	else dim2 = this->MOT.OV_Opair(dim1, L2, L1, N2, N1);
		
	// Find the needed Ylm(dR) dimension 
	const int nlm = dim3 * dim3; //(2lmax+1)*(2lmax+!)

	//Gaunt Index
	const int gindex1 = L1*L1+m1;
	const int gindex2 = L2*L2+m2;

	assert(nlm < 400);
	// Peize Lin change rly, grly 2016-08-26
	vector<double> rly;			
	vector<vector<double>> grly;
	
//	double *ylm = new double[nlm];
//	dR = R1 - R2;
	double arr_dR[3];
	arr_dR[0] = Numerical_Orbital::getX() * this->lat0;
	arr_dR[1] = Numerical_Orbital::getY() * this->lat0;
	arr_dR[2] = Numerical_Orbital::getZ() * this->lat0;
	
	//double xdr = arr_dR[0] / distance;
	//double ydr = arr_dR[1] / distance;
	//double zdr = arr_dR[2] / distance;
	
	//=======================
	// *r**l*Ylm_real
	// include its derivations
	//=======================
	if (job == 0) 
	{
//		Ylm::rlylm(dim3, arr_dR[0], arr_dR[1], arr_dR[2], rly);
//		Ylm::sph_harm (dim3-1, xdr, ydr, zdr, rly);
		Ylm::rl_sph_harm (dim3-1, arr_dR[0], arr_dR[1], arr_dR[2], rly);
	}
	else 
	{
//		Ylm::rlylm(dim3, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
		Ylm::grad_rl_sph_harm (dim3-1, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
	}

	switch( dtype )
	{
		case 'S':
		for (int L = 0; L < dim3; L++) //maxL = dim3-1
		{
			//===========================================================
			// triangle rule for L and sum of L, L1, L2 should be even
			//===========================================================
			int AL = L1 + L2;
			int SL = abs (L1 - L2);

			if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) continue;
			
			double Interp_Slm = 0.0;
			double Interp_dSlm = 0.0;
			double tmpOlm0 = 0.0;
			double tmpOlm1 = 0.0;
			
			// prefactor
			double i_exp = pow(-1.0, (L1 - L2 - L) / 2);
			double rl = pow (distance, L);

			if (distance > tiny2)
			{
				Interp_Slm = i_exp * Mathzone::Polynomial_Interpolation(
						MOT.Table_SR[0][dim1][dim2][L],	rmesh, MOT.dr, distance );
				Interp_Slm /= rl;
			}
			else // distance = 0.0; 
			{
				Interp_Slm = i_exp * MOT.Table_SR[0][dim1][dim2][L][0];
			}
				
			if (job == 1)//calculate the derivative.
			{
				if (distance > tiny2)
				{
					Interp_dSlm = i_exp * Mathzone::Polynomial_Interpolation(
						MOT.Table_SR[1][dim1][dim2][L], rmesh, MOT.dr, distance );
					Interp_dSlm = Interp_dSlm / pow (distance, L) - Interp_Slm * L / distance;
				}
				else 
				{
					Interp_dSlm = 0.0;
				}
			}
			
			for (int m = 0; m < 2*L+1; m++)
			{
				int gindex = L*L+m;
	//			double tmpGaunt1 = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
				double tmpGaunt = MGT.Gaunt_Coefficients (gindex1, gindex2, gindex);	
							
				tmpOlm0 = Interp_Slm * tmpGaunt;
	
				if (job == 1) 
				{
					tmpOlm1 = Interp_dSlm * tmpGaunt;
				}
				
				switch( job )
				{
					case 0: // calculate overlap.
					{	
						if(NSPIN!=4) olm[0] += tmpOlm0 * rly[ MGT.get_lm_index(L, m) ] ;
						else if(olm1!= NULL)
						{
							olm1[0] += tmpOlm0 * rly[ MGT.get_lm_index(L,m) ] ;
							olm1[1] += 0;//tmpOlm0 * (tmp(0,0)+tmp(0,1));
							olm1[2] += 0;//tmpOlm0 * (tmp(1,0)+tmp(1,1));
							olm1[3] += tmpOlm0 * rly[ MGT.get_lm_index(L,m) ] ;
							
						}
						else
						{
							WARNING_QUIT("ORB_gen_tables::snap_psipsi","something wrong!");
							
						}
					
						/*		
						if( abs ( tmpOlm0 * rly[ MGT.get_lm_index(L, m) ] ) > 1.0e-3 )
						{
						cout << " L=" << L << " m=" << m << " tmpOlm0=" << tmpOlm0 
						<< " rly=" << rly[ MGT.get_lm_index(L, m) ] 
						<< " r=" << olm[0]
						<< endl;
						}
						*/
						break;
					}
					case 1: // calculate gradient.
					{
						for(int ir = 0; ir < 3; ir++)
						{
							olm[ir] += tmpOlm0 * grly[ MGT.get_lm_index(L, m) ][ir]
									 + tmpOlm1 * rly[ MGT.get_lm_index(L, m) ] * arr_dR[ir] / distance;
						}
						break;
					}
					default: break;
				}
			}//!m
		}
		break;

		case 'T':
		for (int L = 0; L < dim3; L++)
		{
			int AL = L1 + L2;
			int SL = abs (L1 - L2);

			if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) continue;

			double Interp_Tlm, Interp_dTlm, tmpKem0, tmpKem1;
			Interp_Tlm = Interp_dTlm = tmpKem0 = tmpKem1 = 0.0;
			
			//pre-fac
			double i_exp = pow(-1.0, (L1 - L2 - L) / 2);

			double rl = pow (distance, L);
			if (distance > tiny2)
			{
				Interp_Tlm = i_exp * Mathzone::Polynomial_Interpolation(
						MOT.Table_TR[0][dim1][dim2][L],	rmesh, MOT.dr, distance );	
				Interp_Tlm /= rl;
			}
			else Interp_Tlm = i_exp * MOT.Table_TR[0][dim1][dim2][L][0];
				
			
			if (job == 1)
			{
				if (distance > tiny2)
				{
					Interp_dTlm = i_exp * Mathzone::Polynomial_Interpolation(
						MOT.Table_TR[1][dim1][dim2][L], rmesh, MOT.dr, distance );
					Interp_dTlm = Interp_dTlm / rl - Interp_Tlm * L / distance;
				}
				else Interp_dTlm = 0.0;
			}
			
			for (int m = 0; m < 2*L+1; m++)
			{
				int gindex = L*L+m;
			//	double tmpGaunt = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
				double tmpGaunt = MGT.Gaunt_Coefficients (gindex1, gindex2, gindex);
					
				tmpKem0 = Interp_Tlm * tmpGaunt;
				if (job == 1) 
				{
					tmpKem1 = Interp_dTlm * tmpGaunt;
				}
				
				switch( job )
				{
					case 0:
					{
						if(NSPIN!=4) olm[0] += tmpKem0 * rly[ MGT.get_lm_index(L, m) ];
						else if(olm1 != NULL)
						{
							olm1[0] += tmpKem0 * rly[ MGT.get_lm_index(L,m) ];
							olm1[1] += 0;//tmpKem0 * (tmp(0,0)+tmp(0,1));
							olm1[2] += 0;//tmpKem0 * (tmp(1,0)+tmp(1,1));
							olm1[3] += tmpKem0 * rly[ MGT.get_lm_index(L,m) ];
						}
						else
						{
							WARNING_QUIT("ORB_gen_tables::snap_psipsi","something wrong in T.");
						}
						break;
					}
					case 1: 
					{
						for(int ir = 0; ir < 3; ir++)
						{
							olm[ir] += tmpKem0 * grly[ MGT.get_lm_index(L, m) ][ir]
								    + tmpKem1 * rly[ MGT.get_lm_index(L, m) ] * arr_dR[ir] / distance;
						}
						break;
					}
					default: break;
				}
			}// end T: m
		}// end T: :
		break;
	}
//	timer::tick ("ORB_gen_tables", "snap_psipsi");
	return;
}

double ORB_gen_tables::get_distance( const Vector3<double> &R1, const Vector3<double> &R2)const
{
	assert( this->lat0 > 0.0);
	Vector3<double> dR = R1 - R2;
	return dR.norm() * this->lat0;	
}

//caoyu add 2021-03-17
void ORB_gen_tables::snap_psialpha(
	double olm[],
	const int& job,
	const Vector3<double>& R1,
	const int& T1,
	const int& L1,
	const int& m1,
	const int& N1,
	const Vector3<double>& R2,
	const int& T2,
	const int& L2,
	const int& m2,
	const int& N2,
	complex<double>* olm1,
	const int is) const
{

	if (job != 0 && job != 1)
	{
		WARNING_QUIT("ORB_gen_tables::snap_psialpha", "job must be equal to 0 or 1!");
	}

	Numerical_Orbital::set_position(R1, R2);
	assert(this->lat0 > 0.0);

	// (1) get distance between R1 and R2 (a.u.)
	// judge if there exist overlap
	double distance = Numerical_Orbital::get_distance() * this->lat0;

	const double Rcut1 = ORB.Phi[T1].getRcut();
	const double Rcut2 = ORB.Alpha[0].getRcut();

	if (job == 0) ZEROS(olm, 1);
	else if (job == 1) ZEROS(olm, 3);

	if (distance > (Rcut1 + Rcut2)) return;

	//if distance == 0
	//\int psi(r) psi(r-R) dr independent of R if R == 0
	//distance += tiny1 avoid overflow during calculation
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;
	if (distance < tiny1) distance += tiny1;

	// (2) if there exist overlap, calculate the mesh number
	// between two atoms
	const int rmesh = this->talpha.get_rmesh(Rcut1, Rcut2);

	// (3) Find three dimension of 'Table_DS'
	// dim1 : type pairs, equal to T1 here 
	// dim2 : radial orbital pairs,
	// dim3 : find lmax between T1 and T2, and get lmax*2+1
	const int dim1 = T1;
	int dim2 = this->talpha.DS_Opair(dim1, L1, L2, N1, N2);
	int dim3 = this->talpha.DS_2Lplus1[T1];

	//Gaunt Index
		const int gindex1 = L1 * L1 + m1;
	const int gindex2 = L2 * L2 + m2;

	// Peize Lin change rly, grly 2016-08-26
	vector<double> rly;
	vector<vector<double>> grly;

	double arr_dR[3];
	arr_dR[0] = Numerical_Orbital::getX() * this->lat0;
	arr_dR[1] = Numerical_Orbital::getY() * this->lat0;
	arr_dR[2] = Numerical_Orbital::getZ() * this->lat0;

	//double xdr = arr_dR[0] / distance;
	//double ydr = arr_dR[1] / distance;
	//double zdr = arr_dR[2] / distance;

	//=======================
	// *r**l*Ylm_real
	// include its derivations
	//=======================
	if (job == 0)
	{
		Ylm::rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly);
	}
	else
	{
		Ylm::grad_rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
	}

	for (int L = 0; L < dim3; L++) //maxL = dim3-1
	{
		//===========================================================
		// triangle rule for L and sum of L, L1, L2 should be even
		//===========================================================
		int AL = L1 + L2;
		int SL = abs(L1 - L2);

		if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1)) continue;

		double Interp_Slm = 0.0;
		double Interp_dSlm = 0.0;
		double tmpOlm0 = 0.0;
		double tmpOlm1 = 0.0;

		// prefactor
		double i_exp = pow(-1.0, (L1 - L2 - L) / 2);
		double rl = pow(distance, L);

		if (distance > tiny2)
		{
			Interp_Slm = i_exp * Mathzone::Polynomial_Interpolation(
				talpha.Table_DSR[0][dim1][dim2][L], rmesh, MOT.dr, distance);
			Interp_Slm /= rl;
		}
		else // distance = 0.0; 
		{
			Interp_Slm = i_exp * talpha.Table_DSR[0][dim1][dim2][L][0];
		}

		if (job == 1)//calculate the derivative.
		{
			if (distance > tiny2)
			{
				Interp_dSlm = i_exp * Mathzone::Polynomial_Interpolation(
					talpha.Table_DSR[1][dim1][dim2][L], rmesh, MOT.dr, distance);
				Interp_dSlm = Interp_dSlm / pow(distance, L) - Interp_Slm * L / distance;
			}
			else
			{
				Interp_dSlm = 0.0;
			}
		}

		for (int m = 0; m < 2 * L + 1; m++)
		{
			int gindex = L * L + m;
			//			double tmpGaunt1 = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
			double tmpGaunt = MGT.Gaunt_Coefficients(gindex1, gindex2, gindex);

			tmpOlm0 = Interp_Slm * tmpGaunt;

			if (job == 1)
			{
				tmpOlm1 = Interp_dSlm * tmpGaunt;
			}

			switch (job)
			{
			case 0: // calculate overlap.
			{
				if (NSPIN != 4) olm[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
				else if (olm1 != NULL)
				{
					olm1[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
					olm1[1] += 0;//tmpOlm0 * (tmp(0,0)+tmp(0,1));
					olm1[2] += 0;//tmpOlm0 * (tmp(1,0)+tmp(1,1));
					olm1[3] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];

				}
				else
				{
					WARNING_QUIT("ORB_gen_tables::snap_psialpha", "something wrong!");

				}

				/*
				if( abs ( tmpOlm0 * rly[ MGT.get_lm_index(L, m) ] ) > 1.0e-3 )
				{
				cout << " L=" << L << " m=" << m << " tmpOlm0=" << tmpOlm0
				<< " rly=" << rly[ MGT.get_lm_index(L, m) ]
				<< " r=" << olm[0]
				<< endl;
				}
				*/
				break;
			}
			case 1: // calculate gradient.
			{
				for (int ir = 0; ir < 3; ir++)
				{
					olm[ir] += tmpOlm0 * grly[MGT.get_lm_index(L, m)][ir]
						+ tmpOlm1 * rly[MGT.get_lm_index(L, m)] * arr_dR[ir] / distance;
				}
				break;
			}
			default: break;
			}
		}//!m
	}

	return;
}

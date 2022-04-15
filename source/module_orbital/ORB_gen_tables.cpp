#include "ORB_read.h"
#include "ORB_gen_tables.h"
#include "../module_base/ylm.h"
#include "../module_base/math_polyint.h"
#include "../module_base/timer.h"

namespace GlobalC
{
///here is a member of ORB_gen_tables class
ORB_gen_tables UOT;
}

ORB_gen_tables::ORB_gen_tables() {}
ORB_gen_tables::~ORB_gen_tables() {}

/// call in hamilt_linear::init_before_ions.
void ORB_gen_tables::gen_tables(
	std::ofstream &ofs_in,
	const int &job0,
	LCAO_Orbitals &orb,
	const int &Lmax_exx,
	const bool &deepks_setorb,
	const int &nprojmax, 
	const int* nproj,
	const Numerical_Nonlocal* beta_)
{
	ModuleBase::TITLE("ORB_gen_tables", "gen_tables");
	ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");

	ofs_in << "\n SETUP THE TWO-CENTER INTEGRATION TABLES" << std::endl;

	//////////////////////////////
	/// (1) MOT: make overlap table.
	//////////////////////////////
	MOT.allocate(
		orb.get_ntype(),
		orb.get_lmax(),	 
		orb.get_kmesh(), 
		orb.get_Rmax(),	
		orb.get_dR(),	 
		orb.get_dk());	 

	tbeta.allocate(
		orb.get_ntype(),
		orb.get_lmax(),	
		orb.get_kmesh(), 
		orb.get_Rmax(),	
		orb.get_dR(),
		orb.get_dk());

	//caoyu add 2021-03-18
	//mohan update 2021-04-22
	if (deepks_setorb)
	{
		talpha.allocate(
			orb.get_ntype(), 
			orb.get_lmax(),	
			orb.get_kmesh(),
			orb.get_Rmax(),
			orb.get_dR(),
			orb.get_dk());
	}

	// OV: overlap
	MOT.init_OV_Tpair(orb);
	MOT.init_OV_Opair(orb);

	// NL: nonlocal
	tbeta.init_NL_Tpair(orb.Phi, beta_);
	tbeta.init_NL_Opair(orb, nprojmax, nproj); // add 2009-5-8

	//caoyu add 2021-03-18
	// DS: Descriptor
	if (deepks_setorb)
	{
		talpha.init_DS_Opair(orb);
		talpha.init_DS_2Lplus1(orb);
	}

	//////////////////////////////
	/// (2) init Ylm Coef
	//////////////////////////////
	//liaochen add 2010/4/29
	ModuleBase::Ylm::set_coefficients();

	// PLEASE add explanations for all options of 'orb_num' and 'mode'
	// mohan add 2021-04-03
	// Peize Lin update 2016-01-26
#ifdef __ORBITAL
	int orb_num = 4;
#else
	int orb_num = 2; //
#endif
	int mode = 1;	 // 1: <phi|phi> and <phi|beta>
	int Lmax_used = 0;
	int Lmax = 0;


	MOT.init_Table_Spherical_Bessel(orb_num, mode, Lmax_used, Lmax, Lmax_exx, orb, beta_);

	//calculate S(R) for interpolation
	MOT.init_Table(job0, orb);
	tbeta.init_Table_Beta(MOT.pSB, orb.Phi, beta_, nproj); // add 2009-5-8

	//caoyu add 2021-03-18
	if (deepks_setorb)
	{
		talpha.init_Table_Alpha(MOT.pSB, orb);
		if(GlobalV::deepks_out_unittest) talpha.print_Table_DSR(orb);
	}

	/////////////////////////////
	/// (3) make Gaunt coefficients table
	/////////////////////////////

	const int lmax = (Lmax_used - 1) / 2;
	//MGT.init_Ylm_Gaunt(orb.get_lmax()+1, 0.0,PI,0.0,ModuleBase::TWO_PI);
	MGT.init_Gaunt_CH(lmax);
	//MGT.init_Gaunt(orb.get_lmax()+1);
	MGT.init_Gaunt(lmax);

	ModuleBase::timer::tick("ORB_gen_tables", "gen_tables");
	return;
}

void ORB_gen_tables::snap_psibeta_half(
	const LCAO_Orbitals &orb,
	const InfoNonlocal &infoNL_,
	std::vector<std::vector<double>> &nlm,
	const ModuleBase::Vector3<double> &R1,
	const int &T1,
	const int &L1,
	const int &m1,
	const int &N1,
	const ModuleBase::Vector3<double> &R0, // The projector.
	const int &T0,
	const bool &calc_deri)const // mohan add 2021-04-25)
{
	ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta_half");

	//find number of projectors on atom R0
	const int nproj = infoNL_.nproj[T0];
	assert(nproj>0); // mohan add 2021-04-25

	std::vector<bool> calproj;
	calproj.resize(nproj);
	std::vector<int> rmesh1;
	rmesh1.resize(nproj);

	if(calc_deri)
	{
		nlm.resize(4);
	}
	else
	{
		nlm.resize(1);
	}

	//Count number of projectors (l,m)
	int natomwfc = 0;
	for (int ip = 0;ip < nproj;ip++)
	{
		//============================
		// Use pseudo-atomic orbitals
		//============================
		
		const int L0 = infoNL_.Beta[T0].Proj[ip].getL(); // mohan add 2021-05-07
		natomwfc += 2* L0 +1;
	}

	for(int dim=0;dim<nlm.size();dim++)
	{
		nlm[dim].resize(natomwfc);
		for (auto &x : nlm[dim])
		{
    		x=0.0;
		}
	}

	//rcut of orbtials and projectors
	//in our calculation, we always put orbital phi at the left side of <phi|beta>
	//because <phi|beta> = <beta|phi>
	const double Rcut1 = orb.Phi[T1].getRcut();
	const ModuleBase::Vector3<double> dRa = (R0 - R1) * this->lat0;
	double distance10 = dRa.norm();

	bool all_out = true;
	for (int ip = 0; ip < nproj; ip++)
	{
		const double Rcut0 = infoNL_.Beta[T0].Proj[ip].getRcut();
		if (distance10 > (Rcut1 + Rcut0))
		{
			calproj[ip] = false;
		}
		else
		{
			all_out = false;
			calproj[ip] = true;
			//length of table for interpolation
			rmesh1[ip] = tbeta.get_rmesh(Rcut1, Rcut0);
		}
	}

	if (all_out)
	{
		ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta_half");
		return;
	}

	//FOR INTERPOLATION
	double *curr; //current pointer

	double psa = distance10 / tbeta.dr;
	int iqa = static_cast<int>(psa);
	double x0a = psa - static_cast<double>(iqa);
	double x1a = 1.0 - x0a;
	double x2a = 2.0 - x0a;
	double x3a = 3.0 - x0a;
	double x123a = x1a * x2a * x3a / 6.0;
	double x120a = x1a * x2a * x0a / 6.0;
	double x032a = x0a * x3a * x2a / 2.0;
	double x031a = x0a * x3a * x1a / 2.0;

	double unit_vec_dRa[3];
	unit_vec_dRa[0] = dRa.x;
	unit_vec_dRa[1] = dRa.y;
	unit_vec_dRa[2] = dRa.z;

	//special case for R = 0;
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;

	if (distance10 < tiny1)
	{
		distance10 += tiny1;
	}

	// Find three dimension of 'Table_NR' '
	// Notice!!! T1 must be orbital,
	// T0 must be nonlocal orbital
	// usage : pairs_nonlocal_type(T1 : orbital, T0 : projector);
	const int Tpair1 = tbeta.NL_Tpair(T1, T0);
	const int T1_2Lplus1 = tbeta.NL_L2plus1(T1, T0);

	//gaunt index
	const int gindex1 = L1 * L1 + m1;

	// Peize Lin change rlya, rlyb, grlyb 2016-08-26
	std::vector<double> rlya;
	std::vector<std::vector<double>> grlya;

	if(!calc_deri)
	{
		ModuleBase::Ylm::rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya);
	}
	else
	{
		ModuleBase::Ylm::grad_rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya, grlya);
	}

	//////////////////////////////////////////////////////////////////////////
	/// Formula :               T1            T0
	/// \f[
	/// 			<\psi1_{L1,N1}|\Beta_{L0,m0}>
	///\f]
	//////////////////////////////////////////////////////////////////////////

	int index = 0; //(l,m index of projector)
	for (int nb = 0; nb < nproj; nb++)
	{
		const int L0 = infoNL_.Beta[T0].Proj[nb].getL(); // mohan add 2021-05-07
		if (!calproj[nb])
		{
			index += 2*L0 + 1;
			continue;
		}

		//const int L0 = orb.Beta[T0].getL_Beta(nb); // mohan delete the variable 2021-05-07
		//const int next_ip = 2* L0 +1;

		// <psi1 | Beta>
		const int Opair1 = tbeta.NL_Opair(Tpair1, L1, N1, nb);

		for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
		{
			int gindex0 = L0 * L0 + m0;

			//loop of {lmn}
			double term_a = 0.0;
			double term_a_gr[3] = {0,0,0};

			for (int L = 0; L < T1_2Lplus1; L++)
			{
				//triangle rule for gaunt coefficients
				int AL = L1 + L0;
				int SL = abs(L1 - L0);
				if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				{
					continue;
				}

				//prefac = (i)^{lphi - lbeta - l}
				//R0-R1 ==> <phi|beta>
				double i_exp = pow(-1.0, (L1 - L0 - L) / 2);
				double rl1 = pow(distance10, L);
				double Interp_Vnla = 0.0;
				double Interp_Vnla_gr = 0.0;

//this part is for both deri and not deri
				if (distance10 > tiny2)
				{
					curr = tbeta.Table_NR[0][Tpair1][Opair1][L];
					if (iqa >= rmesh1[nb] - 4)
					{
						Interp_Vnla = 0.0;
					}
					else
					{
						Interp_Vnla = i_exp * (x123a * curr[iqa] 
						+ x120a * curr[iqa + 3] 
						+ x032a * curr[iqa + 1] 
						- x031a * curr[iqa + 2]);
					}
					Interp_Vnla /= rl1;
				}
				else
				{
					Interp_Vnla = i_exp * tbeta.Table_NR[0][Tpair1][Opair1][L][0];
				}

//this part is for deri only
				if(calc_deri)
				{
					if (distance10 > tiny2)
					{
						curr = tbeta.Table_NR[1][Tpair1][Opair1][L];

						if (iqa >= rmesh1[nb] - 4)
						{
							Interp_Vnla_gr = 0.0;
						}
						else
						{
							Interp_Vnla_gr = i_exp * (x123a * curr[iqa] 
							+ x120a * curr[iqa + 3] 
							+ x032a * curr[iqa + 1] 
							- x031a * curr[iqa + 2]);
						}
						Interp_Vnla_gr = Interp_Vnla_gr / pow(distance10, L) - Interp_Vnla * L / distance10;
					}
					else
					{
						Interp_Vnla_gr = 0.0;
					}
				}

				/////////////////////////////////////
				///  Overlap value = S_from_table * G * Ylm
				////////////////////////////////////
				for (int m = 0; m < 2 * L + 1; m++)
				{
					int gindexa = L * L + m;
					//double tmpGaunt = this->MGT.Get_Gaunt_SH(L1, m1, L0, m0, L, m);
					double tmpGaunt, tmpGaunt1;
					if(calc_deri)
					{
						tmpGaunt = this->MGT.Gaunt_Coefficients(gindex1, gindex0, gindexa);
						tmpGaunt1= this->MGT.Gaunt_Coefficients(gindex0, gindex1, gindexa);
					}
					else
					{
						tmpGaunt = this->MGT.Gaunt_Coefficients(gindex0, gindex1, gindexa);
					}
					const int lm = MGT.get_lm_index(L, m);
					
					term_a += tmpGaunt * Interp_Vnla * rlya[lm];
					if(calc_deri)
					{
						double tt1 = tmpGaunt1 * Interp_Vnla_gr * rlya[lm] / distance10;
						double tt2 = tmpGaunt1 * Interp_Vnla;

						for (int ir = 0; ir < 3; ir++)
						{
							term_a_gr[ir] += tt1 * unit_vec_dRa[ir] + tt2 * grlya[lm][ir];
						}
					}
				}
			} //end L

			//===============================================
			// THIRD PART: SAVE THE VALUE FOR ALL PROJECTS.
			//===============================================

			if(!calc_deri)
			{
				nlm[0][index] = term_a;
			}
			else
			{
				nlm[0][index] = term_a;
				for(int dim=1;dim<4;dim++)
				{
					nlm[dim][index] = term_a_gr[dim-1];
				}
			}

			index += 1;
		} // end m0
	}// end nb

	assert(index == natomwfc);
	ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta_half");

	return;
}	

void ORB_gen_tables::snap_psibeta(
	const LCAO_Orbitals &orb,
	const InfoNonlocal& infoNL_,
	double nlm[],
	const int &job,
	const ModuleBase::Vector3<double> &R1,
	const int &T1,
	const int &L1,
	const int &m1,
	const int &N1,
	const ModuleBase::Vector3<double> &R2,
	const int &T2,
	const int &L2,
	const int &m2,
	const int &N2,
	const ModuleBase::Vector3<double> &R0, // The projector.
	const int &T0,
	const ModuleBase::matrix &dion, // mohan add 2021-04-25
	const int &nspin,
	const ModuleBase::ComplexArray &d_so, // mohan add 2021-05-07
	const int &count_soc, // mohan add 2021-05-07
	const int* index1_soc, // mohan add 2021-05-07
	const int* index2_soc, // mohan add 2021-05-07
	const int &nproj_in, // mohan add 2021-05-07
	std::complex<double> *nlm1,
	const int is) const
{
	//ModuleBase::TITLE ("ORB_gen_tables","snap_psibeta");

	//optimized by zhengdy-soc
	if (nspin == 4 && count_soc == 0)
	{
		return;
	}

	ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta");

	bool has_so = 0;
	if (count_soc > 0)
	{
		has_so = 1;
	}

	const int nproj = infoNL_.nproj[T0];
	assert(nproj>0); // mohan add 2021-04-25
	
	bool *calproj = new bool[nproj];
	int *rmesh1 = new int[nproj];
	int *rmesh2 = new int[nproj];

	//rcut of orbtials and projectors
	const double Rcut1 = orb.Phi[T1].getRcut();
	const double Rcut2 = orb.Phi[T2].getRcut();

	//in our calculation, we always put orbital phi at the left side of <phi|beta>
	//because <phi|beta> = <beta|phi>
	const ModuleBase::Vector3<double> dRa = (R0 - R1) * this->lat0;
	const ModuleBase::Vector3<double> dRb = (R0 - R2) * this->lat0;

	double distance10 = dRa.norm();
	double distance20 = dRb.norm();

	// mohan add 2011-03-10
	// because the table length is different accordint to each length
	// of projector, so sometimes some shorter projectors need not be
	// calculated.
	bool all_out = true;
	for (int ip = 0; ip < nproj; ip++)
	{
		const double Rcut0 = infoNL_.Beta[T0].Proj[ip].getRcut();
		if (distance10 > (Rcut1 + Rcut0) || distance20 > (Rcut2 + Rcut0))
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

	if (all_out)
	{
		delete[] calproj;
		delete[] rmesh1;
		delete[] rmesh2;
		ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta");
		return;
	}

	//FOR INTERPOLATION
	double *curr; //current pointer

	double psa = distance10 / tbeta.dr;
	int iqa = static_cast<int>(psa);
	double x0a = psa - static_cast<double>(iqa);
	double x1a = 1.0 - x0a;
	double x2a = 2.0 - x0a;
	double x3a = 3.0 - x0a;
	double x123a = x1a * x2a * x3a / 6.0;
	double x120a = x1a * x2a * x0a / 6.0;
	double x032a = x0a * x3a * x2a / 2.0;
	double x031a = x0a * x3a * x1a / 2.0;

	double psb = distance20 / tbeta.dr;
	int iqb = (int)psb;
	double x0b = psb - (double)iqb;
	double x1b = 1.0 - x0b;
	double x2b = 2.0 - x0b;
	double x3b = 3.0 - x0b;

	double x123b = x1b * x2b * x3b / 6.0;
	double x120b = x1b * x2b * x0b / 6.0;
	double x032b = x0b * x3b * x2b / 2.0;
	double x031b = x0b * x3b * x1b / 2.0;

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

	if (distance10 < tiny1)
	{
		distance10 += tiny1;
	}
	if (distance20 < tiny1)
	{
		distance20 += tiny1;
	}

	// Find three dimension of 'Table_NR' '
	// Notice!!! T1 must be orbital,
	// T0 must be nonlocal orbital
	// usage : pairs_nonlocal_type(T1 : orbital, T0 : projector);
	const int Tpair1 = tbeta.NL_Tpair(T1, T0);
	const int Tpair2 = tbeta.NL_Tpair(T2, T0);
	const int T1_2Lplus1 = tbeta.NL_L2plus1(T1, T0);
	const int T2_2Lplus1 = tbeta.NL_L2plus1(T2, T0);

	//gaunt index
	const int gindex1 = L1 * L1 + m1;
	const int gindex2 = L2 * L2 + m2;

	// Peize Lin change rlya, rlyb, grlyb 2016-08-26
	std::vector<double> rlya;
	std::vector<double> rlyb;
	std::vector<std::vector<double>> grlyb;

	ModuleBase::Ylm::rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya);
	if (job == 0)
	{
		ModuleBase::Ylm::rl_sph_harm(T2_2Lplus1 - 1, dRb.x, dRb.y, dRb.z, rlyb);
	}
	else
	{
		ModuleBase::Ylm::grad_rl_sph_harm(T2_2Lplus1 - 1, dRb.x, dRb.y, dRb.z, rlyb, grlyb);
	}
	//////////////////////////////////////////////////////////////////////////
	/// Formula :                         T1       T0          T0        T2
	/// \f[
	///	\sum_{ L0 }sum_{ m0 }
	/// 			D_{L0,L0} <\psi1_{L1,N1}|\Beta_{L0,m0}><\Beta_{L0,m0}|\psi2_{L2,N2}>
	///\f]
	//////////////////////////////////////////////////////////////////////////
	//double v = 0.0;

	// mohan update 2011-03-07
	int nprojections = 1;
	if (has_so)
	{
//		nprojections = orb.Beta[T0].get_nproj_soc();
		nprojections = nproj_in; // mohan add 2021-05-07 
	}

	std::vector<std::complex<double>> term_a_nc(nprojections, {0, 0}); // Peize Lin change ptr to std::vector at 2020.01.31
	std::vector<std::complex<double>> term_b_nc(nprojections, {0, 0}); // Peize Lin change ptr to std::vector at 2020.01.31
	int ip = -1;

	for (int nb = 0; nb < nproj; nb++)
	{
		if (!calproj[nb])
		{
			continue;
		}

		//const int L0 = orb.Beta[T0].getL_Beta(nb); // mohan delete the variable 2021-05-07
		const int L0 = infoNL_.Beta[T0].Proj[nb].getL(); // mohan add 2021-05-07
		//const int next_ip = 2* L0 +1;

		//////////////////////////////////////////////////////
		/// we should consider move iterations for psi1 and psi2 from cal_fvnl_dbeta
		/// to here --- 2021/03/20 mohan chen
		//////////////////////////////////////////////////////

		// <psi1 | Beta>
		const int Opair1 = tbeta.NL_Opair(Tpair1, L1, N1, nb);
		// <psi2 | Beta>
		const int Opair2 = tbeta.NL_Opair(Tpair2, L2, N2, nb);

		for (int m0 = 0; m0 < 2 * L0 + 1; m0++)
		{
			++ip;
			int gindex0 = L0 * L0 + m0;

			//loop of {lmn}
			double term_a = 0.0;
			double term_b = 0.0;
			double term_c[3] = {0, 0, 0};

			//=============
			// FIRST PART
			//=============
			for (int L = 0; L < T1_2Lplus1; L++)
			{
				//triangle rule for gaunt coefficients
				int AL = L1 + L0;
				int SL = abs(L1 - L0);
				if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				{
					continue;
				}

				//prefac = (i)^{lphi - lbeta - l}
				//R0-R1 ==> <phi|beta>
				double i_exp = pow(-1.0, (L1 - L0 - L) / 2);
				double rl1 = pow(distance10, L);
				double Interp_Vnla = 0.0;
				if (distance10 > tiny2)
				{
					curr = tbeta.Table_NR[0][Tpair1][Opair1][L];
					if (iqa >= rmesh1[nb] - 4)
					{
						Interp_Vnla = 0.0;
					}
					else
					{
						Interp_Vnla = i_exp * (x123a * curr[iqa] 
						+ x120a * curr[iqa + 3] 
						+ x032a * curr[iqa + 1] 
						- x031a * curr[iqa + 2]);
					}
					Interp_Vnla /= rl1;
				}
				else
				{
					Interp_Vnla = i_exp * tbeta.Table_NR[0][Tpair1][Opair1][L][0];
				}

				/////////////////////////////////////
				///  Overlap value = S_from_table * G * Ylm
				////////////////////////////////////
				for (int m = 0; m < 2 * L + 1; m++)
				{
					int gindexa = L * L + m;
					//double tmpGaunt = this->MGT.Get_Gaunt_SH(L1, m1, L0, m0, L, m);
					double tmpGaunt = this->MGT.Gaunt_Coefficients(gindex1, gindex0, gindexa);
					term_a += tmpGaunt * Interp_Vnla * rlya[MGT.get_lm_index(L, m)];
				}
			} //end L

			//=============
			// SECOND PART
			//=============
			for (int L = 0; L < T2_2Lplus1; L++)
			{
				//triangle rule for gaunt coefficients
				int AL = L2 + L0;
				int SL = abs(L2 - L0);
				if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				{
					continue;
				}

				double Interp_Vnlb = 0.0;
				double Interp_Vnlc = 0.0;

				//prefac
				double i_exp = pow(-1.0, (L2 - L0 - L) / 2);
				double rl2 = pow(distance20, L);

				if (distance20 > tiny2)
				{
					curr = tbeta.Table_NR[0][Tpair2][Opair2][L];

					if (iqb >= rmesh2[nb] - 4)
					{
						Interp_Vnlb = 0.0;
					}
					else
					{
						Interp_Vnlb = i_exp * (x123b * curr[iqb] 
						+ x120b * curr[iqb + 3] 
						+ x032b * curr[iqb + 1] 
						- curr[iqb + 2] * x031b);
					}

					Interp_Vnlb /= rl2;
				}
				else
				{
					Interp_Vnlb = i_exp * tbeta.Table_NR[0][Tpair2][Opair2][L][0];
				}// end if(distance20)

				if (job == 1) // 1 means calculate the derivative part.
				{
					if (distance20 > tiny2)
					{
						curr = tbeta.Table_NR[1][Tpair2][Opair2][L];

						if (iqb >= rmesh2[nb] - 4)
						{
							Interp_Vnlc = 0.0;
						}
						else
						{
							Interp_Vnlc = i_exp * (x123b * curr[iqb] 
							+ x120b * curr[iqb + 3] 
							+ x032b * curr[iqb + 1] 
							- curr[iqb + 2] * x031b);
						}
						Interp_Vnlc = Interp_Vnlc / pow(distance20, L) - Interp_Vnlb * L / distance20;
					}
					else
					{
						Interp_Vnlc = 0.0;
					}
				} // end job==1

				// sum up the second part.
				for (int m = 0; m < 2 * L + 1; m++)
				{
					int gindexb = L * L + m;
					//double tmpGaunt = this->MGT.Get_Gaunt_SH(L0, m0, L2, m2, L, m);
					double tmpGaunt = this->MGT.Gaunt_Coefficients(gindex0, gindex2, gindexb);
					const int lm = MGT.get_lm_index(L, m);

					switch (job)
					{
						case 0: // calculate the overlap part.
						{
							term_b += tmpGaunt * Interp_Vnlb * rlyb[lm];
							break;
						}
						case 1: // calculate the derivative part.
						{
							double tt1 = tmpGaunt * Interp_Vnlc * rlyb[lm] / distance20;
							double tt2 = tmpGaunt * Interp_Vnlb;

							for (int ir = 0; ir < 3; ir++)
							{
								term_c[ir] += tt1 * unit_vec_dRb[ir] + tt2 * grlyb[lm][ir];
							}

							break;
						}
						default:
							break;
					}
				} // end m of SECOND PART
			} // end L of SECOND PART

			//added by zhengdy-soc, store them for soc case
			if (has_so)
			{
				term_a_nc[ip] = term_a;
				term_b_nc[ip] = term_b;
			}

			//===============================================
			// THIRD PART: SUM THE VALUE FROM ALL PROJECTS.
			//===============================================
			switch (job)
			{
				case 0: //calculate the overlap part.
				{
					if (!has_so)
					{
						nlm[0] += term_a * term_b * dion(nb, nb); //LiuXh 2016-01-14
					}
					break;
				}
				case 1: //calculate the derivative part.
				{
					for (int jr = 0; jr < 3; jr++)
					{
						if (!has_so)
						{
							nlm[jr] += term_c[jr] * term_a * dion(nb, nb); //LiuXh 2016-01-14
						}
					}
					break;
				}
				default:
					break;
			}

		} // end m0
	}// end nb

	//zhengdy-soc, calculate non-local term
	if (has_so)
	{
		switch (job)
		{
			case 0: //overlap part
				for (int no = 0; no < count_soc; no++)
				{
					const int p1 = index1_soc[no];
					const int p2 = index2_soc[no]; 
					if (nspin == 4 && nlm1 != NULL)
					{
						nlm1[is] += term_a_nc[p1] * term_b_nc[p2] * d_so(is, p2, p1);
					}
					else if (nspin != 4)
					{
						nlm[0] += (term_a_nc[p1] * term_b_nc[p2] * d_so(0, p2, p1)).real();
					}
					else
					{
						ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psibeta", "Conflict! Didn't count non-local part");
					}
				}
				break;
			case 1: //need to be added later
				{
					break;
				}
			default:
				break;
		}
	}

	delete[] calproj;
	delete[] rmesh1;
	delete[] rmesh2;

	ModuleBase::timer::tick("ORB_gen_tables", "snap_psibeta");
	return;
}

void ORB_gen_tables::snap_psipsi(
	const LCAO_Orbitals &orb,
	double olm[],
	const int &job,	   //0, 1
	const char &dtype, // derivative type: S or T
	const ModuleBase::Vector3<double> &R1,
	const int &T1,
	const int &L1,
	const int &m1,
	const int &N1,
	const ModuleBase::Vector3<double> &R2,
	const int &T2,
	const int &L2,
	const int &m2,
	const int &N2,
	const int &nspin,
	std::complex<double> *olm1) const
{
	//ModuleBase::TITLE("ORB_gen_tables","snap_psipsi");
	//ModuleBase::timer::tick ("ORB_gen_tables", "snap_psipsi");
	if (job != 0 && job != 1)
	{
		ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psipsi", "job must be equal to 0 or 1!");
	}

	Numerical_Orbital::set_position(R1, R2);
	assert(this->lat0 > 0.0);

	/// (1) get distance between R1 and R2 (a.u.)
	/// judge if there exist overlap
	double distance = Numerical_Orbital::get_distance() * this->lat0;

	const double Rcut1 = orb.Phi[T1].getRcut();
	const double Rcut2 = (dtype == 'D' ? orb.Alpha[0].getRcut() : orb.Phi[T2].getRcut());	//caoyu modified 2021-05-08

	if (job == 0)
	{
		ModuleBase::GlobalFunc::ZEROS(olm, 1);
	}
	else if (job == 1)
	{
		ModuleBase::GlobalFunc::ZEROS(olm, 3);
	}

	if (distance > (Rcut1 + Rcut2))
		return;

	/// if distance == 0, 
	/// \f[ \int psi(r) psi(r-R)\f] dr independent of R if R == 0. 
	/// distance += tiny1 avoid overflow during calculation.
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;
	if (distance < tiny1)
		distance += tiny1;

	/// (2) if there exist overlap, calculate the mesh number
	/// between two atoms
	const int rmesh = (dtype == 'D' ? this->talpha.get_rmesh(Rcut1, Rcut2) : this->MOT.get_rmesh(Rcut1, Rcut2));	//caoyu modified 2021-05-08

	/// (3) Find three dimension of 'Table_S' or 'Table_T'.
	/// -dim1 : type pairs,
	/// -dim2 : radial orbital pairs,
	/// -dim3 : find lmax between T1 and T2, and get lmax*2+1
	int dim1, dim2, dim3;
	if (dtype == 'D')	//caoyu modified 2021-05-08
	{
		dim1 = T1;
		dim2 = this->talpha.DS_Opair(dim1, L1, L2, N1, N2);
		dim3 = this->talpha.DS_2Lplus1[T1];
	}
	else
	{
		dim1 = this->MOT.OV_Tpair(T1, T2);
		dim3 = this->MOT.OV_L2plus1(T1, T2); //2*lmax+1
		if (T1 <= T2)
		{
			dim2 = this->MOT.OV_Opair(dim1, L1, L2, N1, N2);
		}
		else
		{
			dim2 = this->MOT.OV_Opair(dim1, L2, L1, N2, N1);
		}
	}
	

	//Gaunt Index
	const int gindex1 = L1 * L1 + m1;
	const int gindex2 = L2 * L2 + m2;

	// Peize Lin change rly, grly 2016-08-26
	std::vector<double> rly;
	std::vector<std::vector<double>> grly;

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
		ModuleBase::Ylm::rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly);
	}
	else
	{
		//		Ylm::rlylm(dim3, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
		ModuleBase::Ylm::grad_rl_sph_harm(dim3 - 1, arr_dR[0], arr_dR[1], arr_dR[2], rly, grly);
	}

	switch (dtype)
	{
	case 'S':
		for (int L = 0; L < dim3; L++) //maxL = dim3-1
		{
			//===========================================================
			// triangle rule for L and sum of L, L1, L2 should be even
			//===========================================================
			int AL = L1 + L2;
			int SL = abs(L1 - L2);

			if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				continue;

			double Interp_Slm = 0.0;
			double Interp_dSlm = 0.0;
			double tmpOlm0 = 0.0;
			double tmpOlm1 = 0.0;

			// prefactor
			double i_exp = pow(-1.0, (L1 - L2 - L) / 2);
			double rl = pow(distance, L);

			if (distance > tiny2)
			{
				Interp_Slm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
					MOT.Table_SR[0][dim1][dim2][L], rmesh, MOT.dr, distance);
				Interp_Slm /= rl;
			}
			else // distance = 0.0;
			{
				Interp_Slm = i_exp * MOT.Table_SR[0][dim1][dim2][L][0];
			}

			if (job == 1) //calculate the derivative.
			{
				if (distance > tiny2)
				{
					Interp_dSlm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
						MOT.Table_SR[1][dim1][dim2][L], rmesh, MOT.dr, distance);
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
					if (nspin != 4)
					{
						olm[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
					}
					else if (olm1 != NULL)
					{
						olm1[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
						olm1[1] += 0; //tmpOlm0 * (tmp(0,0)+tmp(0,1));
						olm1[2] += 0; //tmpOlm0 * (tmp(1,0)+tmp(1,1));
						olm1[3] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
					}
					else
					{
						ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psipsi", "something wrong!");
					}

					/*
						if( abs ( tmpOlm0 * rly[ MGT.get_lm_index(L, m) ] ) > 1.0e-3 )
						{
						std::cout << " L=" << L << " m=" << m << " tmpOlm0=" << tmpOlm0
						<< " rly=" << rly[ MGT.get_lm_index(L, m) ]
						<< " r=" << olm[0]
						<< std::endl;
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
				default:
					break;
				}
			} //m
		}
		break;

	case 'T':
		for (int L = 0; L < dim3; L++)
		{
			int AL = L1 + L2;
			int SL = abs(L1 - L2);

			if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				continue;

			double Interp_Tlm, Interp_dTlm, tmpKem0, tmpKem1;
			Interp_Tlm = Interp_dTlm = tmpKem0 = tmpKem1 = 0.0;

			//pre-fac
			double i_exp = pow(-1.0, (L1 - L2 - L) / 2);

			double rl = pow(distance, L);
			if (distance > tiny2)
			{
				Interp_Tlm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
					MOT.Table_TR[0][dim1][dim2][L], rmesh, MOT.dr, distance);
				Interp_Tlm /= rl;
			}
			else
				Interp_Tlm = i_exp * MOT.Table_TR[0][dim1][dim2][L][0];

			if (job == 1)
			{
				if (distance > tiny2)
				{
					Interp_dTlm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
						MOT.Table_TR[1][dim1][dim2][L], rmesh, MOT.dr, distance);
					Interp_dTlm = Interp_dTlm / rl - Interp_Tlm * L / distance;
				}
				else
					Interp_dTlm = 0.0;
			}

			for (int m = 0; m < 2 * L + 1; m++)
			{
				int gindex = L * L + m;
				//	double tmpGaunt = MGT.Get_Gaunt_SH(L1, m1, L2, m2, L, m);
				double tmpGaunt = MGT.Gaunt_Coefficients(gindex1, gindex2, gindex);

				tmpKem0 = Interp_Tlm * tmpGaunt;
				if (job == 1)
				{
					tmpKem1 = Interp_dTlm * tmpGaunt;
				}

				switch (job)
				{
				case 0:
				{
					if (nspin != 4)
					{
						olm[0] += tmpKem0 * rly[MGT.get_lm_index(L, m)];
					}
					else if (olm1 != NULL)
					{
						olm1[0] += tmpKem0 * rly[MGT.get_lm_index(L, m)];
						olm1[1] += 0; //tmpKem0 * (tmp(0,0)+tmp(0,1));
						olm1[2] += 0; //tmpKem0 * (tmp(1,0)+tmp(1,1));
						olm1[3] += tmpKem0 * rly[MGT.get_lm_index(L, m)];
					}
					else
					{
						ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psipsi", "something wrong in T.");
					}
					break;
				}
				case 1:
				{
					for (int ir = 0; ir < 3; ir++)
					{
						olm[ir] += tmpKem0 * grly[MGT.get_lm_index(L, m)][ir] + tmpKem1 * rly[MGT.get_lm_index(L, m)] * arr_dR[ir] / distance;
					}
					break;
				}
				default:
					break;
				}
			} // end T: m
		}	  // end T: :
		break;
	case 'D'://caoyu add 2021-05-08
		for (int L = 0; L < dim3; L++) //maxL = dim3-1
		{
			//===========================================================
			// triangle rule for L and sum of L, L1, L2 should be even
			//===========================================================
			int AL = L1 + L2;
			int SL = abs(L1 - L2);

			if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
				continue;

			double Interp_Slm = 0.0;
			double Interp_dSlm = 0.0;
			double tmpOlm0 = 0.0;
			double tmpOlm1 = 0.0;

			// prefactor
			double i_exp = pow(-1.0, (L1 - L2 - L) / 2);
			double rl = pow(distance, L);

			if (distance > tiny2)
			{
				Interp_Slm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
					talpha.Table_DSR[0][dim1][dim2][L], rmesh, MOT.dr, distance);
				Interp_Slm /= rl;
			}
			else // distance = 0.0;
			{
				Interp_Slm = i_exp * talpha.Table_DSR[0][dim1][dim2][L][0];
			}

			if (job == 1) //calculate the derivative.
			{
				if (distance > tiny2)
				{
					Interp_dSlm = i_exp * ModuleBase::PolyInt::Polynomial_Interpolation(
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
					int nspin = 1; // mohan add 2021-05-07, currently deepks works only for nspin=1
					if (nspin != 4)
					{
						olm[0] += tmpOlm0 * rly[MGT.get_lm_index(L, m)];
					}
					else
					{
						ModuleBase::WARNING_QUIT("ORB_gen_tables::snap_psialpha", "deepks with GlobalV::NSPIN>1 has not implemented yet!");
					}
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
				default:
					break;
				}
			} //m
		}
	}
	//	ModuleBase::timer::tick ("ORB_gen_tables", "snap_psipsi");
	return;
}

double ORB_gen_tables::get_distance(const ModuleBase::Vector3<double> &R1, const ModuleBase::Vector3<double> &R2) const
{
	assert(this->lat0 > 0.0);
	ModuleBase::Vector3<double> dR = R1 - R2;
	return dR.norm() * this->lat0;
}

#ifdef __DEEPKS

void ORB_gen_tables::snap_psialpha_half(
		const LCAO_Orbitals& orb,
		std::vector<std::vector<double>> &nlm,
		const int& job,
		const ModuleBase::Vector3<double>& R1,
		const int& T1,
		const int& L1,
		const int& m1,
		const int& N1,
		const ModuleBase::Vector3<double>& R0, // The projector.
		const int& T0,
		const int& I0
	) const
{
	ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha_half");

    const int ln_per_atom = orb.Alpha[0].getTotal_nchi();
    assert(ln_per_atom > 0); 
	
	std::vector<bool> calproj;
	calproj.resize(ln_per_atom);
	std::vector<int> rmesh1;
	rmesh1.resize(ln_per_atom);

	if(job==0)
	{
		nlm.resize(1); //only energy
	}
	else if(job==1)
	{
		nlm.resize(4); //energy+force
	}

	int nproj = 0;
    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
    {
        for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
        {
			nproj += (2 * L0 + 1);
		}
	}

	for(int dim=0;dim<nlm.size();dim++)
	{
		nlm[dim].resize(nproj);
		for (auto &x : nlm[dim])
		{
    		x=0.0;
		}
	}

	//rcut of orbtials and projectors
	const double Rcut1 = orb.Phi[T1].getRcut();

	//in our calculation, we always put orbital phi at the left side of <phi|alpha>
	const ModuleBase::Vector3<double> dRa = (R0 - R1) * this->lat0;
	double distance10 = dRa.norm();

	bool all_out = true;
	for (int ip = 0; ip < ln_per_atom; ip++)
	{
		const double Rcut0 = orb.Alpha[0].getRcut();
		if (distance10 > (Rcut1 + Rcut0))
		{
			calproj[ip] = false;
		}
		else
		{
			all_out = false;
			calproj[ip] = true;
			//length of table for interpolation
			rmesh1[ip] = talpha.get_rmesh(Rcut1, Rcut0);
		}
	}

	if (all_out)
	{
		ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha_half");
		return;
	}

	//FOR INTERPOLATION
	double *curr; //current pointer

	double psa = distance10 / talpha.dr;
	int iqa = static_cast<int>(psa);
	double x0a = psa - static_cast<double>(iqa);
	double x1a = 1.0 - x0a;
	double x2a = 2.0 - x0a;
	double x3a = 3.0 - x0a;
	double x123a = x1a * x2a * x3a / 6.0;
	double x120a = x1a * x2a * x0a / 6.0;
	double x032a = x0a * x3a * x2a / 2.0;
	double x031a = x0a * x3a * x1a / 2.0;

	//UNIT VECTOR

	double unit_vec_dRa[3];
	unit_vec_dRa[0] = dRa.x;
	unit_vec_dRa[1] = dRa.y;
	unit_vec_dRa[2] = dRa.z;

	//special case for R = 0;
	const double tiny1 = 1e-12;
	const double tiny2 = 1e-10;

	if (distance10 < tiny1)
	{
		distance10 += tiny1;
	}

	// Find three dimension of 'Table_DSR' '
	const int Tpair1 = T1;
	const int T1_2Lplus1 = talpha.DS_2Lplus1[T1];

	//gaunt index
	const int gindex1 = L1 * L1 + m1;
	std::vector<double> rlya;
	std::vector<std::vector<double>> grlya;

	if (job == 0)
	{
		ModuleBase::Ylm::rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya);
	}
	else
	{
		ModuleBase::Ylm::grad_rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya, grlya);
	}

	//////////////////////////////////////////////////////////////////////////
	/// Formula :               T1            T0
	/// \f[
	/// 			<\psi1_{L1,N1}|\Beta_{L0,m0}>
	///\f]
	//////////////////////////////////////////////////////////////////////////

	int ip = 0; //for L0,N0,m0
    int nb = 0; //for L0, N0

    for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
    {
        for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
        {
            if (!calproj[nb])
            {
				ip += 2*L0 + 1;
				++nb;
                continue;
            }
            
            // <psi1 | Beta>
            const int Opair1 = talpha.DS_Opair(Tpair1, L1, L0, N1, N0);
 
            for (int m0 = 0;m0 < 2 * L0 + 1;++m0)
            {
				int gindex0 = L0 * L0 + m0;
				double term_a = 0.0;
				double term_a_gr[3] = {0,0,0};

				for (int L = 0; L < T1_2Lplus1; L++)
				{
					//triangle rule for gaunt coefficients
					int AL = L1 + L0;
					int SL = abs(L1 - L0);
					if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
					{
						continue;
					}

					//prefac = (i)^{lphi - lbeta - l}
					//R0-R1 ==> <phi|beta>
					double i_exp = pow(-1.0, (L1 - L0 - L) / 2);
					double rl1 = pow(distance10, L);
					double Interp_Vnla = 0.0;
					double Interp_Vnla_gr = 0.0;

					//this part is for both deri and not deri
					if (distance10 > tiny2)
					{
						curr = talpha.Table_DSR[0][Tpair1][Opair1][L];
						if (iqa >= rmesh1[nb] - 4)
						{
							Interp_Vnla = 0.0;
						}
						else
						{
							Interp_Vnla = i_exp * (x123a * curr[iqa] 
							+ x120a * curr[iqa + 3] 
							+ x032a * curr[iqa + 1] 
							- x031a * curr[iqa + 2]);
						}
						Interp_Vnla /= rl1;
					}
					else
					{
						Interp_Vnla = i_exp * talpha.Table_DSR[0][Tpair1][Opair1][L][0];
					}

					//this part is for deri only
					if(job==1)
					{
						if (distance10 > tiny2)
						{
							curr = talpha.Table_DSR[1][Tpair1][Opair1][L];

							if (iqa >= rmesh1[nb] - 4)
							{
								Interp_Vnla_gr = 0.0;
							}
							else
							{
								Interp_Vnla_gr = i_exp * (x123a * curr[iqa] 
								+ x120a * curr[iqa + 3] 
								+ x032a * curr[iqa + 1] 
								- x031a * curr[iqa + 2]);
							}
							Interp_Vnla_gr = Interp_Vnla_gr / pow(distance10, L) - Interp_Vnla * L / distance10;
						}
						else
						{
							Interp_Vnla_gr = 0.0;
						}
					}
						/////////////////////////////////////
						///  Overlap value = S_from_table * G * Ylm
						////////////////////////////////////
					for (int m = 0; m < 2 * L + 1; m++)
					{
						int gindexa = L * L + m;
						//double tmpGaunt = this->MGT.Get_Gaunt_SH(L1, m1, L0, m0, L, m);
						double tmpGaunt, tmpGaunt1;
						if(job==1)
						{
							tmpGaunt = this->MGT.Gaunt_Coefficients(gindex1, gindex0, gindexa);
							tmpGaunt1= this->MGT.Gaunt_Coefficients(gindex0, gindex1, gindexa);
						}
						else
						{
							tmpGaunt = this->MGT.Gaunt_Coefficients(gindex0, gindex1, gindexa);
						}
						const int lm = MGT.get_lm_index(L, m);
						
						term_a += tmpGaunt * Interp_Vnla * rlya[lm];
						if(job==1)
						{
							double tt1 = tmpGaunt1 * Interp_Vnla_gr * rlya[lm] / distance10;
							double tt2 = tmpGaunt1 * Interp_Vnla;
							for (int ir = 0; ir < 3; ir++)
							{
								term_a_gr[ir] += tt1 * unit_vec_dRa[ir] + tt2 * grlya[lm][ir];
							}
						}
					}
				}//end L

				//===============================================
				// THIRD PART: SAVE THE VALUE FOR ALL PROJECTS.
				//===============================================

				if(job==0)
				{
					nlm[0][ip] = term_a;
				}
				else if(job==1)
				{
					nlm[0][ip] = term_a;
					for(int dim=1;dim<4;dim++)
					{
						nlm[dim][ip] = term_a_gr[dim-1];
					}
				}

				ip+=1;
			}//end m0
			++nb;
		}//end N0
	}//end L0

	ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha_half");
	return;
}
//caoyu add 2021-08-30
void ORB_gen_tables::snap_psialpha(
    double nlm[],
    const int& job,
    const ModuleBase::Vector3<double>& R1,
    const int& T1,
    const int& L1,
    const int& m1,
    const int& N1,
    const ModuleBase::Vector3<double>& R2,
    const int& T2,
    const int& L2,
    const int& m2,
    const int& N2,
    const ModuleBase::Vector3<double>& R0, // The projector.
    const int& T0,
    const int& A0,  //gedm is related to specific atom
    ModuleBase::IntArray* inl_index,
    double** gedm    //Coefficient Matrix (non-diagonal)
    ) const
{
	//TITLE ("ORB_gen_tables","snap_psialpha")
	ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha");

    const int ln_per_atom = GlobalC::ORB.Alpha[0].getTotal_nchi();
    assert(ln_per_atom > 0); 
	
	bool *calproj = new bool[ln_per_atom];
	int *rmesh1 = new int[ln_per_atom];
	int *rmesh2 = new int[ln_per_atom];

	//rcut of orbtials and projectors
	const double Rcut1 = GlobalC::ORB.Phi[T1].getRcut();
	const double Rcut2 = GlobalC::ORB.Phi[T2].getRcut();

	//in our calculation, we always put orbital phi at the left side of <phi|alpha>
	const ModuleBase::Vector3<double> dRa = (R0 - R1) * this->lat0;
	const ModuleBase::Vector3<double> dRb = (R0 - R2) * this->lat0;

	double distance10 = dRa.norm();
	double distance20 = dRb.norm();

	// mohan add 2011-03-10
	// because the table length is different accordint to each length
	// of projector, so sometimes some shorter projectors need not be
	// calculated.
	bool all_out = true;
	for (int ip = 0; ip < ln_per_atom; ip++)
	{
		const double Rcut0 = GlobalC::ORB.Alpha[0].getRcut();
		if (distance10 > (Rcut1 + Rcut0) || distance20 > (Rcut2 + Rcut0))
		{
			calproj[ip] = false;
		}
		else
		{
			all_out = false;
			calproj[ip] = true;
			//length of table for interpolation
			rmesh1[ip] = talpha.get_rmesh(Rcut1, Rcut0);
			rmesh2[ip] = talpha.get_rmesh(Rcut2, Rcut0);
		}
	}

	if (all_out)
	{
		delete[] calproj;
		delete[] rmesh1;
		delete[] rmesh2;
		ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha");
		return;
	}

	//FOR INTERPOLATION
	double *curr; //current pointer

	double psa = distance10 / talpha.dr;
	int iqa = static_cast<int>(psa);
	double x0a = psa - static_cast<double>(iqa);
	double x1a = 1.0 - x0a;
	double x2a = 2.0 - x0a;
	double x3a = 3.0 - x0a;
	double x123a = x1a * x2a * x3a / 6.0;
	double x120a = x1a * x2a * x0a / 6.0;
	double x032a = x0a * x3a * x2a / 2.0;
	double x031a = x0a * x3a * x1a / 2.0;

	double psb = distance20 / talpha.dr;
	int iqb = (int)psb;
	double x0b = psb - (double)iqb;
	double x1b = 1.0 - x0b;
	double x2b = 2.0 - x0b;
	double x3b = 3.0 - x0b;

	double x123b = x1b * x2b * x3b / 6.0;
	double x120b = x1b * x2b * x0b / 6.0;
	double x032b = x0b * x3b * x2b / 2.0;
	double x031b = x0b * x3b * x1b / 2.0;

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

	if (distance10 < tiny1)
	{
		distance10 += tiny1;
	}
	if (distance20 < tiny1)
	{
		distance20 += tiny1;
	}

	// Find three dimension of 'Table_DSR' '
    const int Tpair1 = T1;
    const int Tpair2 = T2;
	const int T1_2Lplus1 = talpha.DS_2Lplus1[T1];
    const int T2_2Lplus1 = talpha.DS_2Lplus1[T2];

	//gaunt index
	const int gindex1 = L1 * L1 + m1;
	const int gindex2 = L2 * L2 + m2;

	// Peize Lin change rlya, rlyb, grlyb 2016-08-26
	vector<double> rlya;
	vector<double> rlyb;
	vector<vector<double>> grlyb;

	ModuleBase::Ylm::rl_sph_harm(T1_2Lplus1 - 1, dRa.x, dRa.y, dRa.z, rlya);
	if (job == 0)
	{
		ModuleBase::Ylm::rl_sph_harm(T2_2Lplus1 - 1, dRb.x, dRb.y, dRb.z, rlyb);
	}
	else
	{
		ModuleBase::Ylm::grad_rl_sph_harm(T2_2Lplus1 - 1, dRb.x, dRb.y, dRb.z, rlyb, grlyb);
	}
	//////////////////////////////////////////////////////////////////////////
	/// Formula :                         T1       T0          T0        T2
	/// \f[
	///	\sum_{ L0 }sum_{ m0 }
	/// 			D_{L0,L0} <\psi1_{L1,N1}|\alpha_{L0,m0}><\alpha _{L0,m0}|\psi2_{L2,N2}>
	///\f]
	//////////////////////////////////////////////////////////////////////////

    int ip = -1;
    int nb = 0; //for L0, N0

    for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
    {
        for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
        {
            if (!calproj[nb])
            {
                continue;
            }
            ++nb;
            
            // <psi1 | Beta>
            const int Opair1 = talpha.DS_Opair(Tpair1, L1, L0, N1, N0);
            // <psi2 | Beta>
            const int Opair2 = talpha.DS_Opair(Tpair2, L2, L0, N2, N0);
            const int inl = inl_index[T0](A0, L0, N0);
            for (int m01 = 0;m01 < 2 * L0 + 1;++m01)
            {
                for (int m02 = 0; m02 < 2 * L0 + 1; ++m02)
                {
                    ++ip;   //radial*angular
                    int gindex01 = L0 * L0 + m01;
                    int gindex02 = L0 * L0 + m02;
                    

                    //loop of {lmn}
                    double term_a = 0.0;
                    double term_b = 0.0;
                    double term_c[3] = {0, 0, 0};

                    //=============
                    // FIRST PART
                    //=============
                    for (int L = 0; L < T1_2Lplus1; L++)
                    {
                        //triangle rule for gaunt coefficients
                        int AL = L1 + L0;
                        int SL = abs(L1 - L0);
                        if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
                        {
                            continue;
                        }

                        //prefac = (i)^{lphi - lbeta - l}
                        //R0-R1 ==> <phi|beta>
                        double i_exp = pow(-1.0, (L1 - L0 - L) / 2);
                        double rl1 = pow(distance10, L);
                        double Interp_Vnla = 0.0;
                        if (distance10 > tiny2)
                        {
                            curr = talpha.Table_DSR[0][Tpair1][Opair1][L];
                            if (iqa >= rmesh1[nb] - 4)
                            {
                                Interp_Vnla = 0.0;
                            }
                            else
                            {
                                Interp_Vnla = i_exp * (x123a * curr[iqa] 
                                + x120a * curr[iqa + 3] 
                                + x032a * curr[iqa + 1] 
                                - x031a * curr[iqa + 2]);
                            }
                            Interp_Vnla /= rl1;
                        }
                        else
                        {
                            Interp_Vnla = i_exp * talpha.Table_DSR[0][Tpair1][Opair1][L][0];
                        }

                        /////////////////////////////////////
                        ///  Overlap value = S_from_table * G * Ylm
                        ////////////////////////////////////
                        for (int m = 0; m < 2 * L + 1; m++)
                        {
                            int gindexa = L * L + m;
                            //double tmpGaunt = this->MGT.Get_Gaunt_SH(L1, m1, L0, m0, L, m);
                            double tmpGaunt = this->MGT.Gaunt_Coefficients(gindex1, gindex01, gindexa);
                            term_a += tmpGaunt * Interp_Vnla * rlya[MGT.get_lm_index(L, m)];
                        }
                    } //end L

                    //=============
                    // SECOND PART
                    //=============
                    for (int L = 0; L < T2_2Lplus1; L++)
                    {
                        //triangle rule for gaunt coefficients
                        int AL = L2 + L0;
                        int SL = abs(L2 - L0);
                        if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
                        {
                            continue;
                        }

                        double Interp_Vnlb = 0.0;
                        double Interp_Vnlc = 0.0;

                        //prefac
                        double i_exp = pow(-1.0, (L2 - L0 - L) / 2);
                        double rl2 = pow(distance20, L);

                        if (distance20 > tiny2)
                        {
                            curr = talpha.Table_DSR[0][Tpair2][Opair2][L];

                            if (iqb >= rmesh2[nb] - 4)
                            {
                                Interp_Vnlb = 0.0;
                            }
                            else
                            {
                                Interp_Vnlb = i_exp * (x123b * curr[iqb] 
                                + x120b * curr[iqb + 3] 
                                + x032b * curr[iqb + 1] 
                                - curr[iqb + 2] * x031b);
                            }

                            Interp_Vnlb /= rl2;
                        }
                        else
                        {
                            Interp_Vnlb = i_exp * talpha.Table_DSR[0][Tpair2][Opair2][L][0];
                        }// end if(distance20)

                        if (job == 1) // 1 means calculate the derivative part.
                        {
                            if (distance20 > tiny2)
                            {
                                curr = talpha.Table_DSR[1][Tpair2][Opair2][L];

                                if (iqb >= rmesh2[nb] - 4)
                                {
                                    Interp_Vnlc = 0.0;
                                }
                                else
                                {
                                    Interp_Vnlc = i_exp * (x123b * curr[iqb] 
                                    + x120b * curr[iqb + 3] 
                                    + x032b * curr[iqb + 1] 
                                    - curr[iqb + 2] * x031b);
                                }
                                Interp_Vnlc = Interp_Vnlc / pow(distance20, L) - Interp_Vnlb * L / distance20;
                            }
                            else
                            {
                                Interp_Vnlc = 0.0;
                            }
                        } // end job==1

                        // sum up the second part.
                        for (int m = 0; m < 2 * L + 1; m++)
                        {
                            int gindexb = L * L + m;
                            //double tmpGaunt = this->MGT.Get_Gaunt_SH(L0, m0, L2, m2, L, m);
                            double tmpGaunt = this->MGT.Gaunt_Coefficients(gindex02, gindex2, gindexb);
                            const int lm = MGT.get_lm_index(L, m);

                            switch (job)
                            {
                                case 0: // calculate the overlap part.
                                {
                                    term_b += tmpGaunt * Interp_Vnlb * rlyb[lm];
                                    break;
                                }
                                case 1: // calculate the derivative part.
                                {
                                    double tt1 = tmpGaunt * Interp_Vnlc * rlyb[lm] / distance20;
                                    double tt2 = tmpGaunt * Interp_Vnlb;

                                    for (int ir = 0; ir < 3; ir++)
                                    {
                                        term_c[ir] += tt1 * unit_vec_dRb[ir] + tt2 * grlyb[lm][ir];
                                    }

                                    break;
                                }
                                default:
                                    break;
                            }
                        } // end m of SECOND PART
                    } // end L of SECOND PART

                    //===============================================
                    // THIRD PART: SUM THE VALUE FROM ALL PROJECTS.
                    //===============================================
                    const int nm = 2 * L0 + 1;
                    switch (job)
                    {
                        case 0: //calculate the overlap part.
                        {
                            nlm[0] += term_a * term_b * gedm[inl][m01*nm+m02]; 
                            break;
                        }
                        case 1: //calculate the derivative part.
                        {
                            for (int jr = 0; jr < 3; jr++)
                            {
                                nlm[jr] += term_c[jr] * term_a * gedm[inl][m01*nm+m02];
                            }
                            break;
                        }
                        default:
                            break;
                    }
                } //end m02
            }// end m01
        }//end N0
	}// end L0

	delete[] calproj;
	delete[] rmesh1;
	delete[] rmesh2;

	ModuleBase::timer::tick("ORB_gen_tables", "snap_psialpha");
	return;
}
#endif

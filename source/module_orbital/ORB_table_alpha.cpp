//caoyu add 2021-03-17
#include "ORB_table_alpha.h"
#include "ORB_read.h"
#include "../module_base/math_integral.h"
#include <stdexcept>
#include "../module_base/memory.h"
#include "../module_base/timer.h"

double ORB_table_alpha::dr = -1.0;

ORB_table_alpha::ORB_table_alpha()
{
	destroy_nr = false;

	ntype = 0;
	lmax = 0;
	kmesh = 0;
	Rmax = 0.0;
	dr = 0.0;
	dk = 0.0;

	nlm = 0;
	Rmesh = 0;

	kpoint = new double[1];
	r = new double[1];
	rab = new double[1];
	kab = new double[1];
	DS_2Lplus1 = new int[1];
}

ORB_table_alpha::~ORB_table_alpha()
{
	delete[] kpoint;
	delete[] r;
	delete[] rab;
	delete[] kab;
	delete[] DS_2Lplus1;
}

void ORB_table_alpha::allocate(
	const int &ntype_in,
	const int &lmax_in,
	const int &kmesh_in,
	const double &Rmax_in,
	const double &dr_in,
	const double &dk_in)
{
	ModuleBase::TITLE("ORB_table_alpha", "allocate");

	this->ntype = ntype_in; // type of elements.
	this->lmax = lmax_in;
	this->kmesh = kmesh_in;
	this->Rmax = Rmax_in;
	this->dr = dr_in;
	this->dk = dk_in;

	assert(ntype > 0);
	assert(lmax >= 0);
	assert(kmesh > 0.0);
	assert(Rmax >= 0.0);
	assert(dr > 0.0);
	assert(dk > 0.0);

	// calculated from input parameters
	this->nlm = (2 * lmax + 1) * (2 * lmax + 1);
	this->Rmesh = static_cast<int>(Rmax / dr) + 4;
	if (Rmesh % 2 == 0)
	{
		++Rmesh;
	}

	//	OUT(GlobalV::ofs_running,"lmax",lmax);
	//	OUT(GlobalV::ofs_running,"Rmax (Bohr)",Rmax);
	//	OUT(GlobalV::ofs_running,"dr (Bohr)",dr);
	//	OUT(GlobalV::ofs_running,"dk",dk);
	//	OUT(GlobalV::ofs_running,"nlm",nlm);
	//	OUT(GlobalV::ofs_running,"kmesh",kmesh);

	delete[] kpoint;
	delete[] r;
	kpoint = new double[kmesh];
	r = new double[Rmesh];

	delete[] rab;
	delete[] kab;
	kab = new double[kmesh];
	rab = new double[Rmesh];

	for (int ik = 0; ik < kmesh; ik++)
	{
		kpoint[ik] = ik * dk_in;
		kab[ik] = dk_in;
	}

	for (int ir = 0; ir < Rmesh; ir++)
	{
		r[ir] = ir * dr;
		rab[ir] = dr;
	}

	//	OUT(GlobalV::ofs_running,"allocate kpoint, r, rab, kab","Done");
	return;
}

int ORB_table_alpha::get_rmesh(const double &R1, const double &R2)
{
	int rmesh = static_cast<int>((R1 + R2) / ORB_table_alpha::dr) + 5;

	//mohan update 2009-09-08 +1 ==> +5
	//considering interpolation or so on...
	if (rmesh % 2 == 0)
	{
		rmesh++;
	}

	if (rmesh <= 0)
	{
//		GlobalV::ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
//		GlobalV::ofs_warning << "\n rmesh = " << rmesh;
		std::cout << "\n R1 = " << R1 << " R2 = " << R2;
		std::cout << "\n rmesh = " << rmesh;
		ModuleBase::WARNING_QUIT("ORB_table_alpha::get_rmesh", "rmesh <= 0");
	}
	return rmesh;
}

#include "../module_base/mathzone_add1.h"
#include "../module_base/constants.h"

void ORB_table_alpha::cal_S_PhiAlpha_R(
	ModuleBase::Sph_Bessel_Recursive::D2 *pSB, // mohan add 2021-03-06
	const int &l,
	const Numerical_Orbital_Lm &n1,
	const Numerical_Orbital_Lm &n2,
	const int &rmesh,
	double *rs,
	double *drs)
{
	ModuleBase::timer::tick("ORB_table_alpha", "S_PhiAlpha_R");

	assert(kmesh > 0);

	//start calc
	double *k1_dot_k2 = new double[kmesh];

	for (int ik = 0; ik < kmesh; ik++)
	{
		k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
	}

	//previous version
	double *integrated_func = new double[kmesh];

	const vector<vector<double>> &jlm1 = pSB->get_jlx()[l - 1];
	const vector<vector<double>> &jl = pSB->get_jlx()[l];
	const vector<vector<double>> &jlp1 = pSB->get_jlx()[l + 1];

	for (int ir = 0; ir < rmesh; ir++)
	{
		ModuleBase::GlobalFunc::ZEROS(integrated_func, kmesh);
		double temp = 0.0;

		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jl[ir][ik] * k1_dot_k2[ik];
		}
		// Call simpson integration
		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func, kab, temp);
		rs[ir] = temp * ModuleBase::FOUR_PI;

		//drs
		double temp1, temp2;

		if (l > 0)
		{
			for (int ik = 0; ik < kmesh; ik++)
			{
				integrated_func[ik] = jlm1[ir][ik] * k1_dot_k2[ik] * kpoint[ik];
			}

			ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func, kab, temp1);
		}

		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jlp1[ir][ik] * k1_dot_k2[ik] * kpoint[ik];
		}

		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func, kab, temp2);

		if (l == 0)
		{
			drs[ir] = -ModuleBase::FOUR_PI * temp2;
		}
		else
		{
			drs[ir] = ModuleBase::FOUR_PI * (temp1 * l - (l + 1) * temp2) / (2.0 * l + 1);
		}
	}

	//liaochen modify on 2010/4/22
	//special case for R=0
	//we store Slm(R) / R**l at the fisrt point, rather than Slm(R)
	if (l > 0)
	{
		ModuleBase::GlobalFunc::ZEROS(integrated_func, kmesh);
		double temp = 0.0;

		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = k1_dot_k2[ik] * pow(kpoint[ik], l);
		}

		// Call simpson integration
		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func, kab, temp);
		rs[0] = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac(2 * l + 1) * temp;
	}

	delete[] integrated_func;
	delete[] k1_dot_k2;

	ModuleBase::timer::tick("ORB_table_alpha", "S_PhiAlpha_R");
	return;
}

void ORB_table_alpha::init_Table_Alpha(
	ModuleBase::Sph_Bessel_Recursive::D2 *pSB,
	LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("ORB_table_alpha", "init_Table_Alpha");
	ModuleBase::timer::tick("ORB_table_alpha", "init_Table_Alpha");

	assert(ntype > 0);

	// (1) allocate 1st dimension ( overlap, derivative)
	this->Table_DSR = new double ****[2];
	// (2) allocate 2nd dimension ( overlap, derivative)
	this->Table_DSR[0] = new double ***[this->ntype];
	this->Table_DSR[1] = new double ***[this->ntype];

	// <1Phi|2Alpha>
	for (int T1 = 0; T1 < ntype; T1++) // type 1 is orbital
	{
		const int Lmax1 = orb.Phi[T1].getLmax();
		const int Lmax2 = orb.Alpha[0].getLmax();
		const int lmax_now = std::max(Lmax1, Lmax2);
		int L2plus1 = 2 * lmax_now + 1;
		//-------------------------------------------------------------
		// how many <psi|alpha_l>
		// here we count all possible psi with (L,N) index for type T1.
		//-------------------------------------------------------------
		const int pairs_chi = orb.Phi[T1].getTotal_nchi() * orb.Alpha[0].getTotal_nchi();

		if (pairs_chi == 0)
		{
			continue;
		}

		// init 2nd dimension
		this->Table_DSR[0][T1] = new double **[pairs_chi];
		this->Table_DSR[1][T1] = new double **[pairs_chi];

		const double Rcut1 = orb.Phi[T1].getRcut();
		for (int L1 = 0; L1 < Lmax1 + 1; L1++)
		{
			for (int N1 = 0; N1 < orb.Phi[T1].getNchi(L1); N1++)
			{
				for (int L2 = 0; L2 < Lmax2 + 1; L2++)
				{
					for (int N2 = 0; N2 < orb.Alpha[0].getNchi(L2); N2++)
					{
						// get the second index.
						const int Opair = this->DS_Opair(T1, L1, L2, N1, N2);

						// init 3rd dimension
						this->Table_DSR[0][T1][Opair] = new double *[L2plus1];
						this->Table_DSR[1][T1][Opair] = new double *[L2plus1];

						const double Rcut1 = orb.Phi[T1].getRcut();
						const double Rcut2 = orb.Alpha[0].getRcut();
						assert(Rcut1 > 0.0 && Rcut1 < 100);
						assert(Rcut2 > 0.0 && Rcut2 < 100);

						const int rmesh = this->get_rmesh(Rcut1, Rcut2);
						assert(rmesh < this->Rmesh);

						//L=|L1-L2|,|L1-L2|+2,...,L1+L2
						const int SL = abs(L1 - L2);
						const int AL = L1 + L2;

						for (int L = 0; L < L2plus1; L++)
						{
							//Allocation
							this->Table_DSR[0][T1][Opair][L] = new double[rmesh];
							this->Table_DSR[1][T1][Opair][L] = new double[rmesh];

							ModuleBase::Memory::record("ORB_table_alpha", "Table_DSR",
										   2 * this->ntype * pairs_chi * rmesh, "double");

							//for those L whose Gaunt Coefficients = 0, we
							//assign every element in Table_DSR as zero
							if ((L > AL) || (L < SL) || ((L - SL) % 2 == 1))
							{
								ModuleBase::GlobalFunc::ZEROS(Table_DSR[0][T1][Opair][L], rmesh);
								ModuleBase::GlobalFunc::ZEROS(Table_DSR[1][T1][Opair][L], rmesh);

								continue;
							}

							this->cal_S_PhiAlpha_R(
								pSB, // mohan add 2021-03-06
								L,
								orb.Phi[T1].PhiLN(L1, N1),
								orb.Alpha[0].PhiLN(L2, N2), // mohan update 2011-03-07
								rmesh,
								this->Table_DSR[0][T1][Opair][L],
								this->Table_DSR[1][T1][Opair][L]);
						} // end L2plus1
					}	  // end N2
				}		  // end L2
			}			  // end N1
		}				  // end L1
	}					  // end T1
	destroy_nr = true;

	//	OUT(GlobalV::ofs_running,"allocate non-local potential matrix","Done");
	ModuleBase::timer::tick("ORB_table_alpha", "init_Table_Alpha");
	return;
}

void ORB_table_alpha::Destroy_Table_Alpha(LCAO_Orbitals &orb)
{
	if (!destroy_nr)
	{
		return;
	}

	const int ntype = orb.get_ntype();
	for (int ir = 0; ir < 2; ir++)
	{
		for (int T1 = 0; T1 < ntype; T1++)
		{
			const int Lmax1 = orb.Phi[T1].getLmax();
			const int Lmax2 = orb.Alpha[0].getLmax();
			const int lmax_now = std::max(Lmax1, Lmax2);
			const int pairs = orb.Phi[T1].getTotal_nchi() * orb.Alpha[0].getTotal_nchi();

			// mohan fix bug 2011-03-30
			if (pairs == 0)
			{
				continue;
			}

			for (int dim2 = 0; dim2 < pairs; dim2++)
			{
				for (int L = 0; L < 2 * lmax_now + 1; L++)
				{
					delete[] Table_DSR[ir][T1][dim2][L];
				}
				delete[] Table_DSR[ir][T1][dim2];
			}
			delete[] Table_DSR[ir][T1];
		}
		delete[] Table_DSR[ir];
	}
	delete[] Table_DSR;
	return;
}

void ORB_table_alpha::init_DS_2Lplus1(LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("Make_Overlap_Table", "init_DS_2Lplus1");
	assert(this->ntype > 0);
	delete[] DS_2Lplus1;
	DS_2Lplus1 = new int[ntype]; // 2Lmax+1 for each T1

	int index = 0;
	for (int T1 = 0; T1 < ntype; T1++)
	{
		this->DS_2Lplus1[T1] = max(orb.Phi[T1].getLmax(), orb.Alpha[0].getLmax()) * 2 + 1;
	}
	return;
}

void ORB_table_alpha::init_DS_Opair(LCAO_Orbitals &orb)
{
	const int lmax = orb.get_lmax();
	const int nchimax = orb.get_nchimax();
	const int lmax_d = orb.get_lmax_d();
	const int nchimax_d = orb.get_nchimax_d();
	assert(lmax + 1 > 0);
	assert(lmax_d + 1 > 0);
	assert(nchimax > 0);
	assert(nchimax_d > 0);

	this->DS_Opair.create(this->ntype, lmax + 1, lmax_d + 1, nchimax, nchimax_d);

	// <1psi|2beta>
	// 1. orbital
	for (int T1 = 0; T1 < ntype; T1++) //alpha is not related to atom type !
	{
		int index = 0;
		for (int L1 = 0; L1 < orb.Phi[T1].getLmax() + 1; L1++)
		{
			for (int N1 = 0; N1 < orb.Phi[T1].getNchi(L1); N1++)
			{
				for (int L2 = 0; L2 < orb.Alpha[0].getLmax() + 1; L2++)
				{
					for (int N2 = 0; N2 < orb.Alpha[0].getNchi(L2); N2++)
					{
						this->DS_Opair(T1, L1, L2, N1, N2) = index;
						++index;
					}
				}
			}
		}
	}
	return;
}


//caoyu add 2021-03-20
void ORB_table_alpha::print_Table_DSR(LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("ORB_table_alpha", "print_Table_DSR");
	//NEW_PART("Overlap table S between lcao orbital and descriptor basis : S_{I_mu_alpha}");

	ofstream ofs;
	stringstream ss;
	// the parameter 'winput::spillage_outdir' is read from INPUTw.
	ss << "./S_I_mu_alpha.dat";
	if (GlobalV::MY_RANK == 0)
	{
		ofs.open(ss.str().c_str());
	}
	ofs << std::setprecision(10);

	for (int T1 = 0; T1 < this->ntype; T1++)	//T1
	{
		const int Lmax1 = orb.Phi[T1].getLmax();
		const int Lmax2 = orb.Alpha[0].getLmax();
		for (int L1 = 0; L1 < Lmax1 + 1; L1++)
		{
			for (int N1 = 0; N1 < orb.Phi[T1].getNchi(L1); N1++)
			{
				for (int L2 = 0; L2 < Lmax2 + 1; L2++)
				{
					for (int N2 = 0; N2 < orb.Alpha[0].getNchi(L2); N2++)
					{
						const int Opair = this->DS_Opair(T1, L1, L2, N1, N2);	//Opair
						//ofs <<setw(20)<< "atom_type: " << label << endl;
						ofs <<setw(20)<< "lcao basis: " << "L1=" << L1 << ", N1=" << N1 << endl;
						ofs <<setw(20)<< "descriptor basis: " << "L2=" << L2 << ", N2=" << N2 << endl;
						for (int il = 0; il < this-> DS_2Lplus1[T1]; il++)
						{
							ofs << "L=" << il << endl;
							const double Rcut1 = orb.Phi[T1].getRcut();
							const double Rcut2 = orb.Alpha[0].getRcut();
							const int rmesh = this->get_rmesh(Rcut1, Rcut2);

							if (Table_DSR[0][T1][Opair][il][1]==0)	//remain to be discussed
							{
								ofs << "S(R)=0"<<endl<<endl;
								continue;
							}
							ofs << "Rcut1="<<Rcut1<<", Rcut2="<<Rcut2<<", rmesh="<<rmesh<<", dr="<<this->dr<<";"<<endl;
							for (int ir = 0; ir < rmesh; ir++)
							{
								ofs << Table_DSR[0][T1][Opair][il][ir] << " ";
								if ( (ir+1) % 8 == 0) ofs << endl;
							}
							ofs << endl <<endl;
						}// il
					}// N2
				}// L2
			}// N1
		}// L1
	}// T1
	return;
}


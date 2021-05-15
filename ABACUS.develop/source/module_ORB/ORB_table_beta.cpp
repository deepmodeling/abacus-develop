#include <stdexcept>
#include "ORB_table_beta.h"
#include "ORB_read.h"
#include "../src_global/math_integral.h"

double ORB_table_beta::dr = -1.0;

ORB_table_beta::ORB_table_beta()
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
	r=new double[1];
	rab=new double[1];
	kab=new double[1];
}

ORB_table_beta::~ORB_table_beta()
{
	delete[] kpoint;
	delete[] r;
	delete[] rab;
	delete[] kab;
}

void ORB_table_beta::allocate
(
 	const int &ntype_in,
    const int &lmax_in,
    const int &kmesh_in,
	const double &Rmax_in,
    const double &dr_in,
    const double &dk_in
)
{
	TITLE("ORB_table_beta", "allocate");

	this->ntype = ntype_in;// type of elements.
	this->lmax = lmax_in;
	this->kmesh = kmesh_in;
	this->Rmax = Rmax_in;
	this->dr = dr_in;
	this->dk = dk_in;

	assert(ntype > 0);
	assert(lmax >= 0);
	assert(kmesh > 0.0);
	assert(Rmax >= 0.0);
	assert(dr>0.0);
	assert(dk>0.0);

	// calculated from input parameters
	this->nlm = (2*lmax+1) * (2*lmax+1);
	this->Rmesh = static_cast<int>( Rmax/dr ) + 4;
	if(Rmesh%2==0) 
	{
		++Rmesh;
	}

//	OUT(ofs_running,"lmax",lmax);
//	OUT(ofs_running,"Rmax (Bohr)",Rmax);
//	OUT(ofs_running,"dr (Bohr)",dr);
//	OUT(ofs_running,"dk",dk);
//	OUT(ofs_running,"nlm",nlm);
//	OUT(ofs_running,"kmesh",kmesh);
	
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

//	OUT(ofs_running,"allocate kpoint, r, rab, kab","Done");
	return;
}


int ORB_table_beta::get_rmesh(const double &R1, const double &R2)
{
	int rmesh = static_cast<int>((R1+R2)/ ORB_table_beta::dr) + 5;
	//mohan update 2009-09-08 +1 ==> +5
	//considering interpolation or so on...
	if (rmesh % 2 == 0) rmesh ++;
	
	if(rmesh <= 0)
	{
		//ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
		//ofs_warning << "\n rmesh = " << rmesh;
		cout << "\n R1 = " << R1 << " R2 = " << R2;
		cout << "\n rmesh = " << rmesh;
		WARNING_QUIT("ORB_table_beta::get_rmesh", "rmesh <= 0");
	}
	return rmesh;
}



void ORB_table_beta::cal_VNL_PhiBeta_R(
		Sph_Bessel_Recursive::D2 *pSB, // mohan add 2021-03-06
		const int &l,
		const Numerical_Orbital_Lm &n1,
		const Numerical_Nonlocal_Lm &n2,
		const int &rmesh,
		double *rs,
		double *drs)
{
	timer::tick ("ORB_table_beta", "VNL_PhiBeta_R");

	assert(kmesh > 0);

	//start calc	
    double *k1_dot_k2 = new double[kmesh];

	for (int ik = 0; ik < kmesh; ik++)
	{
		k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getBeta_k(ik);
	}

	//previous version
	double* integrated_func = new double[kmesh];
	
	const vector<vector<double>> &jlm1 = pSB->get_jlx()[l-1];
	const vector<vector<double>> &jl = pSB->get_jlx()[l];
	const vector<vector<double>> &jlp1 = pSB->get_jlx()[l+1];	
	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func,kmesh);
		double temp = 0.0;
		
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jl[ir][ik] * k1_dot_k2[ik];
		}
		// Call simpson integration
		Integral::Simpson_Integral(kmesh,integrated_func,kab,temp);
		rs[ir] = temp * FOUR_PI;
		
		//drs
		double temp1, temp2;
		
		if (l > 0)
		{
			for (int ik = 0; ik < kmesh; ik++)
			{
				integrated_func[ik] = jlm1[ir][ik] * k1_dot_k2[ik] * kpoint[ik];
			}

			Integral::Simpson_Integral(kmesh,integrated_func,kab,temp1);
		}
		
				
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jlp1[ir][ik] * k1_dot_k2[ik] * kpoint[ik];
		}
		
		Integral::Simpson_Integral(kmesh,integrated_func,kab,temp2);
		
		if (l == 0)
		{
			drs[ir] = -FOUR_PI*temp2;
		}
		else
		{
			drs[ir] = FOUR_PI*(temp1*l-(l+1)*temp2)/(2.0*l+1);
		}
	}
	
	//liaochen modify on 2010/4/22
	//special case for R=0
	//we store Slm(R) / R**l at the fisrt point, rather than Slm(R)
	if (l > 0)
	{
		ZEROS(integrated_func,kmesh);
		double temp = 0.0;
	
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = k1_dot_k2[ik] * pow (kpoint[ik], l);
		}
		
		// Call simpson integration
		Integral::Simpson_Integral(kmesh,integrated_func,kab,temp);
		rs[0] = FOUR_PI / Mathzone_Add1::dualfac (2*l+1) * temp;
	}
	
	delete [] integrated_func;
	delete[] k1_dot_k2;

	timer::tick ("ORB_table_beta", "VNL_PhiBeta_R");
	return;
}


void ORB_table_beta::init_Table_Beta(Sph_Bessel_Recursive::D2 *pSB)
{
	TITLE("ORB_table_beta", "init_Table_Beta");
	timer::tick("ORB_table_beta", "init_Table_Beta",'D');

	// (1) allocate 1st dimension ( overlap, derivative)
	this->Table_NR = new double****[2];
	// (2) allocate 2nd dimension ( overlap, derivative)
	this->Table_NR[0] = new double*** [this->NL_nTpairs];
	this->Table_NR[1] = new double*** [this->NL_nTpairs];
	
	// <1Phi|2Beta> 
	for (int T1 = 0;  T1 < ntype ; T1++) // type 1 is orbital
	{
		for (int T2 = 0 ; T2 < ntype ; T2++)// type 2 is non-local projector
		{
			// Tpair: type pair.
			const int Tpair=this->NL_Tpair(T1,T2);
			const int Lmax1 = ORB.Phi[T1].getLmax();			
			const int NBeta = ORB.nproj[T2];
			
			//-------------------------------------------------------------
			// how many <psi|beta_l>
			// here we count all possible psi with (L,N) index for type T1.
			//-------------------------------------------------------------
			const int pairs_chi = ORB.Phi[T1].getTotal_nchi() * NBeta;

			// CAUTION!!!
			// no matter nchi = 0 or NBeta = 0,
			// means the Tpair in this table is never used!
			if(pairs_chi == 0)continue;

			// init 2nd dimension
			this->Table_NR[0][Tpair] = new double** [ pairs_chi ];
			this->Table_NR[1][Tpair] = new double** [ pairs_chi ];

            const int T12_2Lplus1 = this->NL_L2plus1(T1,T2);

			const double Rcut1 = ORB.Phi[T1].getRcut();
			for (int L1 = 0; L1 < Lmax1 + 1; L1++)
            {
                for (int N1 = 0; N1 < ORB.Phi[T1].getNchi(L1); N1++)
				{
					// number of projectors.
					for (int nb = 0; nb < NBeta; nb ++)
					{
						//const int L2 = ORB.Beta[T2].getL_Beta(nb); // mohan delete the variable 2021-05-07
						const int L2 = ORB.Beta[T2].Proj[nb].getL(); // mohan add 2021-05-07

						const double Rcut2 = ORB.Beta[T2].Proj[nb].getRcut();

						const int Opair = this->NL_Opair(Tpair,L1,N1,nb);
						assert( Opair < pairs_chi );

						// init 3rd dimension
						this->Table_NR[0][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
						this->Table_NR[1][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
						
						const int rmesh = this->get_rmesh( Rcut1, Rcut2);
						assert( rmesh < this->Rmesh );

						//not all L in T12_2Lplus1 would function
						const int SL = abs(L1-L2);
						const int AL = L1+L2;
							
						for (int L=0; L < T12_2Lplus1 ; L++)
						{
							//Allocation
							this->Table_NR[0][Tpair][Opair][L] = new double[rmesh];
							this->Table_NR[1][Tpair][Opair][L] = new double[rmesh];

							Memory::record("ORB_table_beta","Table_NR",
							2*NL_nTpairs*pairs_chi*rmesh,"double");

							//for those L whose Gaunt Coefficients = 0, we
							//assign every element in Table_NR as zero
							if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) 
							{
								ZEROS (Table_NR[0][Tpair][Opair][L], rmesh);
								ZEROS (Table_NR[1][Tpair][Opair][L], rmesh);
								
								continue;
							}

							assert(nb < ORB.nproj[T2]);	

							this->cal_VNL_PhiBeta_R(
								pSB, // mohan add 2021-03-06
								L,
                                ORB.Phi[T1].PhiLN(L1,N1),
                                ORB.Beta[T2].Proj[nb], // mohan update 2011-03-07
                                rmesh,
								this->Table_NR[0][Tpair][Opair][L],
								this->Table_NR[1][Tpair][Opair][L]);
						}// end T12_2Lplus1
					}// end L2
				}// end N1
			}// end L1
		}// end T2
	}// end T1
	destroy_nr = true;


//	OUT(ofs_running,"allocate non-local potential matrix","Done");
	timer::tick("ORB_table_beta", "init_Table_Beta",'D');
	return;
}


void ORB_table_beta::Destroy_Table_Beta(LCAO_Orbitals &orb)
{
	if(!destroy_nr) return;

	const int ntype = orb.get_ntype();
	for(int ir = 0; ir < 2; ir ++)
	{
		for(int T1=0; T1<ntype; T1++)
		{
			for(int T2=0; T2<ntype; T2++)
			{
				const int Tpair = this->NL_Tpair(T1,T2); 
				const int L2plus1 = this->NL_L2plus1(T1,T2);
				const int pairs = orb.Phi[T1].getTotal_nchi() * orb.nproj[T2]; 

				// mohan fix bug 2011-03-30
				if(pairs ==0) continue;
				for(int dim2=0; dim2<pairs; dim2++)
				{
					for(int L=0; L<L2plus1; L++)
					{
						delete[] Table_NR[ir][Tpair][dim2][L];
					}
					delete[] Table_NR[ir][Tpair][dim2];
				}
				delete[] Table_NR[ir][Tpair];
			}
		}
		delete[] Table_NR[ir];
	}
	delete[] Table_NR;
	return;
}


void ORB_table_beta::init_NL_Tpair(void)
{
	TITLE("ORB_table_beta","init_NL_index");
	assert(ntype>0);
	this->NL_nTpairs = this->ntype * this->ntype;	
	this->NL_Tpair.create( this->ntype, this->ntype);
	this->NL_L2plus1.create( this->ntype, this->ntype); // mohan fix bug 2011-03-14

//	OUT(ofs_running,"Number of Nonlocal Pairs",NL_nTpairs);

	int index = 0;
	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		for (int T0 = 0 ; T0 < ntype ; T0++)
		{
			 this->NL_Tpair(T1,T0) = index;
			 ++index;

			 // the pair < psi | beta >
			 // be careful! This is not a symmetry matrix.
			 this->NL_L2plus1(T1,T0) = std::max(ORB.Phi[T1].getLmax(), ORB.Beta[T0].getLmax() )*2+1;
			 
			 // there are special situations:
			 // for example, two H atom without projector.
			 // if we use s orbital, 
			 // Phi.getLmax = 0,
			 // Beta.getLmax < 0, 
			 // so the value is 1.
			 // however, there are no projectors.
			 if(NL_L2plus1(T1,T0) <= 0)
			 {
				WARNING_QUIT("ORB_table_beta::init_paris_nonlocal_type","NL_L2plus1<=0");
			 }
		}
	}
	return;
}



void ORB_table_beta::init_NL_Opair(LCAO_Orbitals &orb)
{
	const int lmax = orb.get_lmax();
	const int nchimax = orb.get_nchimax();
	const int nprojmax = orb.nprojmax;
	
	// may have bug if we use all H!
	if( nprojmax == 0)
	{
		WARNING("ORB_table_beta","nproj for nonlocal pseudopotetials are zero, it must be all H atoms");
		return;
	}
	assert( NL_nTpairs > 0);
	
	this->NL_Opair.create( this->NL_nTpairs, lmax+1, nchimax, nprojmax);
	
	// <1psi|2beta>
	// 1. orbital
	for(int T1=0; T1<ntype; T1++)
	{
		// 2. NL projector
		for(int T0=0; T0<ntype; T0++)
		{
			const int nlpair = this->NL_Tpair(T1, T0);
			int index = 0;
			for(int L1=0; L1<orb.Phi[T1].getLmax()+1; L1++)
			{
				for(int N1=0; N1<orb.Phi[T1].getNchi(L1); N1++)
				{
					// notice !! T0 must be Beta( Nonlocal projector)
					// mohan update 2011-03-07
					for(int ip=0; ip<orb.nproj[T0]; ip++)
					{
						assert( nlpair < NL_nTpairs );
						assert( L1 < lmax+1 );
						assert( N1 < nchimax );
						assert( ip < nprojmax );
						this->NL_Opair(nlpair, L1, N1, ip) = index;
						++index;
					}
				}
			}
		}
	}

	return;
}

#include <stdexcept>
#include "ORB_table_phi.h"
#include "module_base/math_integral.h"
#include "module_base/memory.h"
#include "module_base/constants.h"
#include "module_base/timer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

ORB_table_phi::ORB_table_phi()
{
	overlap_table_allocated = false;
	kinetic_table_allocated = false;

	ntype = 0;
	lmax = 0;
	kmesh = 0;
	Rmax = 0.0;
	dr = -1.0;
	dk = 0.0;

	nlm = 0;
	Rmesh = 0;

	kpoint = nullptr;
	r=nullptr;
	rab=nullptr;
	kab=nullptr;

	Table_SR = nullptr;
	Table_TR = nullptr;
}

ORB_table_phi::~ORB_table_phi()
{
	delete[] kpoint;
	delete[] r;
	delete[] rab;
	delete[] kab;

	// pSB does not own memory

	_destroy_table();
}

void ORB_table_phi::allocate
(
 	const int &ntype_in,
    const int &lmax_in,
    const int &kmesh_in,
	const double &Rmax_in,
    const double &dr_in,
    const double &dk_in
)
{
	ModuleBase::TITLE("ORB_table_phi", "allocate");

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

	return;
}

int ORB_table_phi::get_rmesh(const double &R1, const double &R2) const
{
	int rmesh = static_cast<int>((R1+R2)/ this->dr) + 5;
	//mohan update 2009-09-08 +1 ==> +5
	//considering interpolation or so on...
	if (rmesh % 2 == 0) rmesh ++;

	if(rmesh <= 0)
	{
		//GlobalV::ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
		//GlobalV::ofs_warning << "\n rmesh = " << rmesh;
		std::cout << "\n R1 = " << R1 << " R2 = " << R2;
		std::cout << "\n rmesh = " << rmesh;
		ModuleBase::WARNING_QUIT("ORB_table_phi::get_rmesh", "rmesh <= 0");
	}
	return rmesh;
}

#include "module_base/mathzone_add1.h"
// Peize Lin accelerate 2017-10-02
void ORB_table_phi::cal_ST_Phi12_R
(
 	const int &job,
    const int &l,
    const Numerical_Orbital_Lm &n1,
    const Numerical_Orbital_Lm &n2,
    const int &rmesh,
    double* rs,
	double* drs
) const
{
	ModuleBase::timer::tick("ORB_table_phi", "cal_ST_Phi12_R");

	std::vector<double> k1_dot_k2(kmesh);
	// Peize Lin change 2017-12-12
	switch(job)
	{
		case 1: // calculate overlap
			if( !n1.get_psif().empty() && !n2.get_psi_k2().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsif(ik) * n2.getPsi_k2(ik);
				}
			}
			else if( !n1.get_psi_k().empty() && !n2.get_psi_k().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
				}
			}
			else if( !n1.get_psi_k2().empty() && !n2.get_psif().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsif(ik);
				}
			}
			break;

		case 2: // calculate kinetic energy
			for (int ik = 0; ik < kmesh; ik++)
			{
				k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsi_k2(ik);
			}
			break;
	}

	std::vector<double> k1_dot_k2_dot_kpoint(kmesh);
	for (int ik = 0; ik < kmesh; ik++)
	{
		k1_dot_k2_dot_kpoint[ik] = k1_dot_k2[ik] * this->kpoint[ik];
	}

	//Drs
	//djl = (l*j(l-1) - (l+1)j(l+1))/(2l+1)

	//previous version

	//double* integrated_func = new double[kmesh];

	int ll =0;
	if(l==0) ll=0;
	else 	 ll=l-1;

	const std::vector<std::vector<double>> &jlm1 = pSB->get_jlx()[ll];
	const std::vector<std::vector<double>> &jl = pSB->get_jlx()[l];
	const std::vector<std::vector<double>> &jlp1 = pSB->get_jlx()[l+1];

#ifdef _OPENMP
	#pragma omp parallel for schedule(static)
#endif
	for (int ir = 0; ir < rmesh; ir++)
	{
		std::vector<double> integrated_func(kmesh);
		const std::vector<double> &jl_r = jl[ir];
		for (int ik=0; ik<kmesh; ++ik)
		{
			integrated_func[ik] = jl_r[ik] * k1_dot_k2[ik];
		}
		// Call simpson integration
		double temp = 0.0;

		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), dk, temp);
		rs[ir] = temp * ModuleBase::FOUR_PI ;

		// Peize Lin accelerate 2017-10-02
		const std::vector<double> &jlm1_r = jlm1[ir];
		const std::vector<double> &jlp1_r = jlp1[ir];
		const double fac = l/(l+1.0);
		if( l==0 )
		{
			for (int ik=0; ik<kmesh; ++ik)
			{
				integrated_func[ik] = jlp1_r[ik] * k1_dot_k2_dot_kpoint[ik];
			}
		}
		else
		{
			for (int ik=0; ik<kmesh; ++ik)
			{
				integrated_func[ik] = (jlp1_r[ik]-fac*jlm1_r[ik]) * k1_dot_k2_dot_kpoint[ik];
			}
		}

		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), dk, temp);
		drs[ir] = -ModuleBase::FOUR_PI*(l+1)/(2.0*l+1) * temp;
	}

	//liaochen modify on 2010/4/22
	//special case for R=0
	//we store Slm(R) / R**l at the fisrt point, rather than Slm(R)

	if (l > 0)
	{
		std::vector<double> integrated_func(kmesh);
		double temp = 0.0;

		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = k1_dot_k2[ik] * std::pow (kpoint[ik], l);
		}

		ModuleBase::Integral::Simpson_Integral(kmesh, integrated_func.data(), kab, temp);
		rs[0] = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac (2*l+1) * temp;
	}

	ModuleBase::timer::tick("ORB_table_phi", "cal_ST_Phi12_R");
}

#include "module_base/constants.h"

// Peize Lin add 2017-10-27
void ORB_table_phi::cal_ST_Phi12_R
(
 	const int &job,
    const int &l,
    const Numerical_Orbital_Lm &n1,
    const Numerical_Orbital_Lm &n2,
	const std::set<size_t> &radials,
    double* rs,
	double* drs
) const
{
//	ModuleBase::TITLE("ORB_table_phi","cal_ST_Phi12_R");
	ModuleBase::timer::tick("ORB_table_phi", "cal_ST_Phi12_R");

	std::vector<double> k1_dot_k2(kmesh);
	switch(job)
	{
		case 1: // calculate overlap
			if( !n1.get_psif().empty() && !n2.get_psi_k2().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsif(ik) * n2.getPsi_k2(ik);
				}
			}
			else if( !n1.get_psi_k().empty() && !n2.get_psi_k().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
				}
			}
			else if( !n1.get_psi_k2().empty() && !n2.get_psif().empty() )
			{
				for (int ik = 0; ik < kmesh; ik++)
				{
					k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsif(ik);
				}
			}
			break;

		case 2: // calculate kinetic energy
			for (int ik = 0; ik < kmesh; ik++)
			{
				k1_dot_k2[ik] = n1.getPsi_k2(ik) * n2.getPsi_k2(ik);
			}
			break;
	}

	std::vector<double> k1_dot_k2_dot_kpoint(kmesh);
	for (int ik = 0; ik < kmesh; ik++)
	{
		k1_dot_k2_dot_kpoint[ik] = k1_dot_k2[ik] * this->kpoint[ik];
	}

	std::vector<double> integrated_func(kmesh);

	const std::vector<std::vector<double>> &jlm1 = pSB->get_jlx()[l-1];
	const std::vector<std::vector<double>> &jl = pSB->get_jlx()[l];
	const std::vector<std::vector<double>> &jlp1 = pSB->get_jlx()[l+1];

	for( const size_t &ir : radials )
	{
		// if(rs[ir])  => rs[ir]  has been calculated
		// if(drs[ir]) => drs[ir] has been calculated
		// Actually, if(ir[ir]||dr[ir]) is enough. Double insurance for the sake of avoiding numerical errors
		if( rs[ir] && drs[ir] )	continue;

		const std::vector<double> &jl_r = jl[ir];
		for (int ik=0; ik<kmesh; ++ik)
		{
			integrated_func[ik] = jl_r[ik] * k1_dot_k2[ik];
		}
		double temp = 0.0;

		ModuleBase::Integral::Simpson_Integral(kmesh,ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func),dk,temp);
		rs[ir] = temp * ModuleBase::FOUR_PI ;

		const std::vector<double> &jlm1_r = jlm1[ir];
		const std::vector<double> &jlp1_r = jlp1[ir];
		const double fac = l/(l+1.0);
		if( l==0 )
		{
			for (int ik=0; ik<kmesh; ++ik)
			{
				integrated_func[ik] = jlp1_r[ik] * k1_dot_k2_dot_kpoint[ik];
			}
		}
		else
		{
			for (int ik=0; ik<kmesh; ++ik)
			{
				integrated_func[ik] = (jlp1_r[ik]-fac*jlm1_r[ik]) * k1_dot_k2_dot_kpoint[ik];
			}
		}

		ModuleBase::Integral::Simpson_Integral(kmesh,ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func),dk,temp);
		drs[ir] = -ModuleBase::FOUR_PI*(l+1)/(2.0*l+1) * temp;
	}

	// cal rs[0] special
	if (l > 0)
	{
		if( radials.find(0)!=radials.end() )
		{
			for (int ik = 0; ik < kmesh; ik++)
			{
				integrated_func[ik] = k1_dot_k2[ik] * pow (kpoint[ik], l);
			}
			double temp = 0.0;

			ModuleBase::Integral::Simpson_Integral(kmesh,ModuleBase::GlobalFunc::VECTOR_TO_PTR(integrated_func),dk,temp);

			// PLEASE try to make dualfac function as input parameters
			// mohan note 2021-03-23
			rs[0] = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac (2*l+1) * temp;
		}
	}

	ModuleBase::timer::tick("ORB_table_phi", "cal_ST_Phi12_R");

	return;
}



void ORB_table_phi::init_Table(LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("ORB_table_phi", "init_Table");
	ModuleBase::timer::tick("ORB_table_phi", "init_Table");
	const int ntype = orb.get_ntype();
	assert( this->dr > 0.0);
	assert( OV_nTpairs>0);

	// record necessary information for the sizes of tables
	nelem_ = ntype;
	lmax_.resize(nelem_);
	nchi_tot_.resize(nelem_);
	for (int ielem = 0; ielem != nelem_; ++ielem) {
		lmax_[ielem] = orb.Phi[ielem].getLmax();
		nchi_tot_[ielem] = orb.Phi[ielem].getTotal_nchi();
	}

	// init 1st dimension
	this->Table_SR = new double****[2];
	this->Table_TR = new double****[2];
	for(int ir = 0; ir < 2; ir++)
	{
		this->Table_SR[ir] = new double***[ this->OV_nTpairs ];
		this->Table_TR[ir] = new double***[ this->OV_nTpairs ];
	}

	size_t memory_cost = 0;
	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		// Notice !! T2 start from T1
		// means that T2 >= T1
		for (int T2 = T1 ; T2 < ntype ; T2++)
		{
			// get the bigger lmax between two types
			const int Tpair=this->OV_Tpair(T1,T2);
			const int Lmax1 = orb.Phi[T1].getLmax();
			const int Lmax2 = orb.Phi[T2].getLmax();

			//L2plus1 could be reduced by considering Gaunt Coefficient
			//remain to be modified
			//??????
			const int lmax_now = std::max( Lmax1, Lmax2 );


			///////////////////////////////////
			// mohan add 2011-03-07
			// I think the lmax_now should be judged from two
			// orbitals, not atom type!!!!!!!!!!!!!
			// there are space that can imporve the efficiency.
			//////////////////////////////////

			const int L2plus1 =  2*lmax_now + 1;

			const int nchi1 = orb.Phi[T1].getTotal_nchi();
			const int nchi2 = orb.Phi[T2].getTotal_nchi();
			const int pairs_chi = nchi1 * nchi2;

			// init 2nd dimension
			for(int ir = 0; ir < 2; ir++)
			{
				this->Table_SR[ir][ Tpair ] = new double**[pairs_chi];
				this->Table_TR[ir][ Tpair ] = new double**[pairs_chi];
			}

			const double Rcut1 = orb.Phi[T1].getRcut();
			const double Rcut2 = orb.Phi[T2].getRcut();
			assert(Rcut1>0.0 && Rcut1<100);
			assert(Rcut2>0.0 && Rcut2<100);

			const int rmesh = this->get_rmesh( Rcut1, Rcut2);
			assert( rmesh < this->Rmesh );
#ifdef __ORBITAL
			ModuleBase::GlobalFunc::MAKE_DIR("Table_SR0");
			ModuleBase::GlobalFunc::MAKE_DIR("Table_TR0");
#endif
			for (int L1 = 0; L1 < Lmax1 + 1; L1++)
			{
				for (int N1 = 0; N1 < orb.Phi[T1].getNchi(L1); N1++)
				{
					for (int L2 = 0; L2 < Lmax2 + 1; L2 ++)
					{
						for (int N2 = 0; N2 < orb.Phi[T2].getNchi(L2); N2++)
						{
							// get the second index.
							const int Opair = this->OV_Opair(Tpair,L1,L2,N1,N2);

							// init 3rd dimension
							for(int ir = 0; ir < 2; ir++)
							{
								this->Table_SR[ir][ Tpair ][ Opair ] = new double *[L2plus1];
								this->Table_TR[ir][ Tpair ][ Opair ] = new double *[L2plus1];
							}

							//L=|L1-L2|,|L1-L2|+2,...,L1+L2
							const int SL = std::abs(L1-L2);
							const int AL = L1+L2;

							for (int L=0; L < L2plus1 ; L++)
							{
								//Allocation
								Table_SR[0][Tpair][Opair][L] = new double[rmesh];
								Table_SR[1][Tpair][Opair][L] = new double[rmesh];
								Table_TR[0][Tpair][Opair][L] = new double[rmesh];
								Table_TR[1][Tpair][Opair][L] = new double[rmesh];

								memory_cost += rmesh * 4;

								//for those L whose Gaunt Coefficients = 0, we
								//assign every element in Table_SR or Table_TR as zero
								if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1))
								{
									ModuleBase::GlobalFunc::ZEROS (Table_SR[0][Tpair][Opair][L], rmesh);
									ModuleBase::GlobalFunc::ZEROS (Table_SR[1][Tpair][Opair][L], rmesh);
									ModuleBase::GlobalFunc::ZEROS (Table_TR[0][Tpair][Opair][L], rmesh);
									ModuleBase::GlobalFunc::ZEROS (Table_TR[1][Tpair][Opair][L], rmesh);

									continue;
								}

								this->cal_ST_Phi12_R(1,L,
										orb.Phi[T1].PhiLN(L1,N1),
										orb.Phi[T2].PhiLN(L2,N2),
										rmesh,
										Table_SR[0][Tpair][Opair][L],
										Table_SR[1][Tpair][Opair][L]);

								this->cal_ST_Phi12_R(2,L,
										orb.Phi[T1].PhiLN(L1,N1),
										orb.Phi[T2].PhiLN(L2,N2),
										rmesh,
										Table_TR[0][Tpair][Opair][L],
										Table_TR[1][Tpair][Opair][L]);

#ifdef __ORBITAL
								int plot_length = 20;

								std::stringstream ss_sr;
								ss_sr << "Table_SR0/"<<Tpair<<Opair<<L<<".dat";
								std::string filename1 = ss_sr.str();
								plot_table(filename1,plot_length,Table_SR[0][Tpair][Opair][L]);
								std::stringstream ss_tr;
								ss_tr << "Table_TR0/"<<Tpair<<Opair<<L<<".dat";
								std::string filename2 = ss_tr.str();
								plot_table(filename2,plot_length,Table_TR[0][Tpair][Opair][L]);
#endif
							}//end m
						}
					}//end jl
				}
			}// end il
		}// end jt
	}// end it

	overlap_table_allocated = true;
	kinetic_table_allocated = true;
	ModuleBase::Memory::record("ORB::Table_SR&TR", sizeof(double) * memory_cost);

	ModuleBase::timer::tick("ORB_table_phi", "init_Table");
	return;
}


void ORB_table_phi::Destroy_Table(LCAO_Orbitals &orb)
{
	if(!overlap_table_allocated && !kinetic_table_allocated) return;

	const int ntype = orb.get_ntype();
	int dim1 = 0;
	for (int ir = 0; ir < 2; ir++)
	{
    	for (int T1 = 0; T1 < ntype; T1++)
		{
			// Notice !! T2 start from T1
			// means that T2 >= T1
    	    for (int T2 = T1; T2 < ntype; T2++)
        	{
				const int Lmax1 = orb.Phi[T1].getLmax();
				const int Lmax2 = orb.Phi[T2].getLmax();
				const int lmax_now = std::max(Lmax1, Lmax2);
				const int pairs = orb.Phi[T1].getTotal_nchi() * orb.Phi[T2].getTotal_nchi();

				for (int dim2 = 0; dim2 < pairs; dim2++)
				{
					for (int L = 0; L < 2*lmax_now + 1; L++)
					{
						if(overlap_table_allocated) delete [] Table_SR[ir][dim1][dim2][L];
						if(kinetic_table_allocated) delete [] Table_TR[ir][dim1][dim2][L];
                	}
                	if(overlap_table_allocated) delete [] Table_SR[ir][dim1][dim2];
					if(kinetic_table_allocated) delete [] Table_TR[ir][dim1][dim2];
				}
            	if(overlap_table_allocated) delete [] Table_SR[ir][dim1];
				if(kinetic_table_allocated) delete [] Table_TR[ir][dim1];
            	dim1++;

			}
        }

		dim1 = 0;
		if(overlap_table_allocated) delete [] Table_SR[ir];
		if(kinetic_table_allocated) delete [] Table_TR[ir];
	}

	if(overlap_table_allocated) delete[] Table_SR;
	if(kinetic_table_allocated) delete[] Table_TR;

	Table_SR = nullptr;
	Table_TR = nullptr;

	overlap_table_allocated = false;
	kinetic_table_allocated = false;

	return;
}


void ORB_table_phi::_destroy_table() {
	if(!overlap_table_allocated && !kinetic_table_allocated) {
		return;
	}

	// below is almost the same as Destroy_Table
	int dim1 = 0;
	for (int ir = 0; ir < 2; ir++)
	{
		for (int T1 = 0; T1 < ntype; T1++)
		{
			// Notice !! T2 start from T1
			// means that T2 >= T1
			for (int T2 = T1; T2 < ntype; T2++)
			{
				//const int Lmax1 = orb.Phi[T1].getLmax();
				//const int Lmax2 = orb.Phi[T2].getLmax();
				const int Lmax1 = lmax_[T1];
				const int Lmax2 = lmax_[T2];

				const int lmax_now = std::max(Lmax1, Lmax2);

				//const int pairs = orb.Phi[T1].getTotal_nchi() * orb.Phi[T2].getTotal_nchi();
				const int pairs = nchi_tot_[T1] * nchi_tot_[T2];

				for (int dim2 = 0; dim2 < pairs; dim2++)
				{
					for (int L = 0; L < 2*lmax_now + 1; L++)
					{
						if(overlap_table_allocated) delete [] Table_SR[ir][dim1][dim2][L];
						if(kinetic_table_allocated) delete [] Table_TR[ir][dim1][dim2][L];
					}
					if(overlap_table_allocated) delete [] Table_SR[ir][dim1][dim2];
					if(kinetic_table_allocated) delete [] Table_TR[ir][dim1][dim2];
				}
				if(overlap_table_allocated) delete [] Table_SR[ir][dim1];
				if(kinetic_table_allocated) delete [] Table_TR[ir][dim1];
				dim1++;

			}
        }

		dim1 = 0;
		if(overlap_table_allocated) delete [] Table_SR[ir];
		if(kinetic_table_allocated) delete [] Table_TR[ir];
	}

	if(overlap_table_allocated) delete[] Table_SR;
	if(kinetic_table_allocated) delete[] Table_TR;

	overlap_table_allocated = false;
	kinetic_table_allocated = false;

	Table_SR = nullptr;
	Table_TR = nullptr;

	return;
}


void ORB_table_phi::init_OV_Tpair(LCAO_Orbitals &orb)
{
	ModuleBase::TITLE("ORB_table_phi","init_OV_Tpair");
    assert(ntype>0);

    this->OV_nTpairs = this->ntype * (this->ntype + 1) / 2;
    this->OV_Tpair.create(ntype, ntype);
	this->OV_L2plus1.create(ntype, ntype); // mohan fix bug 2011-03-14

    int index = 0;
    for (int T1 = 0;  T1 < ntype ; T1++)
    {
		// Notice !! T2 start from T1
		// means that T2 >= T1
        for (int T2 = T1 ; T2 < ntype ; T2++)
        {
			/// (1) pairs about atom types
			//liaochen modify 2010/8/4
			///index for T1>T2 is also needed
            this->OV_Tpair(T2, T1) = index;
			this->OV_Tpair(T1, T2) = this->OV_Tpair(T2, T1);

			++index;
			/// (2) pairs about lmax
			this->OV_L2plus1(T1,T2) = std::max(orb.Phi[T1].getLmax(), orb.Phi[T2].getLmax() )*2+1;
			this->OV_L2plus1(T2,T1) = this->OV_L2plus1(T1,T2);
        }
    }
    return;
}



void ORB_table_phi::init_OV_Opair(LCAO_Orbitals &orb)
{
    const int lmax = orb.get_lmax();
    const int nchimax = orb.get_nchimax();
	assert(lmax+1 > 0);
	assert(nchimax > 0);
	assert(OV_nTpairs > 0);

    this->OV_Opair.create(OV_nTpairs, lmax+1, lmax+1, nchimax, nchimax);

    for(int T1=0; T1<ntype; T1++)
    {
		// Notice !! T2 start from T1
		// means that T2 >= T1
        for(int T2=T1; T2<ntype; T2++)
        {
			const int dim1 = this->OV_Tpair(T1,T2);
			int index=0;
            for(int L1=0; L1<orb.Phi[T1].getLmax()+1; L1++)
            {
                for(int N1=0; N1<orb.Phi[T1].getNchi(L1); N1++)
                {
                    for(int L2=0; L2<orb.Phi[T2].getLmax()+1; L2++)
                    {
                        for(int N2=0; N2<orb.Phi[T2].getNchi(L2); N2++)
                        {
                            this->OV_Opair(dim1, L1, L2, N1, N2) = index;
                            ++index;
                        }// N2
                    }// L2
                }// N1
            }// L1
        }// T2
    }// T1
    return;
}

// Peize Lin update 2016-01-26
void ORB_table_phi::init_Lmax (
	const int orb_num,
	const int mode,
	int &Lmax_used,
	int &Lmax,
	const int &Lmax_exx,
	const LCAO_Orbitals &orb,
	const Numerical_Nonlocal* beta_) const
{

	auto cal_Lmax_Phi = [](int &Lmax,const LCAO_Orbitals &orb)
	{
		//obtain maxL of all type
		const int ntype = orb.get_ntype();
		for (int it = 0; it < ntype; it++)
		{
			Lmax = std::max(Lmax, orb.Phi[it].getLmax());
		}
	};

	auto cal_Lmax_Beta = [](int &Lmax,const LCAO_Orbitals &orb, const Numerical_Nonlocal* beta_)
	{
		// fix bug.
		// mohan add the nonlocal part.
		// 2011-03-07
		const int ntype = orb.get_ntype();
		for(int it=0; it< ntype; it++)
		{
			Lmax = std::max(Lmax, beta_[it].getLmax());
		}
	};
	auto cal_Lmax_Alpha = [](int &Lmax, const LCAO_Orbitals &orb)
	{
		//caoyu add 2021-08-05 for descriptor basis
		Lmax = std::max(Lmax, orb.get_lmax_d());
	};


	Lmax = -1;

	switch( orb_num )
	{
		case 2:
			switch( mode )
			{
				case 1:			// used in <Phi|Phi> or <Beta|Phi>
					cal_Lmax_Phi(Lmax,orb);
					cal_Lmax_Beta(Lmax,orb, beta_);
					//use 2lmax+1 in dS
					Lmax_used = 2*Lmax + 1;
					break;
				case 2:			// used in <jY|jY> or <Abfs|Abfs>
					Lmax = std::max(Lmax, Lmax_exx);
					Lmax_used = 2*Lmax + 1;
					break;
				case 3:                // used in berryphase by jingan
					cal_Lmax_Phi(Lmax,orb);
					Lmax++;
					Lmax_used = 2*Lmax + 1;
					break;
				default:
					throw std::invalid_argument("ORB_table_phi::init_Lmax orb_num=2, mode error");
					break;
			}
			break;
		case 3:
			switch( mode )
			{
				case 1:			// used in <jY|PhiPhi> or <Abfs|PhiPhi>
					cal_Lmax_Phi(Lmax,orb);
					Lmax_used = 2*Lmax + 1;
					Lmax = std::max(Lmax, Lmax_exx);
					Lmax_used += Lmax_exx;
					break;
				default:
					throw std::invalid_argument("ORB_table_phi::init_Lmax orb_num=3, mode error");
					break;
			}
			break;
		case 4:
			switch( mode )
			{
				case 1:			// used in <PhiPhi|PhiPhi>
					cal_Lmax_Phi(Lmax,orb);
					Lmax_used = 2*( 2*Lmax + 1 );
					break;
				default:
					throw std::invalid_argument("ORB_table_phi::init_Lmax orb_num=4, mode error");
					break;
			}
			break;
		default:
			throw std::invalid_argument("ORB_table_phi::init_Lmax orb_num error");
			break;
	}

	assert(Lmax_used >= 1);
}

// Peize Lin update 2016-01-26
void ORB_table_phi::init_Table_Spherical_Bessel (
	const int orb_num,
	const int mode,
	int &Lmax_used,
	int &Lmax,
	const int &Lmax_exx,
	const LCAO_Orbitals &orb,
	const Numerical_Nonlocal* beta_)
{
	ModuleBase::TITLE("ORB_table_phi", "init_Table_Spherical_Bessel");

	this->init_Lmax (orb_num,mode,Lmax_used,Lmax,Lmax_exx,orb, beta_);		// Peize Lin add 2016-01-26

	for( auto & sb : ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool )
	{
		if( this->dr * this->dk == sb.get_dx() )
		{
			pSB = &sb;
			break;
		}
	}

	if(!pSB)
	{
		ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.push_back({});
		pSB = &ModuleBase::Sph_Bessel_Recursive_Pool::D2::sb_pool.back();
	}

	pSB->set_dx( this->dr * this->dk );
	pSB->cal_jlx( Lmax_used, this->Rmesh, this->kmesh );

	ModuleBase::Memory::record ("ORB::Jl(x)", sizeof(double) * (Lmax_used+1) * this->kmesh * this->Rmesh);
}

void ORB_table_phi::plot_table(
	const std::string filename,
	const int rmesh,
	double* column)
{
	std::ofstream ofs;
	ofs.open(filename.c_str());
	ofs << "ir    table_entry" << std::endl;
	for(int ir=0;ir<rmesh;ir++)
	{
		ofs<< std::setw(4) << ir << "  " << column[ir]<< std::endl;
	}
}

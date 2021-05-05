//=========================================================
//AUTHOR : liaochen
//DATE : 2008-03-04
//=========================================================
#include "ORB_nonlocal.h"

Numerical_Nonlocal::Numerical_Nonlocal()
{
	//make pair of new and delete
	//question remains
	this->type = 0;
	this->lmax = 0;
	this->LfromBeta = new int[1];
	this->Proj = new Numerical_Nonlocal_Lm[1];
	this->nproj = -1;
	//zhengdy-soc, for optimize nonlocal part
	for(int is=0;is<4;is++) this->index1_soc[is] = new int[1];
	for(int is=0;is<4;is++) this->index2_soc[is] = new int[1];
}

Numerical_Nonlocal::~Numerical_Nonlocal()
{
	delete[] Proj;
	delete[] LfromBeta;
	for(int is=0;is<4;is++) delete[] this->index1_soc[is];
	for(int is=0;is<4;is++) delete[] this->index2_soc[is];
}

void Numerical_Nonlocal::set_type_info
(
	const int& type_in,
	const string& label_in,
	const string& type_ps_in,
	const int& lmax_in,
	matrix& Coefficient_D_in,
	ComplexMatrix& Coefficient_D_in_so,
	const int& nproj_in,
	const int& nproj_in_so,
	int* LfromBeta_in,
	const Numerical_Nonlocal_Lm* Proj_in,
	const bool has_so
)
{
	if (type_in < 0 || type_in > 2)
	{
		WARNING("Numerical_Nonlocal", "bad input of type_in: not ready yet for type >2");
	}

	this->type = type_in;
	this->label = label_in;
	this->type_ps = type_ps_in;

	if (lmax_in < -1 || lmax_in > 20)
	{
		WARNING_QUIT("Numerical_Nonlocal", "bad input of lmax : should be between -1 and 20");
	}

	this->lmax = lmax_in;
//----------------------------------------------------------
//EXPLAIN : Coefficient D used in calculate elements of NLps
//----------------------------------------------------------
/*2016-07-19, LiuXh
	this->Coefficient_D.create( lmax_in+1, lmax_in+1);
	for (int L1 = 0; L1 < lmax + 1; L1++)
	{
		for (int L2 = 0; L2 < lmax + 1; L2++)
		{
			this->Coefficient_D(L1, L2) = Coefficient_D_in(L1, L2);
		}
	}
2016-07-19, LiuXh*/

//----------------------------------------------------------
//EXPLAIN : LfromBeta
//----------------------------------------------------------
	this->nproj = nproj_in;
	if(has_so){ 
		this->nproj_soc = nproj_in_so;
	}
	//assert(nproj <= lmax_in+1); //LiuXh 2016-01-13, 2016-05-16
	assert(nproj <= nproj_in+1); //LiuXh 2016-01-13, 2016-05-16
	assert(nproj >= 0);

//2016-07-19 begin, LiuXh
	if(!has_so){
		this->Coefficient_D.create( nproj_in+1, nproj_in+1);
		ZEROS(this->non_zero_count_soc, 4);
		if(lmax_in > -1) //LiuXh add 20180328, fix bug of Hydrogen element with single projector pseudopot
		{ //LiuXh add 20180328
//			for (int L1 = 0; L1 < nproj + 1; L1++)
			for (int L1 = 0; L1 < min(this->Coefficient_D.nr, Coefficient_D_in.nr); L1++)
			{
//				for (int L2 = 0; L2 < nproj + 1; L2++)
				for (int L2 = 0; L2 < min(this->Coefficient_D.nc, Coefficient_D_in.nc); L2++)
				{
					this->Coefficient_D(L1, L2) = Coefficient_D_in(L1, L2);
				}
			}
		} //LiuXh add 20180328
	}
	else//zhengdy-soc
	{
		this->Coefficient_D_so.create(NSPIN,  nproj_soc+1,  nproj_soc+1);
		//optimize
		for(int is=0;is<4;is++)
		{
			this->non_zero_count_soc[is] = 0;
			delete[] this->index1_soc[is];
			this->index1_soc[is] = new int[nproj_soc * nproj_soc];
			delete[] this->index2_soc[is];
			this->index2_soc[is] = new int[nproj_soc * nproj_soc];
		}
		if(lmax_in > -1)
		{
			if(LSPINORB)
			{
				int is = 0;
				for (int is1 = 0; is1 < 2; is1++)
				{
					for (int is2 = 0; is2 < 2; is2++)
					{
						for (int L1 = 0; L1 < nproj_soc; L1++)
						{
							for (int L2 = 0; L2 < nproj_soc; L2++)
							{
								this->Coefficient_D_so(is, L1, L2) = Coefficient_D_in_so(L1 + nproj_soc*is1, L2 + nproj_soc*is2);
								if(fabs(this->Coefficient_D_so(is, L1, L2).real())>1.0e-8 || 
									fabs(this->Coefficient_D_so(is, L1, L2).imag())>1.0e-8 )
								{
									this->index1_soc[is][non_zero_count_soc[is]] = L1;
									this->index2_soc[is][non_zero_count_soc[is]] = L2;
									this->non_zero_count_soc[is]++;	
								}
							}
						}
						is++;
					}
				}
			}
			else
			{
				int is = 0;
				for (int is1 = 0; is1 < 2; is1++)
				{
					for (int is2 = 0; is2 < 2; is2++)
					{
						if(is>=NSPIN) break;
						for (int L1 = 0; L1 < nproj_soc; L1++)
						{
							for (int L2 = 0; L2 < nproj_soc; L2++)
							{
								if(is==1||is==2) this->Coefficient_D_so(is, L1, L2) = ZERO;
								else this->Coefficient_D_so(is, L1, L2) = Coefficient_D_in_so(L1 + nproj_soc*is1, L2 + nproj_soc*is2);
								if(abs(this->Coefficient_D_so(is, L1, L2).real())>1.0e-8 
									|| abs(this->Coefficient_D_so(is, L1, L2).imag())>1.0e-8)
								{
									this->index1_soc[is][non_zero_count_soc[is]] = L1;
									this->index2_soc[is][non_zero_count_soc[is]] = L2;
									this->non_zero_count_soc[is]++;	
								}
							}
						}
						is++;
					}
				}
				
			}
		}
	}
//2016-07-19 end, LiuXh

	delete[] LfromBeta;
	this->LfromBeta = new int[nproj];
	ZEROS(LfromBeta, nproj);
	for(int i=0; i<nproj; i++)
	{
		this->LfromBeta[i] = LfromBeta_in[i];
	}

//----------------------------------------------------------
// EXPLAIN : non_local pseudopotential projector for each l
//----------------------------------------------------------
	//only store radial function
	delete[] Proj;
	this->Proj = new Numerical_Nonlocal_Lm[this->nproj];

	for (int p1=0; p1<nproj; p1++)
	{
		this->Proj[p1] = Proj_in[p1];
	}

	this->rcut_max = 0.0;
	for(int p1=0; p1<nproj; p1++)
	{
		this->rcut_max = max( this->Proj[p1].getRcut(), rcut_max ); 
	}
	return;
}


#include <math.h>
#include <cassert>
#include "ORB_gaunt_table.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#include "module_base/mathzone.h"
#include "module_base/global_function.h"
#include "module_base/vector3.h"
#include "module_base/constants.h"
#include "module_base/math_ylmreal.h"

ORB_gaunt_table::ORB_gaunt_table(){}
ORB_gaunt_table::~ORB_gaunt_table(){}

void ORB_gaunt_table::init_Gaunt(const int &lmax)
{
////////////////////////////////////////
/// EXPLAIN : make table of Gaunt Coefficients
////////////////////////////////////////
    ModuleBase::TITLE("ORB_gaunt_table", "init_Gaunt");
    ModuleBase::timer::tick("ORB_gaunt_table", "init_Gaunt");
    
	const int nlm = (lmax * 2 + 1) * (lmax * 2 + 1);
	this->Gaunt_Coefficients.create(nlm, nlm, nlm);

	// Becaful! , L ends at 2*lmax+1
    for (int L = 0; L < 2*lmax + 1; L++)
    {
		// m order is ( 0,0,1,-1,0,1,-1,2,-2...)
        for (int m = 0; m < 2*L + 1 ; m++)
        {
            const int dim = this->get_lm_index(L,m);
            for (int L1 = 0; L1 < lmax + 1; L1++)
            {
                for (int m1 = 0; m1 < 2*L1+1 ; m1++)
                {
                    const int dim1 = this->get_lm_index(L1,m1);
                    for (int L2 = 0; L2 < lmax + 1; L2++)
                    {
                        for (int m2 = 0; m2 < 2*L2+1 ; m2++)
                        {
                            const int dim2 = this->get_lm_index(L2,m2);
                            /////////////////////
                            ///call Calculate::Cal_G
                            ////////////////////
					
							Gaunt_Coefficients(dim1, dim2, dim) = 
								this->Get_Gaunt_SH (L1, m1, L2, m2, L, m);	
						}// m2
                    }// L2
                }// m1
            }// L1
        }// m
    }// L

    ModuleBase::timer::tick("ORB_gaunt_table", "init_Gaunt");
    return;
}

/*
double ORB_gaunt_table::Cal_Gaunt_single
(
    const int &L1, 
	const int &m1,
    const int &L2, 
	const int &m2,
    const int &L, 
	const int &m,
    const double &s1, 
	const double &e1,
    const double &s2, 
	const double &e2
)
{
	ModuleBase::timer::tick("ORB_gaunt_table", "Cal_Gaunt_single");
	if ((L1 - L2 - L) % 2 != 0)
	{
		return 0.0;
	}

	double result = 0.0;
	static double absc[16] = { 
		-0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.755404408355003, 
		-0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744, 
		0.09501250983763744, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 
		0.755404408355003, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499 };

	static double weight[16] = { 
		0.02715245941175406, 0.06225352393864778, 0.0951585116824929, 0.1246289712555339, 
		0.1495959888165768, 0.1691565193950026, 0.1826034150449236, 0.1894506104550685, 
		0.1894506104550685, 0.1826034150449236, 0.1691565193950026, 0.1495959888165768, 
		0.1246289712555339, 0.0951585116824929, 0.06225352393864778, 0.02715245941175406 };

	for (int i = 0;i < 16;i++)
	{
		for (int j = 0;j < 16;j++)
		{
			double theta = ((s1 + e1) + (e1 - s1) * absc[i]) / 2;

			result += weight[i] * weight[j] * sin(theta) *
			          this->Ylm_Gaunt( this->get_lm_index(L1, m1), 16 * i + j) *
			          this->Ylm_Gaunt( this->get_lm_index(L2, m2), 16 * i + j) *
			          this->Ylm_Gaunt( this->get_lm_index(L,  m), 16 * i + j);
		}
	}

	result *= ((e1 - s1) / 2) * ((e2 - s2) / 2);
	ModuleBase::timer::tick("ORB_gaunt_table", "Cal_Gaunt_single");
	return result;
}
*/
/*
void ORB_gaunt_table::init_Ylm_Gaunt
(
 	const int &lmax,
    const double &s1, 
	const double &e1,
    const double &s2, 
	const double &e2
)
{
	ModuleBase::TITLE("ORB_gaunt_table", "init_Ylm_Gaunt");
	ModuleBase::timer::tick("ORB_gaunt_table", "inite_Ylm_Gaunt");

	const int nlm = (2*lmax+1) * (2*lmax+1);

	static double absc[16] = { 
		-0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.755404408355003, 
		-0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.09501250983763744, 
		0.09501250983763744, 0.2816035507792589, 0.4580167776572274, 0.6178762444026438, 
		0.755404408355003, 0.8656312023878318, 0.9445750230732326, 0.9894009349916499 };

	//initialization of ylm_map

	ModuleBase::Vector3<double> g_gaunt[256];

	this->Ylm_Gaunt.create(nlm , 256);

	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			const double theta = ((s1 + e1) + (e1 - s1) * absc[i]) / 2;
			const double phi = ((s2 + e2) + (e2 - s2) * absc[j]) / 2;
			ModuleBase::Vector3<double> u(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
			g_gaunt[16*i+j] = u;
		}
	}

	ModuleBase::YlmReal::Ylm_Real(nlm, 256, &g_gaunt[0], this->Ylm_Gaunt);

	ModuleBase::timer::tick("ORB_gaunt_table", "init_Ylm_Gaunt");
	return;
}
*/

int ORB_gaunt_table::get_lm_index(
	const int l, 
	const int m)
{
	return l*l+m;
}


/**********************//**
  * Rasch and Yu's Method
  ***********************/
///total pointers
//int ORB_gaunt_table::P_EL(const int& L)
//{
//	return (L+1) * (L+2) * (L+3) * (L+4) / 24;
//}


///effective pointers
int ORB_gaunt_table::EP_EL(const int& L)
{
	if(L % 2 == 0) return (L+2) * (L+4) * (3*L*L+14*L+24) / 192;
	else return (L+1) * (L+3) * (L+5) * (3*L+5) / 192;
}


int ORB_gaunt_table::index_func
(
 	const int& l1,
	const int& l2,
	const int& l3,
	const int& m3
)
{
	const int aux1 = l1*(l1*l1*l1+6*l1*l1+11*l1+6)/24;
	const int aux2 = l2*(l2*l2+3*l2+2)/6;
	const int aux3 = l3*(l3+1)/2;
	
	return aux1 + aux2 + aux3 + m3;
}


void ORB_gaunt_table::init_Gaunt_CH(const int& Lmax)
{
	ModuleBase::TITLE("ORB_gaunt_table","init_Gaunt_CH");
	ModuleBase::timer::tick("ORB_gaunt_table","init_Gaunt_CH");

//	assert(Lmax <= 6);			// Peize Lin delete 2016-08-26. why?

	int L = 2*Lmax + 1;
	
	//int Np = this->P_EL(L);
//	assert(Np <= 5000);			// Peize Lin delete 2016-08-26. why?

	int Eff_Np = this->EP_EL(L);

	ModuleBase::Memory::record("ORB::Gaunt_CH", sizeof(double) * Eff_Np * 30);
	
	int ic1 = 0;
	for(int l1 = 0; l1 <= L; l1++)
	{
		for(int l2 = 0; l2 <= l1; l2++)
		{
			for(int l3 = 0; l3 <= l2; l3++)
			{
				for(int m3 = 0; m3 <= l3; m3++)
				{
					int idx = index_func(l1, l2, l3, m3);
					assert(ic1 == idx);
					
					int l_sum = l1 + l2 + l3;
					if((l_sum % 2 == 0) && (l2 + l3 >= l1))
					{
						int uplmt_m2 = l1 - m3 > l2 ? l2 : l1 - m3;
						//int dim = l2 + uplmt_m2 + 1;
						
//						assert(dim <= 30);			// Peize Lin delete 2016-08-26. why?

						int ic2 = 0;
						for(int m2 = -l2; m2 <= uplmt_m2; m2++)
						{
							//m1 + m2 + m3 == 0
							int m1 = -m2 - m3;
							assert(std::abs(m1) <= l1);

							// Peize Lin delete assert 2016-08-26
//							assert(ic1 < 5000);
//							assert(ic2 < 30);
							Gaunt_CH[ic1][ic2] = Calc_Gaunt_CH(l1, m1, l2, m2, l3, m3);
							ic2++;
						}
					}

					ic1++;
				}// m3
			}// l3
		}// l2
	} // l1

	ModuleBase::timer::tick("ORB_gaunt_table","init_Gaunt_CH");
	return;
}


//using wigner 3j expression
double ORB_gaunt_table::Calc_Gaunt_CH
(
 	const int& l1,
	const int& m1,
	const int& l2, 
	const int& m2,
	const int& l3,
	const int& m3
)
{
//	ModuleBase::TITLE("ORB_gaunt_table","Calc_Gaunt_CH");
	ModuleBase::timer::tick("ORB_gaunt_table","Calc_Gaunt_CH");
	
	double fac = sqrt((2*l1+1) * (2*l2+1) * (2*l3+1) / ModuleBase::FOUR_PI);

	int g = (l1+l2+l3)/2;
	double triangle_f = sqrt( Fact(l1+l2-l3) * Fact(l1-l2+l3) * Fact(-l1+l2+l3) / Fact(2*g+1) );

	fac *= pow(-1.0, g) * triangle_f * Fact(g) / Fact(g-l1) / Fact(g-l2) / Fact(g-l3);

	double aux1 = sqrt(Fact(l1+m1) * Fact(l1-m1) * Fact(l2+m2) * Fact(l2-m2) * Fact(l3+m3) * Fact(l3-m3));

	int kmin, kmax;
	
	kmin = (l2-l3-m1) > (l1-l3+m2) ? (l2-l3-m1) : (l1-l3+m2);
	kmin = kmin > 0 ? kmin : 0;

	kmax = (l1+l2-l3) > (l1-m1) ? (l1-m1) : (l1+l2-l3);
	kmax = kmax > (l2+m2) ? (l2+m2) : kmax;

	double aux2 = 0.0;
	for(int k = kmin; k <= kmax; k++)
	{
		aux2 += pow(-1.0, k) / Fact(k) / Fact(l1+l2-l3-k) / Fact(l1-m1-k) / Fact(l2+m2-k)
					/ Fact(l3-l2+m1+k) / Fact(l3-l1+k-m2);
	}

	return fac * pow(-1.0, l1-l2-m3) * triangle_f * aux1 * aux2;

	ModuleBase::timer::tick("ORB_gaunt_table","Calc_Gaunt_CH");
}
	

double ORB_gaunt_table::Get_Gaunt_CH
(
 	const int& l1,
	const int& m1,
	const int& l2,
	const int& m2,
	const int& l3,
	const int& m3
)
{
	assert(l1 >= 0);
	assert(l2 >= 0);
	assert(l3 >= 0);
	
	int l_sum = l1 + l2 + l3;
	
	if(l_sum % 2 == 1) return 0.0;
	
	if(std::abs(m1) > l1 || std::abs(m2) > l2 || std::abs(m3) > l3) return 0.0;

	if( (m1 + m2 + m3) != 0) return 0.0;

	int L1, M1, L2, M2, L3, M3;
	
	L1 = l1;
	M1 = m1;
	L2 = l2;
	M2 = m2;
	Swap(L1, M1, L2, M2);

	L3 = l3;
	M3 = m3;
	Swap(L1, M1, L3, M3);

	Swap(L2, M2, L3, M3);

	if(M3 < 0)
	{
		M1 = -M1;
		M2 = -M2;
		M3 = -M3;
	}

	/*	
	if(l1 == 2 && m1 == -1 && l2 == 2 && m2 == 2 && l3 == 2 && m3 == -1)
	{
		std::cout << L1 << " " << L2 << " " << L3 << std::endl;
		std::cout << M1 << " " << M2 << " " << M3 <<std::endl;
	}
	*/
	
	int ic1 = index_func(L1, L2, L3, M3);
	int ic2 = M2 + L2;

	try{ return Gaunt_CH.at(ic1).at(ic2); }			// Peize Lin add 2016-08-26
	catch( std::out_of_range ){ return 0; }
}
	

///Input value, 
///m1, m2, m3 are restricted within 0 to 2l+1, 
///and should be transformed first.
double ORB_gaunt_table::Get_Gaunt_SH
(
 	const int& l1,
	const int& mm1,
	const int& l2,
	const int& mm2,
	const int& l3,
	const int& mm3
)
{
//	ModuleBase::TITLE("ORB_gaunt_table","Get_Gaunt_SH");
	ModuleBase::timer::tick("ORB_gaunt_table","Get_Gaunt_SH");
	
	//Tranform M index
	int m1 = Index_M(mm1);
	int m2 = Index_M(mm2);
	int m3 = Index_M(mm3);

	if(m1 >= 0 && m2 >= 0 && m3 >= 0)
	{
		if(m1 * m2 * m3 > 0)
		{
			if(m1 == m2 + m3) return pow(-1.0, m1) * sqrt(2.0) / 2.0 
											* Get_Gaunt_CH(l1, -m1, l2, m2, l3, m3);
			else if(m2 == m1 + m3) return pow(-1.0, m2) * sqrt(2.0) / 2.0
												* Get_Gaunt_CH(l1, m1, l2, -m2, l3, m3);
			else if(m3 == m1 + m2) return pow(-1.0, m3) * sqrt(2.0) / 2.0 
												* Get_Gaunt_CH(l1, m1, l2, m2, l3, -m3);
			else return 0.0;
		}
		else
		{
			if(m1 == 0 && m2 == 0 && m3 == 0) return Get_Gaunt_CH(l1, 0, l2, 0, l3, 0);
			else if( (m1 == 0) && (m2 == m3)) return pow(-1.0, m2) * Get_Gaunt_CH(l1, 0, l2, m2, l3, -m2);
			else if( (m2 == 0) && (m3 == m1)) return pow(-1.0, m1) * Get_Gaunt_CH(l2, 0, l1, m1, l3, -m1);
			else if( (m3 == 0) && (m1 == m2)) return pow(-1.0, m2) * Get_Gaunt_CH(l3, 0, l2, m2, l1, -m2);
			else return 0.0;
		}
	}
	else
	{
		if(m1 >= 0 && m2 < 0 && m3 < 0) 
		{
			if((m1 == 0) && (m2 == m3)) return pow(-1.0, m2) * Get_Gaunt_CH(l1, 0, l2, m2, l3, -m2);
			else if(m1 > 0 && (m2 == m1+m3)) 
				return pow(-1.0, m3) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l1, m1, l2, -m2, l3, m3);
			else if(m1 > 0 && (m3 == m1+m2)) 
				return pow(-1.0, m2) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l1, m1, l2, m2, l3, -m3);
			else if(m1 > 0 && ( (m1+m2+m3) == 0)) 
				return pow(-1.0, m1+1) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l1, m1, l2, m2, l3, m3);
			else return 0.0;
		}
		else if(m2 >= 0 && m1 < 0 && m3 < 0)
		{
			if((m2 == 0) && (m1 == m3)) return pow(-1.0, m1) * Get_Gaunt_CH(l2, 0, l1, m1, l3, -m1);
			else if(m2 > 0 && (m1 == (m2 + m3))) 
				return pow(-1.0, m3) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l2, m2, l1, -m1, l3, m3);
			else if(m2 > 0 && (m3 == (m2 + m1))) 
				return pow(-1.0, m1) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l2, m2, l1, m1, l3, -m3);
			else if(m2 > 0 && ((m1+m2+m3) == 0)) 
				return pow(-1.0, m2+1) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l2, m2, l1, m1, l3, m3);
			else return 0.0;
			
		}
		else if(m3 >= 0 && m1 < 0 && m2 < 0)
		{
			if((m3 == 0) && (m1 == m2)) return pow(-1.0, m1) * Get_Gaunt_CH(l3, 0, l1, m1, l2, -m1);
			else if(m3 > 0 && (m1 == m3+m2)) 
				return pow(-1.0, m2) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l3, m3, l1, -m1, l2, m2);
			else if(m3 > 0 && (m2 == m3+m1)) 
				return pow(-1.0, m1) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l3, m3, l1, m1, l2, -m2);
			else if(m3 > 0 && ( (m1+m2+m3) == 0)) 
				return pow(-1.0, m3+1) * sqrt(2.0) / 2.0 * Get_Gaunt_CH(l3, m3, l1, m1, l2, m2);
			else return 0.0;
		}	
		else return 0.0;
	}

	ModuleBase::timer::tick("ORB_gaunt_table","Get_Gaunt_SH");
}


double ORB_gaunt_table::Fact(const int& n)
{
	double val = 1.0;
	for(int i = 1; i <= n; i++)
	{
		val *= static_cast<double>(i);
	}
	return val;
}


void ORB_gaunt_table::Swap(
	int& l1, 
	int& m1, 
	int& l2, 
	int & m2)
{
	int tmp1=0, tmp2=0;
	if(l1 >= l2) return;
	else
	{
		tmp1 = l2;
		tmp2 = m2;
		
		l2 = l1;
		m2 = m1;

		l1 = tmp1;
		m1 = tmp2;
	}
	return;
}


int ORB_gaunt_table::Index_M(const int& m)
{
	if(m % 2 == 0) return (- m / 2);
	else return ((m+1) / 2);
}

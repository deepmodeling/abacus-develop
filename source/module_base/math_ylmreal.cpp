#include "math_ylmreal.h"
#include "timer.h"
#include "constants.h"
#include "tool_quit.h"
#include "realarray.h"
#include <cassert>
#include "ylm.h"
#include "module_base/kernels/math_op.h"
#include "module_psi/kernels/memory_op.h"
#include "module_base/libm/libm.h"

namespace ModuleBase
{

YlmReal::YlmReal(){}
YlmReal::~YlmReal(){}

void YlmReal::rlylm
(
    const int lmax, 	
    const double& x,				
    const double& y,
	const double& z, // g_cartesian_vec(x,y,z)
    double* rly 	 // output
)
{
	ModuleBase::timer::tick("YlmReal","rlylm");

	assert(lmax >= 0);

	//get xy_dependence
	assert(lmax <= 19);
	
	double Am[20];
	double Bm[20];

	// mohan add 2021-05-07
	for(int i=0; i<20; ++i)
	{
		Am[i]=0.0;
		Bm[i]=0.0;
	}
	
	//ZEROS(Am, 20);
	//ZEROS(Bm, 20);
	
	double x2, x3, x4, x5;
	double y2, y3, y4, y5;
	
	x2 = x * x;
	x3 = x2 * x;
	x4 = x3 * x;
	x5 = x4 * x;

	y2 = y * y;
	y3 = y2 * y;
	y4 = y3 * y;
	y5 = y4 * y;
		
	//x-y dependence
	//Am
	//Bm
	for(int im = 0; im < lmax+1; im++)
	{
		if(im == 0)
		{
			Am[0] = 1.0; 
			Bm[0] = 0.0;
		}
		else if(im == 1)
		{
			Am[1] = x; 
			Bm[1] = y;
		}
		else if(im == 2)
		{
			Am[2] = x2- y2; 
			Bm[2] = 2.0 * x * y;
		}
		else if(im == 3)
		{
			Am[3] = x3 - 3.0 * x * y2;
			Bm[3] = 3.0 * x2 * y - y3;
		}
		else if(im == 4)
		{
			Am[4] = x4 - 6.0 * x2 * y2 + y4;
			Bm[4] = 4.0 * (x3 * y - x * y3);
		}
		else if(im == 5)
		{
			Am[5] = x5 - 10.0 * x3 * y2 + 5.0 * x * y4;
			Bm[5] = 5.0 * x4 * y - 10.0 * x2 * y3 + y5;
		}
		else
		{
			for(int ip = 0; ip <= im; ip++)
			{
				double aux = Fact(im) / Fact(ip) / Fact(im - ip);
				Am[im] += aux * pow(x, ip) * pow(y, im-ip) * cos( (im-ip) * ModuleBase::PI / 2.0 );
				Bm[im] += aux * pow(x, ip) * pow(y, im-ip) * sin( (im-ip) * ModuleBase::PI / 2.0 );
			}
		}
	}
			
	//z dependence
	double zdep[20][20];
	
	for(int il = 0; il < 20; il++)
	{
		for(int jl=0; jl < 20; jl++)
		{
			zdep[il][jl]=0.0; // mohan add 2021-05-07
		}
//		ZEROS(zdep[il], 20);
	}

	double z2 = z * z;
	double z3 = z2 * z;
	double z4 = z3 * z;
	//double z5 = z4 * z;
	
	double r = sqrt(x*x + y*y + z*z);
	double r2 = r * r;
	double r3 = r2 * r;
	double r4 = r3 * r;
	
	for(int il = 0; il < lmax + 1; il++)
	{
		if(il == 0)
		{
			zdep[0][0] = 1.0;
		}
		else if(il == 1)
		{
			zdep[1][0] = z;
			zdep[1][1] = 1.0;
		}
		else if(il == 2)
		{
			zdep[2][0] = 0.5 * (3.0 * z2 - r2);
			zdep[2][1] = sqrt(3.0) * z;
			zdep[2][2] = sqrt(3.0) * 0.5;
		}
		else if(il == 3)
		{
			zdep[3][0] = 2.5 * z3 - 1.5 * z * r2;
			zdep[3][1] = 0.25 * sqrt(6.0) * (5.0 * z2 - r2);
			zdep[3][2] = 0.5 * sqrt(15.0) * z;
			zdep[3][3] = 0.25 * sqrt(10.0);
		}
		else if(il == 4)
		{
			zdep[4][0] = 0.125 * (35.0 * z4 - 30.0 * r2 * z2 + 3.0 * r4);
			zdep[4][1] = sqrt(10.0) * 0.25 * z * (7.0 * z2 - 3.0 * r2);
			zdep[4][2] = sqrt(5.0) * 0.25 * (7.0 * z2 - r2);
			zdep[4][3] = sqrt(70.0) * 0.25 * z;
			zdep[4][4] = sqrt(35.0) * 0.125;
		}
		else if(il == 5)
		{
			zdep[5][0] = 0.125 * z *( 63.0 * z4 - 70.0 * z2 * r2 + 15.0 * r4);
			zdep[5][1] = 0.125 * sqrt(15.0) * (21.0 * z4 - 14.0 * z2 * r2 + r4);
			zdep[5][2] = 0.25 * sqrt(105.0) * z * (3.0 * z2 - r2);
			zdep[5][3] = 0.0625 * sqrt(70.0) * (9.0 * z2 - r2);
			zdep[5][4] = 0.375 * sqrt(35.0) * z;
			zdep[5][5] = 0.1875 * sqrt(14.0);
		}
		else
		{
			for(int im = 0; im <= il; im++)
			{
				int kmax = static_cast<int>( (il - im) / 2 );
				for(int ik = 0; ik <= kmax; ik++)
				{
					int twok = 2 * ik;
				
					double gamma;
					double aux0, aux1, aux2, aux3;
				
					aux0 = pow(-1.0, ik) * pow(2.0, -il);
					aux1 = Fact(il) / Fact(ik) / Fact(il-ik);
					aux2 = Fact(2*il - twok) / Fact(il) / Fact(il - twok);
					aux3 = Fact(il - twok) / Fact(il - twok - im);
				
					gamma = aux0 * aux1 * aux2 * aux3;
					
					assert(il - twok - im >= 0);
					zdep[il][im] += pow(r, twok) * pow(z, il-twok-im) * gamma;
				}

				if(im >= 1)
				{
					zdep[il][im] *= sqrt(2 * Fact(il - im) / Fact(il + im));
					
				}
			}
		}			
	}

	//calc
	int ic = 0;

	//special case for r=0
	double rpi = r;
	const double tiny =  1.0E-10;
	if (rpi < tiny) rpi += tiny;
	
	for(int il = 0; il <= lmax; il++)
	{
		double fac = sqrt( (2.0 * il + 1.0) / ModuleBase::FOUR_PI );

		double rl = pow(rpi, il);
			
		//m=0
		rly[ic] = Am[0] * zdep[il][0] * fac / rl;
		
		ic++;
		
		//m ! = 0
		for(int im = 1; im <= il; im++)
		{
			//m>0
			rly[ic] = Am[im] * zdep[il][im] * pow(-1.0, im) * fac / rl;
			
			ic++;
			
			//m<0
			rly[ic] = Bm[im] * zdep[il][im] * pow(-1.0, im) * fac / rl;

			ic++;
		}
	}

	ModuleBase::timer::tick("YlmReal","rlylm");
	return;
}


void YlmReal::Ylm_Real2
(
    const int lmax2, 			// lmax2 = (lmax+1)^2
    const int ng,				//
    const ModuleBase::Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    matrix &ylm 				// output
)
{
    if (ng<1 || lmax2<1)
    {
        ModuleBase::WARNING("YLM_REAL","ng<1 or lmax2<1");
        return;
    }

//----------------------------------------------------------
// EXPLAIN : find out lmax
//----------------------------------------------------------
    bool out_of_range = true;
    int lmax = 0;
    for (int l= 0; l< 30; l++)
    {
        if ((l+1)*(l+1) == lmax2)
        {
            lmax = l;
            out_of_range = false;
            break;
        }
    }
    if (out_of_range)
    {
        ModuleBase::WARNING_QUIT("YLM_REAL","l>30 or l<0");
    }

//----------------------------------------------------------
//	Start CALC
//----------------------------------------------------------
	double* rly = new double[lmax2];
	
	for (int ig = 0; ig < ng; ig++)
	{
		rlylm (lmax, g[ig].x, g[ig].y, g[ig].z, rly);
		
		for (int lm = 0; lm < lmax2; lm++)
		{
			ylm (lm, ig) = rly[lm];
		}
	}

	delete [] rly;

	return;
}

//==========================================================
// MEMBER FUNCTION :
// NAME : YLM_REAL(Real spherical harmonics ylm(G) up to l=lmax
// Use Numerical recursive algorithm as given in Numerical Recipes
//==========================================================
// from ylmr2.f90

template <typename FPTYPE, typename Device>
void YlmReal::Ylm_Real(Device * ctx, const int lmax2, const int ng, const FPTYPE *g, FPTYPE * ylm)
{
    using resmem_var_op = psi::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = psi::memory::delete_memory_op<FPTYPE, Device>;
    using cal_ylm_real_op = ModuleBase::cal_ylm_real_op<FPTYPE, Device>;

    if (ng < 1 || lmax2 < 1) {
        ModuleBase::WARNING("YLM_REAL","ng<1 or lmax2<1");
        return;
    }

//----------------------------------------------------------
// EXPLAIN : find out lmax
//----------------------------------------------------------
    bool out_of_range = true;
    int lmax = 0;
    for (int l = 0; l < 30; l++) {
        if ((l + 1) * (l + 1) == lmax2) {
            lmax = l;
            out_of_range = false;
            break;
        }
    }
    if (out_of_range) {
        ModuleBase::WARNING_QUIT("YLM_REAL","l>30 or l<0");
    }
    FPTYPE * p = nullptr, * phi = nullptr, * cost = nullptr;
    resmem_var_op()(ctx, p, (lmax + 1) * (lmax + 1) * ng, "YlmReal::Ylm_Real");

    cal_ylm_real_op()(
        ctx,
        ng,
        lmax,
        ModuleBase::SQRT2,
        ModuleBase::PI,
        ModuleBase::PI_HALF,
        ModuleBase::FOUR_PI,
        ModuleBase::SQRT_INVERSE_FOUR_PI,
        g,
        p,
        ylm);

    delmem_var_op()(ctx, p);
    delmem_var_op()(ctx, phi);
    delmem_var_op()(ctx, cost);
} // end subroutine ylmr2

//==========================================================
// MEMBER FUNCTION :
// NAME : YLM_REAL(Real spherical harmonics ylm(G) up to l=lmax
// Use Numerical recursive algorithm as given in Numerical Recipes
//==========================================================
// from ylmr2.f90
void YlmReal::Ylm_Real
(
    const int lmax2, 			// lmax2 = (lmax+1)^2
    const int ng,				//
    const ModuleBase::Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    matrix &ylm 				// output
)
{

    if (ng<1 || lmax2<1)
    {
        ModuleBase::WARNING("YLM_REAL","ng<1 or lmax2<1");
        return;
    }

//----------------------------------------------------------
// EXPLAIN : find out lmax
//----------------------------------------------------------
    bool out_of_range = true;
    int lmax = 0;
    for (int l= 0; l< 30; l++)
    {
        if ((l+1)*(l+1) == lmax2)
        {
            lmax = l;
            out_of_range = false;
            break;
        }
    }
    if (out_of_range)
    {
        ModuleBase::WARNING_QUIT("YLM_REAL","l>30 or l<0");
    }

//----------------------------------------------------------
// EXPLAIN : if lmax = 1,only use Y00 , output result.
//----------------------------------------------------------
    if (lmax == 0)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0;i<ng;i++)
        {
            ylm(0, i) = ModuleBase::SQRT_INVERSE_FOUR_PI;
        }
        return;
    }

//----------------------------------------------------------
// LOCAL VARIABLES :
// NAME : cost = cos(theta),theta and phi are polar angles
// NAME : phi
//----------------------------------------------------------
    double *cost = new double[ng];
    double *phi = new double[ng];

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ig = 0;ig < ng;ig++)
    {
        const double gmod = g[ig].norm();
        if (gmod < 1.0e-9)
        {
            cost[ig] = 0.0;
        }
        else
        {
            cost[ig] = g[ig].z / gmod;
        }// endif

        //  beware the arc tan, it is defined modulo pi
        if (g[ig].x > 1.0e-9)
        {
            phi[ig] = atan(g[ig].y / g[ig].x);
        }
        else if (g[ig].x < -1.e-9)
        {
            phi[ig] = atan(g[ig].y / g[ig].x) + ModuleBase::PI;
        }
        else
        {
            phi[ig] = ModuleBase::PI_HALF * ((g[ig].y >= 0.0) ? 1.0 : -1.0); //HLX: modified on 10/13/2006
        } // end if
    } // enddo

//==========================================================
// NAME : p(Legendre Polynomials) (0 <= m <= l)
//==========================================================
    ModuleBase::realArray p(lmax+1,lmax+1,ng);
    int lm = -1;
    for (int l=0; l<=lmax; l++)
    {
        const double c = sqrt((2*l+1) / ModuleBase::FOUR_PI);
        if (l == 0)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i=0;i<ng;i++)
            {
                p(0,0,i) = 1.0;
            }
        }
        else if (l == 1)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i=0;i<ng;i++)
            {
                p(0,1,i) = cost[i];
                auto x1 = 1.0 - cost[i] * cost[i];
                x1 = std::max(0.0, x1);
                p(1,1,i) = -sqrt(x1);
            }
        }
        else
        {
            const int l1 = l-1;
            const int l2 = l-2;
            const int l3 = 2*l-1;
            //  recursion on l for P(:,l,m)
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int m=0; m<=l2; m++)  // do m = 0, l - 2//mohan modify 2007-10-13
            {
                for (int i=0; i<ng; i++)
                {
                    p(m, l, i) = (cost[i] * l3 * p(m, l1, i) -
                                  (l1 + m ) * p(m, l2, i)) / (l - m);
                }
            } // end do
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i=0;i<ng;i++)
            {
                p(l1, l, i) = cost[i] * l3 * p(l1, l1, i);
                auto x2 = 1.0 - cost[i] * cost[i];
                x2 = std::max(0.0, x2);
                p(l, l, i) = Semi_Fact(l3) * pow(x2, static_cast<double>(l) / 2.0) ;//mohan modify 2007-10-13
                if (l%2 == 1)
                {
                    p(l, l, i) = -p(l, l, i);
                }
            }
        } // end if

        // Y_lm, m = 0
        ++lm;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0;i<ng;i++)
        {
            ylm(lm, i) = c*p(0, l, i);
        }

        for (int m=1;m<=l;m++)
        {
            // Y_lm, m > 0
            const double same = c * sqrt
                                (
                                    static_cast<double>(Fact(l - m)) /
                                    static_cast<double>(Fact(l + m))
                                )
                                *ModuleBase::SQRT2;

            ++lm;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i=0;i<ng;i++)
            {
				double sinp, cosp;
                ModuleBase::libm::sincos(m * phi[i], &sinp, &cosp);
                ylm(lm, i) = same * p(m,l,i) * cosp;
				ylm(lm + 1, i) = same * p(m,l,i) * sinp;
            }

            // Y_lm, m < 0
            ++lm;

            /*
             * mohan test bug 2009-03-03
             *
            if(l==9 && m==8)
            {
            	if(my_rank==0)
            	{
            		std::ofstream ofs("Log2.txt");
            		for(int ig=0; ig<ng; ig++)
            		{
            			if(ig%1==0) ofs << "\n";
            			ofs << std::setw(20) << same
            				<< std::setw(20) << Fact(l - m)
            				<< std::setw(20) << Fact(l + m)
            				<< std::setw(20) << ylm(lm, ig);
            		}
            	}
            	MPI_Barrier(MPI_COMM_WORLD);
            	ModuleBase::QUIT();
            }
            */

        }
    }// end do



    /*	GlobalV::ofs_running<<"\n Unit Condition About Ylm_Real"<<std::endl;
    	int count=0;
    	for(int l=0; l<=lmax; l++)
    	{
    		for(int m=0; m<2*l+1; m++)
    		{
    			//  mohan debug 2009-03-03
    			if(l==9 && m==15)
    			{
    				if(my_rank==0)
    				{
    					std::ofstream ofs("Log1.txt");
    					for(int ig=0; ig<ng; ig++)
    					{
    						if(ig%6==0) ofs << "\n";
    						ofs << std::setw(20) << ylm(count, ig);
    					}
    				}
    				MPI_Barrier(MPI_COMM_WORLD);
    				ModuleBase::QUIT();
    			}
    			double sum_before = 0.0;
    			for(int ig=0; ig<ng; ig++)
    			{
    				sum_before += ylm(count, ig) * ylm(count, ig);
    			}
    			sum_before *= ModuleBase::FOUR_PI/ng;
    			GlobalV::ofs_running<<std::setw(5)<<l<<std::setw(5)<<m<<std::setw(15)<<sum_before;


    //			for(int ig=0; ig<ng; ig++)
    //			{
    //				ylm(count, ig) /= sqrt(sum_before);
    //			}
    //			double sum = 0;
    //			for(int ig=0; ig<ng; ig++)
    //			{
    //				sum += ylm(count, ig) * ylm(count, ig);
    //			}
    //			count++;
    //			GlobalV::ofs_running<<std::setw(15)<<sum*ModuleBase::FOUR_PI/ng;

    			GlobalV::ofs_running<<std::endl;
    		}
    	}
    	GlobalV::ofs_running<<std::endl;
    */


    delete [] cost;
    delete [] phi;

    return;
} // end subroutine ylmr2

void YlmReal::grad_Ylm_Real
    (
        const int lmax2, 		
        const int ng,				
        const ModuleBase::Vector3<double> *g, 
		matrix &ylm,
        matrix &dylmx,
		matrix &dylmy,
		matrix &dylmz	
    )
{
	ModuleBase::Ylm::set_coefficients();
	int lmax = int(sqrt( double(lmax2) ) + 0.1) - 1;

	for (int ig = 0;ig < ng;ig++)
    {
		ModuleBase::Vector3<double> gg = g[ig];
        double gmod = gg.norm();
        if (gmod < 1.0e-9)
        {
			for(int lm = 0 ; lm < lmax2 ; ++lm)
			{
				if(lm == 0) 
					ylm(lm,ig) = ModuleBase::SQRT_INVERSE_FOUR_PI;
				else	
					ylm(lm,ig) = 0;
				dylmx(lm,ig) = dylmy(lm,ig) = dylmz(lm,ig) = 0;
			}
		}
		else
		{
			std::vector<double> tmpylm;
			std::vector<std::vector<double>> tmpgylm;
			Ylm::grad_rl_sph_harm(lmax2, gg.x, gg.y, gg.z, tmpylm,tmpgylm);
			int lm = 0;
			for(int il = 0 ; il <= lmax ; ++il)
			{
				for(int im = 0; im < 2*il+1; ++im, ++lm)
				{
					double rlylm = tmpylm[lm];
					ylm(lm,ig) = rlylm / pow(gmod,il);
					dylmx(lm,ig) = ( tmpgylm[lm][0] - il*rlylm * gg.x / pow(gmod,2) )/pow(gmod,il);
					dylmy(lm,ig) = ( tmpgylm[lm][1] - il*rlylm * gg.y / pow(gmod,2) )/pow(gmod,il);
					dylmz(lm,ig) = ( tmpgylm[lm][2] - il*rlylm * gg.z / pow(gmod,2) )/pow(gmod,il);
				}
			}
			
		}
	}
	return;
}


//==========================================================
// MEMBER FUNCTION :
// NAME : Fact ( n! )
// NAME : Semi_Fact ( n!! )
//==========================================================
long double YlmReal::Fact(const int n)
{
    long double f = 1;
    for (int i=n; i>1; i--)
    {
        f *= i;
    }
    return f;
}

int YlmReal::Semi_Fact(const int n)
{
    int semif = 1;
    for (int i=n; i>2; i -= 2)
    {
        semif *= i;
    }
    return semif;
}

template void YlmReal::Ylm_Real<float, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int, int, const float *, float*);
template void YlmReal::Ylm_Real<double, psi::DEVICE_CPU>(psi::DEVICE_CPU*, int, int, const double *, double*);
#if ((defined __CUDA) || (defined __ROCM))
template void YlmReal::Ylm_Real<float, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int, int, const float *, float*);
template void YlmReal::Ylm_Real<double, psi::DEVICE_GPU>(psi::DEVICE_GPU*, int, int, const double *, double*);
#endif
}  // namespace ModuleBase

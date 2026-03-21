#include "math_ylmreal.h"

#include "constants.h"
#include "source_base/kernels/math_ylm_op.h"
#include "source_base/libm/libm.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/array_pool.h"
#include "realarray.h"
#include "timer.h"
#include "tool_quit.h"
#include "ylm.h"

#include <cassert>
#include <vector>
#include <cmath>

namespace ModuleBase
{

// Helper function: calculate Am, Bm.
static void compute_xy_dependence(int lmax, double x, double y,
                                  double* Am, double* Bm)
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;
    double y2 = y * y;
    double y3 = y2 * y;
    double y4 = y3 * y;
    double y5 = y4 * y;

    for (int im = 0; im <= lmax; ++im) {
        if (im == 0) {
            Am[im] = 1.0;
            Bm[im] = 0.0;
        }
        else if (im == 1) {
            Am[im] = x;
            Bm[im] = y;
        }
        else if (im == 2) {
            Am[im] = x2 - y2;
            Bm[im] = 2.0 * x * y;
        }
        else if (im == 3) {
            Am[im] = x3 - 3.0 * x * y2;
            Bm[im] = 3.0 * x2 * y - y3;
        }
        else if (im == 4) {
            Am[im] = x4 - 6.0 * x2 * y2 + y4;
            Bm[im] = 4.0 * (x3 * y - x * y3);
        }
        else if (im == 5) {
            Am[im] = x5 - 10.0 * x3 * y2 + 5.0 * x * y4;
            Bm[im] = 5.0 * x4 * y - 10.0 * x2 * y3 + y5;
        }
        else {
            for (int ip = 0; ip <= im; ++ip) {
                double aux = Fact(im) / Fact(ip) / Fact(im - ip);
                double term = aux * pow(x, ip) * pow(y, im - ip);
                Am[im] += term * cos((im - ip) * ModuleBase::PI / 2.0);
                Bm[im] += term * sin((im - ip) * ModuleBase::PI / 2.0);
            }
        }
    }
}

// Helper function: calculating zdep
static void compute_z_dependence(int lmax, double x, double y, double z,
                                 double r, double zdep[][20])
{
    double r2 = r * r;
    double r3 = r2 * r;
    double r4 = r3 * r;
    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;

    for (int il = 0; il <= lmax; ++il) {
        if (il == 0) {
            zdep[il][0] = 1.0;
        }
        else if (il == 1) {
            zdep[il][0] = z;
            zdep[il][1] = 1.0;
        }
        else if (il == 2) {
            zdep[il][0] = 0.5 * (3.0 * z2 - r2);
            zdep[il][1] = sqrt(3.0) * z;
            zdep[il][2] = sqrt(3.0) * 0.5;
        }
        else if (il == 3) {
            zdep[il][0] = 2.5 * z3 - 1.5 * z * r2;
            zdep[il][1] = 0.25 * sqrt(6.0) * (5.0 * z2 - r2);
            zdep[il][2] = 0.5 * sqrt(15.0) * z;
            zdep[il][3] = 0.25 * sqrt(10.0);
        }
        else if (il == 4) {
            zdep[il][0] = 0.125 * (35.0 * z4 - 30.0 * r2 * z2 + 3.0 * r4);
            zdep[il][1] = sqrt(10.0) * 0.25 * z * (7.0 * z2 - 3.0 * r2);
            zdep[il][2] = sqrt(5.0) * 0.25 * (7.0 * z2 - r2);
            zdep[il][3] = sqrt(70.0) * 0.25 * z;
            zdep[il][4] = sqrt(35.0) * 0.125;
        }
        else if (il == 5) {
            zdep[il][0] = 0.125 * z * (63.0 * z4 - 70.0 * z2 * r2 + 15.0 * r4);
            zdep[il][1] = 0.125 * sqrt(15.0) * (21.0 * z4 - 14.0 * z2 * r2 + r4);
            zdep[il][2] = 0.25 * sqrt(105.0) * z * (3.0 * z2 - r2);
            zdep[il][3] = 0.0625 * sqrt(70.0) * (9.0 * z2 - r2);
            zdep[il][4] = 0.375 * sqrt(35.0) * z;
        }
        else {
            for (int im = 0; im <= il; ++im) {
                int kmax = static_cast<int>((il - im) / 2);
                for (int ik = 0; ik <= kmax; ++ik) {
                    int twok = 2 * ik;
                    double gamma = 0.0;
                    double aux0 = pow(-1.0, ik) * pow(2.0, -il);
                    double aux1 = Fact(il) / Fact(ik) / Fact(il - ik);
                    double aux2 = Fact(2 * il - twok) / Fact(il) / Fact(il - twok);
                    double aux3 = Fact(il - twok) / Fact(il - twok - im);
                    gamma = aux0 * aux1 * aux2 * aux3;
                    zdep[il][im] += pow(r, twok) * pow(z, il - twok - im) * gamma;
                }
                if (im >= 1) {
                    zdep[il][im] *= sqrt(2.0 * Fact(il - im) / Fact(il + im));
                }
            }
        }
    }
}

YlmReal::YlmReal(){}
YlmReal::~YlmReal(){}

void YlmReal::rlylm(const int lmax,
                    const double& x,
                    const double& y,
                    const double& z,
                    double* rly)
{
    ModuleBase::timer::tick("YlmReal", "rlylm");
    assert(lmax >= 0);
    assert(lmax <= 19);

    constexpr int MAX_L = 20;           // Maximum l + 1 (since lmax <= 19)
    constexpr double TINY = 1.0E-10;    // Small value to avoid division by zero

    // Arrays for x-y angular part
    double Am[MAX_L] = {0.0};
    double Bm[MAX_L] = {0.0};
    // Array for z-dependent part (zdep[l][m])
    double zdep[MAX_L][MAX_L] = {{0.0}};

    // Compute x-y dependence (Am, Bm)
    compute_xy_dependence(lmax, x, y, Am, Bm);

    // Compute radial distance and z dependence
    double r = std::sqrt(x * x + y * y + z * z);
    compute_z_dependence(lmax, x, y, z, r, zdep);

    // Combine results into output array rly
    int ic = 0;
    double rpi = r;
    if (rpi < TINY) rpi += TINY;

    for (int il = 0; il <= lmax; ++il) {
        double fac = std::sqrt((2.0 * il + 1.0) / ModuleBase::FOUR_PI);
        double rl = std::pow(rpi, il);

        // m = 0
        rly[ic] = Am[0] * zdep[il][0] * fac / rl;
        ic++;

        for (int im = 1; im <= il; ++im) {
            // m > 0
            rly[ic] = Am[im] * zdep[il][im] * std::pow(-1.0, im) * fac / rl;
            ic++;
            // m < 0
            rly[ic] = Bm[im] * zdep[il][im] * std::pow(-1.0, im) * fac / rl;
            ic++;
        }
    }

    ModuleBase::timer::tick("YlmReal", "rlylm");
}
	avoid YlmReal::rlylm(const int lmax,
                    const double& x,
                    const double& y,
                    const double& z,
                    double* rly)
{
    ModuleBase::timer::tick("YlmReal", "rlylm");
    assert(lmax >= 0);
    assert(lmax <= 19);

    constexpr int MAX_L = 20;           // Maximum l + 1 (since lmax <= 19)
    constexpr double TINY = 1.0E-10;    // Small value to avoid division by zero

    // Arrays for x-y angular part
    double Am[MAX_L] = {0.0};
    double Bm[MAX_L] = {0.0};
    // Array for z-dependent part (zdep[l][m])
    double zdep[MAX_L][MAX_L] = {{0.0}};

    // Compute x-y dependence (Am, Bm)
    compute_xy_dependence(lmax, x, y, Am, Bm);

    // Compute radial distance and z dependence
    double r = std::sqrt(x * x + y * y + z * z);
    compute_z_dependence(lmax, x, y, z, r, zdep);

    // Combine results into output array rly
    int ic = 0;
    double rpi = r;
    if (rpi < TINY) rpi += TINY;

    for (int il = 0; il <= lmax; ++il) {
        double fac = std::sqrt((2.0 * il + 1.0) / ModuleBase::FOUR_PI);
        double rl = std::pow(rpi, il);

        // m = 0
        rly[ic] = Am[0] * zdep[il][0] * fac / rl;
        ic++;

        for (int im = 1; im <= il; ++im) {
            // m > 0
            rly[ic] = Am[im] * zdep[il][im] * std::pow(-1.0, im) * fac / rl;
            ic++;
            // m < 0
            rly[ic] = Bm[im] * zdep[il][im] * std::pow(-1.0, im) * fac / rl;
            ic++;
        }
    }

    ModuleBase::timer::tick("YlmReal", "rlylm");
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
	std::vector<double> rly(lmax2);
	
	for (int ig = 0; ig < ng; ig++)
	{
		rlylm (lmax, g[ig].x, g[ig].y, g[ig].z, rly.data());
		
		for (int lm = 0; lm < lmax2; lm++)
		{
			ylm (lm, ig) = rly[lm];
		}
	}

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
    using resmem_var_op = base_device::memory::resize_memory_op<FPTYPE, Device>;
    using delmem_var_op = base_device::memory::delete_memory_op<FPTYPE, Device>;
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
    resmem_var_op()(p, (lmax + 1) * (lmax + 1) * ng, "YlmReal::Ylm_Real");

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

    delmem_var_op()(p);
    delmem_var_op()(phi);
    delmem_var_op()(cost);
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
	std::vector <double> cost(ng);
	std::vector <double> phi(ng);

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
	const int lmax = int(sqrt( double(lmax2) ) + 0.1) - 1;
	std::vector<double> tmpylm((lmax2+1) * (lmax2+1));
	Array_Pool<double> tmpgylm((lmax2+1) * (lmax2+1), 3);

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
			Ylm::grad_rl_sph_harm(lmax2, gg.x, gg.y, gg.z, tmpylm.data(),tmpgylm.get_ptr_2D());
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

template void YlmReal::Ylm_Real<float, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                int,
                                                                int,
                                                                const float*,
                                                                float*);
template void YlmReal::Ylm_Real<double, base_device::DEVICE_CPU>(base_device::DEVICE_CPU*,
                                                                 int,
                                                                 int,
                                                                 const double*,
                                                                 double*);
#if ((defined __CUDA) || (defined __ROCM))
template void YlmReal::Ylm_Real<float, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                int,
                                                                int,
                                                                const float*,
                                                                float*);
template void YlmReal::Ylm_Real<double, base_device::DEVICE_GPU>(base_device::DEVICE_GPU*,
                                                                 int,
                                                                 int,
                                                                 const double*,
                                                                 double*);
#endif
}  // namespace ModuleBase

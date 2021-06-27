#include "math_ylmreal.h"
#include "timer.h"
#include "../module_base/constants.h"
#include "tool_quit.h"
#include "realarray.h"
#include <cassert>

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
	timer::tick("YlmReal","rlylm");

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
				Am[im] += aux * pow(x, ip) * pow(y, im-ip) * cos( (im-ip) * PI / 2.0 );
				Bm[im] += aux * pow(x, ip) * pow(y, im-ip) * sin( (im-ip) * PI / 2.0 );
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
		double fac = sqrt( (2.0 * il + 1.0) / FOUR_PI );

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

	timer::tick("YlmReal","rlylm");
	return;
}


void YlmReal::Ylm_Real2
(
    const int lmax2, 			// lmax2 = (lmax+1)^2
    const int ng,				//
    const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    matrix &ylm 				// output
)
{
    if (ng<1 || lmax2<1)
    {
        WARNING("YLM_REAL","ng<1 or lmax2<1");
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
        WARNING_QUIT("YLM_REAL","l>30 or l<0");
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
void YlmReal::Ylm_Real
(
    const int lmax2, 			// lmax2 = (lmax+1)^2
    const int ng,				//
    const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    matrix &ylm 				// output
)
{

    if (ng<1 || lmax2<1)
    {
        WARNING("YLM_REAL","ng<1 or lmax2<1");
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
        WARNING_QUIT("YLM_REAL","l>30 or l<0");
    }

//----------------------------------------------------------
// EXPLAIN : if lmax = 1,only use Y00 , output result.
//----------------------------------------------------------
    if (lmax == 0)
    {
        for (int i=0;i<ng;i++)
        {
            ylm(0, i) = SQRT_INVERSE_FOUR_PI;
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
            phi[ig] = atan(g[ig].y / g[ig].x) + PI;
        }
        else
        {
            phi[ig] = PI_HALF * ((g[ig].y >= 0.0) ? 1.0 : -1.0); //HLX: modified on 10/13/2006
        } // end if
    } // enddo

//==========================================================
// NAME : p(Legendre Polynomials) (0 <= m <= l)
//==========================================================
    realArray p(lmax+1,lmax+1,ng);
    int m;
    int i;
    double x1, x2;
    int lm = -1;
    for (int l=0; l<=lmax; l++)
    {
        const double c = sqrt((2*l+1) / FOUR_PI);
        if (l == 0)
        {
            for (i=0;i<ng;i++)
            {
                p(0,0,i) = 1.0;
            }
        }
        else if (l == 1)
        {
            for (i=0;i<ng;i++)
            {
                p(0,1,i) = cost[i];
                x1 = 1.0 - cost[i] * cost[i];
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
            for (m=0; m<=l2; m++)  // do m = 0, l - 2//mohan modify 2007-10-13
            {
                for (i=0; i<ng; i++)
                {
                    p(m, l, i) = (cost[i] * l3 * p(m, l1, i) -
                                  (l1 + m ) * p(m, l2, i)) / (l - m);
                }
            } // end do
            for (i=0;i<ng;i++)
            {
                p(l1, l, i) = cost[i] * l3 * p(l1, l1, i);
                x2 = 1.0 - cost[i] * cost[i];
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
        for (i=0;i<ng;i++)
        {
            ylm(lm, i) = c*p(0, l, i);
        }

        for (m=1;m<=l;m++)
        {
            // Y_lm, m > 0
            const double same = c * sqrt
                                (
                                    static_cast<double>(Fact(l - m)) /
                                    static_cast<double>(Fact(l + m))
                                )
                                *SQRT2;

            ++lm;
            for (i=0;i<ng;i++)
            {
                ylm(lm, i) = same * p(m,l,i) * cos(m * phi[i]);
            }

            // Y_lm, m < 0
            ++lm;
            for (i=0;i<ng;i++)
            {
                ylm(lm, i) = same * p(m,l,i) * sin(m * phi[i]);
            }


            /*
             * mohan test bug 2009-03-03
             *
            if(l==9 && m==8)
            {
            	if(my_rank==0)
            	{
            		ofstream ofs("Log2.txt");
            		for(int ig=0; ig<ng; ig++)
            		{
            			if(ig%1==0) ofs << "\n";
            			ofs << setw(20) << same
            				<< setw(20) << Fact(l - m)
            				<< setw(20) << Fact(l + m)
            				<< setw(20) << ylm(lm, ig);
            		}
            	}
            	MPI_Barrier(MPI_COMM_WORLD);
            	QUIT();
            }
            */

        }
    }// end do



    /*	ofs_running<<"\n Unit Condition About Ylm_Real"<<endl;
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
    					ofstream ofs("Log1.txt");
    					for(int ig=0; ig<ng; ig++)
    					{
    						if(ig%6==0) ofs << "\n";
    						ofs << setw(20) << ylm(count, ig);
    					}
    				}
    				MPI_Barrier(MPI_COMM_WORLD);
    				QUIT();
    			}
    			double sum_before = 0.0;
    			for(int ig=0; ig<ng; ig++)
    			{
    				sum_before += ylm(count, ig) * ylm(count, ig);
    			}
    			sum_before *= FOUR_PI/ng;
    			ofs_running<<setw(5)<<l<<setw(5)<<m<<setw(15)<<sum_before;


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
    //			ofs_running<<setw(15)<<sum*FOUR_PI/ng;

    			ofs_running<<endl;
    		}
    	}
    	ofs_running<<endl;
    */


    delete [] cost;
    delete [] phi;

    return;
} // end subroutine ylmr2


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

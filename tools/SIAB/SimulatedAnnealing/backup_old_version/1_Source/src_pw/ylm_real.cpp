#include "ylm_real.h"
#include "../src_parallel/parallel_reduce.h"
#include "../src_spillage/tools.h"

//==========================================================
// MEMBER FUNCTION :
// NAME : YLM_REAL(Real spherical harmonics ylm(G) up to l=lmax
// Use Numerical recursive algorithm as given in Numerical Recipes
//==========================================================
void Ylm_Real
(
    const int lmax2, 			// lmax2 = (lmax+1)^2
    const int ng,				//
    const Vector3<double> *g, 	// g_cartesian_vec(x,y,z)
    matrix &ylm 				// output
)
{
//	TITLE("Mathzone","Ylm_Real");
//	timer::tick("Mathzone","Ylm_Real");

    if (ng<1 || lmax2<1)
    {
		cout << "\n YLM_REAL: ng < 1 or lmax2 < 1";
		QUIT();
        timer::tick("Mathzone","Ylm_Real");
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
        //	timer::tick("Mathzone","Ylm_Real");
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

//	timer::tick("Mathzone","Ylm_Real");
    return;
} // end subroutine ylmr2

//==========================================================
// MEMBER FUNCTION :
// NAME : Fact ( n! )
// NAME : Semi_Fact ( n!! )
//==========================================================
long double Fact(const int n)
{
    long double f = 1;
    for (int i=n; i>1; i--)
    {
        f *= i;
//		if(n>16)
//		{
//			cout<<"\n n="<<n<<" "<<f;
//		}
    }
//	if(n>16) QUIT();

    return f;
}

int Semi_Fact(const int n)
{
    int semif = 1;
    for (int i=n; i>2; i -= 2)
    {
        semif *= i;
    }
    return semif;
}


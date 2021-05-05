#include "./stress_func.h"
#include "./H_Ewald_pw.h"

//calcualte the Ewald stress term in PW and LCAO
void Stress_Func::stress_ewa(matrix& sigma, const bool is_pw)
{
    timer::tick("Stress_Func","stress_ew",'F');

    int i,j,l,m;
    double g[3];

    double charge=0;
    for(int it=0; it < ucell.ntype; it++)
	{
		for(int i=0; i<ucell.atoms[it].na; i++)
		{
			charge = charge + ucell.atoms[it].zv;
		}
	}
    //choose alpha in order to have convergence in the sum over G
    //upperbound is a safe upper bound for the error ON THE ENERGY

    double alpha=2.9;
    double upperbound;
    do{
       alpha-=0.1;
       if(alpha==0.0)
          WARNING_QUIT("stres_ew", "optimal alpha not found");
       upperbound =e2 * pow(charge,2) * sqrt( 2 * alpha / (TWO_PI)) * erfc(sqrt(ucell.tpiba2 * pw.ggchg / 4.0 / alpha));
    }
    while(upperbound>1e-7);

    //G-space sum here
    //Determine if this processor contains G=0 and set the constant term 
    double sdewald;
    if(pw.gstart == 1)
	{
       sdewald = (TWO_PI) * e2 / 4.0 / alpha * pow(charge/ucell.omega,2);
    }
    else 
	{
       sdewald = 0.0;
    }

    //sdewald is the diagonal term 

    double fact=1.0;
    if (INPUT.gamma_only && is_pw) fact=2.0;
//    else fact=1.0;

    double g2,g2a;
    double arg;
    complex<double> rhostar;
    double sewald;
    for(int ng=pw.gstart;ng<pw.ngmc;ng++)
	{
		g2 = pw.gg[ng]* ucell.tpiba2;
		g2a = g2 /4.0/alpha;
		rhostar=(0.0,0.0);
		for(int it=0; it < ucell.ntype; it++)
		{
			for(int i=0; i<ucell.atoms[it].na; i++)
			{
				arg = (pw.gcar[ng].x * ucell.atoms[it].tau[i].x + pw.gcar[ng].y * ucell.atoms[it].tau[i].y + pw.gcar[ng].z * ucell.atoms[it].tau[i].z) * (TWO_PI);
				rhostar = rhostar + complex<double>(ucell.atoms[it].zv * cos(arg),ucell.atoms[it].zv * sin(arg));
			}
		}
		rhostar /= ucell.omega;
		sewald = fact* (TWO_PI) * e2 * exp(-g2a) / g2 * pow(abs(rhostar),2);
		sdewald = sdewald - sewald;
		g[0]=pw.gcar[ng].x;
		g[1]=pw.gcar[ng].y;
		g[2]=pw.gcar[ng].z;
		for(l=0;l<3;l++)
		{
			for(m=0;m<l+1;m++)
			{
				sigma(l,m) +=sewald * ucell.tpiba2 * 2.0 * g[l] * g[m] / g2 * (g2a+1);
			}
		}
	}

	for(l=0;l<3;l++)
	{
		sigma(l,l) +=sdewald;
	}
    //R-space sum here (only for the processor that contains G=0) 
    int mxr = 50;
    int *irr;
    Vector3<double> *r;
    double *r2;
    r  = new Vector3<double>[mxr];
    r2 = new double[mxr];
    irr = new int[mxr];
    double rr;
    Vector3<double> d_tau;
    double r0[3];
    double rmax=0.0;
    int nrm=0;
    double fac;
	if(pw.gstart==1)
	{
		rmax = 4.0/sqrt(alpha)/ucell.lat0;
		//with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2*10^-1)
		for(int it=0; it < ucell.ntype; it++)
		{
			for(int i=0; i<ucell.atoms[it].na; i++)
			{
				for(int jt=0; jt < ucell.ntype; jt++)
				{
					for(int j=0; j<ucell.atoms[jt].na; j++)
					{
						//calculate tau[na]-tau[nb]
						d_tau = ucell.atoms[it].tau[i] - ucell.atoms[jt].tau[j];
						//generates nearest-neighbors shells 
						H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);
						for(int nr=0 ; nr<nrm ; nr++)
						{
							rr=sqrt(r2[nr]) * ucell.lat0;
							fac = -e2/2.0/ucell.omega*pow(ucell.lat0,2)*ucell.atoms[it].zv * ucell.atoms[jt].zv / pow(rr,3) * (erfc(sqrt(alpha)*rr)+rr * sqrt(8 * alpha / (TWO_PI)) * exp(-alpha * pow(rr,2)));
							for(l=0; l<3; l++)
							{
								for(m=0; m<l+1; m++)
								{
									r0[0] = r[nr].x;
									r0[1] = r[nr].y;
									r0[2] = r[nr].z;
									sigma(l,m) += fac * r0[l] * r0[m];
								}//end m
							}//end l
						}//end nr
					}//end j
				}//end jt
			}//end i
		}//end it
	}//end if
	for(l=0;l<3;l++)
	{
		for(m=0;m<l+1;m++)
		{
			sigma(m,l)=sigma(l,m);
		}
	}
	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigma(l,m)=-sigma(l,m);
			Parallel_Reduce::reduce_double_pool( sigma(l,m) );
		}
	}

	delete[] r;
	delete[] r2;
	delete[] irr;
	// this->print(ofs_running, "ewald stress", stression);
	timer::tick("Force_Func","stress_ew");

	return;
}

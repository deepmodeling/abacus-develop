#include "mathzone.h"


double Mathzone::PI = 3.14159265;

Mathzone::Mathzone()
{
}

Mathzone::~Mathzone()
{
}

double Mathzone::BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{
	const int XMIN = 2;
	const double FPMIN = 1.0e-30;
	const double EPS = 1.0e-10;
	const int MAXIT = 10000;

	int i, isign, l, nl;
    double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
    fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
    rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
    temp, w, x2, xi, xi2;

    if (x <= 0.0 || xnu < 0.0)
	{
//		WARNING_QUIT("Mathzone::BESSJY","bad arguments");
		exit(0);
	}


    nl = (x < XMIN ? (int)(xnu + 0.5) : IMAX(0, (int)(xnu - x + 1.5)));
    const double xmu = xnu - nl;
    const double xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    w = xi2 / PI;
    isign = 1;
    h = xnu * xi;
    if (h < FPMIN) h = FPMIN;

    b = xi2 * xnu;

    d = 0.0;

    c = h;
	
	for (i = 1;i <= MAXIT;i++)
    {
        b += xi2;
        d = b - d;

        if (fabs(d) < FPMIN) d = FPMIN;

        c = b - 1.0 / c;

        if (fabs(c) < FPMIN) c = FPMIN;

        d = 1.0 / d;

        del = c * d;

        h = del * h;

        if (d < 0.0) isign = -isign;

        if (fabs(del - 1.0) < EPS) break;
    }

    if (i > MAXIT) cout << "x too large in bessjy; try asymptotic expansion" << endl;

    rjl = isign * FPMIN;

    rjpl = h * rjl;

    rjl1 = rjl;

    rjp1 = rjpl;

    fact = xnu * xi;

	for (l = nl;l >= 1;l--)
    {
        rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
    }

    if (rjl == 0.0) rjl = EPS;

    f = rjpl / rjl;

    if (x < XMIN)
    {
        x2 = 0.5 * x;
        pimu = PI * xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
        d = -log(x2);
        e = xmu * d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e) / e);
		// call BESCHB
        BESCHB(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = 2.0 / PI * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = exp(e);
        p = e / (gampl * PI);
        q = 1.0 / (e * PI * gammi);
        pimu2 = 0.5 * pimu;
        fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2) / pimu2);
        r = PI * pimu2 * fact3 * fact3;
        c = 1.0;
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;

		for (i = 1;i <= MAXIT;i++)
        {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * (ff + r * q);
            sum += del;
            del1 = c * p - i * del;
            sum1 += del1;

            if (fabs(del) < (1.0 + fabs(sum))*EPS) break;
        }

        if (i > MAXIT) cout << "bessy series failed to converge";

        rymu = -sum;

        ry1 = -sum1 * xi2;

        rymup = xmu * xi * rymu - ry1;

        rjmu = w / (rymup - f * rymu);
    }

	else
    {
        a = 0.25 - xmu2;
        p = -0.5 * xi;
        q = 1.0;
        br = 2.0 * x;
        bi = 2.0;
        fact = a * xi / (p * p + q * q);
        cr = br + q * fact;
        ci = bi + p * fact;
        den = br * br + bi * bi;
        dr = br / den;
        di = -bi / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;

		for (i = 2;i <= MAXIT;i++)
        {
            a += 2 * (i - 1);
            bi += 2.0;
            dr = a * dr + br;
            di = a * di + bi;

            if (fabs(dr) + fabs(di) < FPMIN) dr = FPMIN;

            fact = a / (cr * cr + ci * ci);

            cr = br + cr * fact;

            ci = bi - ci * fact;

            if (fabs(cr) + fabs(ci) < FPMIN) cr = FPMIN;

            den = dr * dr + di * di;

            dr /= den;

            di /= -den;

            dlr = cr * dr - ci * di;

            dli = cr * di + ci * dr;

            temp = p * dlr - q * dli;

            q = p * dli + q * dlr;

            p = temp;

            if (fabs(dlr - 1.0) + fabs(dli) < EPS) break;
        }

		if (i > MAXIT) cout << "cf2 failed in bessjy";

        gam = (p - f) / q;

        rjmu = sqrt(w / ((p - f) * gam + q));

		if(rjl >=0 ) rjmu = fabs(rjmu);
		else rjmu = -fabs(rjmu);

        rymu = rjmu * gam;

        rymup = rymu * (p + q / gam);

        ry1 = xmu * xi * rymu - rymup;
    }

    fact = rjmu / rjl;

    *rj = rjl1 * fact;
    *rjp = rjp1 * fact;

    for (i = 1;i <= nl;i++)
    {
        rytemp = (xmu + i) * xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
    }

    *ry = rymu;

    *ryp = xnu * xi * rymu - ry1;
}

int Mathzone::IMAX(int a, int b)
{
    if (a > b) return a;
    else return b;
}

void Mathzone::BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	const int NUSE1 = 7;
	const int NUSE2 = 8;
    double xx;
    static double c1[] = {   -1.142022680371168e0, 6.5165112670737e-3,
                             3.087090173086e-4, -3.4706269649e-6,
                             6.9437664e-9, 3.67795e-11, -1.356e-13
                         };
    static double c2[] = {   1.843740587300905e0, -7.68528408447867e-2,
                             1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8,
                             2.423096e-10, -1.702e-13, -1.49e-15
                         };
    xx = 8.0 * x * x - 1.0; //Multiply x by 2 to make range be .1 to 1,and then apply transformation for evaluating even     Chebyshev series.
    *gam1 = CHEBEV(-1.0, 1.0, c1, NUSE1, xx);
    *gam2 = CHEBEV(-1.0, 1.0, c2, NUSE2, xx);
    *gampl = *gam2 - x * (*gam1);
    *gammi = *gam2 + x * (*gam1);
}

double Mathzone::CHEBEV(double a, double b, double c[], int m, double x)
{
    double d = 0.0, dd = 0.0, sv, y, y2;
    int j;

    if ((x - a)*(x - b) > 0.0) cout << "x not in range in routine chebev" << endl;

    y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));

    for (j = m - 1;j >= 1;j--)
    {
        sv = d;
        d = y2 * d - dd + c[j];
        dd = sv;
    }

    return y*d - dd + 0.5*c[0];
}


double Mathzone::Spherical_Bessel_7(const int n, const double &x)
{
	if(x==0)
	{
		if(n!=0) return 0;
		if(n==0) return 1;
	}
	double order, rj, rjp, ry, ryp;

	if (n < 0 || x <= 0.0) 
	{
		exit(0);
//		WARNING_QUIT("Mathzone::Spherical_Bessel_7","bad arguments in sphbes");
	}

	order = n + 0.5;

	// call BESSSJY
	BESSJY(x, order, &rj, &ry, &rjp, &ryp);

	const double RTPIO2=1.2533141; // mohan remian question: here is INT????????????

	const double factor = RTPIO2 / sqrt(x);

	return factor*rj;
}


void Mathzone::Spherical_Bessel_Roots
(
	const int &num,
	const int &l,	
	const double &epsilon,
	double* eigenvalue,
	const double &rcut
)
{
//	TITLE("Mathzone","Spherical_Bessel_Roots");
//	if(num<=0) WARNING_QUIT("Mathzone::Spherical_Bessel_Roots","num<=0");

	const double max = 2*PI + (num + (l+0.5)/2 + 0.75)*PI/2 + 
		sqrt((num + (l+0.5)/2+0.75)*(num + (l+0.5)/2+0.75)*PI*PI/4-(l+0.5)*(l+0.5)/2);

	// magic number !!
	// only need to > 1
	const int msh = 10 * num;
//	cout<<"\n msh = "<<msh;

	// delta don't need to be small,
	// it only needs to make sure can find the eigenstates
	const double delta = max / static_cast<double>(msh);

//	cout<<"\n delta = "<<delta;

	double *r = new double[msh];
	for(int i=0; i<msh; i++) 
	{
		r[i] = i*delta;
	}
//	cout<<"\n rmax(r) = "<<r[msh-1];
//	cout<<"\n Enumber = "<<num;

	double *jl = new double[msh];	
	Mathzone::Spherical_Bessel(msh, r, 1, l, jl);

	int n=0;
	for(int i=0; i<msh && n<num; i++)
	{
		if(jl[i]*jl[i+1] < 0.0)
		{
			double y_1 = jl[i];
			double y_2 = jl[i+1];
			double x_1 = r[i];
			double x_2 = r[i+1];
			double acc = abs( y_2 - y_1 );
			while(acc > epsilon)
			{
				double *rad = new double[100];
				double *jl_new = new double[100];

				// if not enough accurate, divide again.
				const double new_delta = (x_2 - x_1)/99.0;
				for(int j=0;j<100;j++)
				{
					rad[j] = x_1 + j*new_delta;
				}
				Mathzone::Spherical_Bessel(100,rad,1,l,jl_new);
				
				int j=0;
				for(;j<100;j++)
				{	
					if(jl_new[j]*jl_new[j+1]<0)break;
				}

				x_1 = rad[j];
				x_2 = rad[j+1];
				y_1 = jl_new[j];
				y_2 = jl_new[j+1];
				acc = abs( y_2 - y_1 );
				delete[] rad;
				delete[] jl_new;
			}	
			eigenvalue[n]=(x_2 + x_1)*0.5/rcut;
//			cout << "\n e="<< eigenvalue[n] << " n="<<n;
			n++;
		}
	}

//	cout << "\n emax="<< eigenvalue[num-1] << " num="<<num;
//	cout << "\n rmax( rcut * r ) = " << r[msh-1]/rcut << "(a.u.)" << endl;

	delete[] r;
	delete[] jl;
	return;
}

// from sph_bes.f90
void Mathzone::Spherical_Bessel
(
	const int &msh,	//number of grid points
	const double *r,//radial grid
	const double &q,	//
	const int &l,	//angular momentum
	double *jl		//jl(1:msh) = j_l(q*r(i)),spherical bessel function
)
{
//	timer::tick("Mathzone","Spherical_Bessel");
	double x1;
	int i, ir, ir0;

	if(l>=7)
	{
		for(int ir=0; ir<msh; ir++)
		{
			x1 = q * r[ir];
			jl[ir] = Spherical_Bessel_7(l, x1);
		}
		return;
	}
	
	if (abs(q) < 1.0e-8)
	{
		if (l == -1)
		{
			cout << "\n sph_bes, j_{-1}(0) ????";
		}
		else if (l == 0)
		{
			for (i = 0;i < msh;i++)
			{
				jl[i] = 1.0;
			}
		}
		else
		{
			for (i = 0;i < msh;i++)
			{
				jl[i] = 0.0;
			}
		}
	}
	else
	{
		if (abs(q * r [0]) > 1.0e-8)
		{
			ir0 = 0;//mohan modify 2007-10-13
		}
		else
		{
			if (l == -1)
			{
				cout << "\n sph_bes, j_{-1}(0) ?//?";
			}
			else if (l == 0)
			{
				jl [0] = 1.0;//mohan modify 2007-10-13
			}
			else
			{
				jl [0] = 0.0;//mohan modify 2007-10-13
			}
			ir0 = 1;//mohan modify 2007-10-13
		}
		if (l == - 1)
		{
			for (ir = ir0;ir < msh; ir++)
			{
				x1 = q * r[ir];
				jl [ir] = cos(x1) / x1;
			}
		}
		else if (l == 0)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				x1 = q * r[ir];
				jl [ir] = sin(x1) / x1;
			}
		}
		else if (l == 1)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				x1 = q * r[ir];
				const double sinx = sin(x1);
				const double cosx = cos(x1);
				jl [ir] = (sinx / x1 - cosx) / x1;
			}
		}
		else if (l == 2)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				const double x1 = q * r[ir];
				const double sinx = sin(x1);
				const double cosx = cos(x1);
				jl [ir] = ((3.0 / x1  - x1) * sinx
				           - 3.0 * cosx) / (x1 * x1);
			}
		}
		else if (l == 3)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				x1 = q * r[ir];
				jl [ir] = (sin(x1) * (15.0 / x1 - 6.0 * x1) +
				           cos(x1) * (x1 * x1 - 15.0)) / pow(x1, 3);//mohan modify 2007-10-13
			}
		}
		else if (l == 4)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				const double x1 = q * r[ir];
				const double x2 = x1 * x1;
				const double x3 = x1 * x2;
				const double x4 = x1 * x3;
				const double x5 = x1 * x4;
				jl [ir] = (sin(x1) * (105.0 - 45.0 * x2 + x4) +
				           cos(x1)  * (10.0 * x3 - 105.0 * x1)) / x5;   // mohan modify 2007-10-13
			}
		}
		else if (l == 5)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				x1 = q * r[ir];

				if (x1 < 0.14)
				{
					jl[ir] = 0;//mohan add 2007-10-15
				}
				else
				{
					double cx1 = cos(x1);
					double sx1 = sin(x1);
					jl [ir] = (-cx1 -
					           (945.0 * cx1) / pow(x1, 4) +
					           (105.0 * cx1) / (x1 * x1)  +
					           (945.0 * sx1) / pow(x1, 5) -
					           (420.0 * sx1) / pow(x1, 3) +
					           (15.0 * sx1) / x1) / x1;

				}
			}
		}
		else if (l == 6)
		{
			for (ir = ir0;ir < msh;ir++)
			{
				x1 = q * r[ir];

				if (x1 < 0.29)
				{
					jl[ir] = 0;//mohan add 2007-10-15
				}
				else
				{
					double cx1 = cos(x1);
					double sx1 = sin(x1);
					jl [ir] = ((-10395.0 * cx1) / pow(x1, 5) +
					           (1260.0 * cx1) / pow(x1, 3) -
					           (21.0 * cx1) / x1 - sx1 +
					           (10395.0 * sx1) / pow(x1, 6) -
					           (4725.0 * sx1) / pow(x1, 4) +
					           (210.0 * sx1) / (x1 * x1)) / x1;
				}
			}
		}//mohan modify 2007-11-20 reduce cos , sin , q*r[ir] times;
		else
		{
			cout << "\n error in sph_bes, l out of {-1 ... 6},l = " << l ;
		}
	}
//	printr1_d(ofs,"jl[] :",jl,msh);
//	timer::tick("Mathzone","Spherical_Bessel");
	return;
} 

void Mathzone::Simpson_Integral
(
	const int mesh,
	const double *func,
	const double *rab,
	double &asum
)
{
  /*     simpson's rule integration. On input:
  !      mesh = mhe number of grid points (should be odd)
  !      func(i)= function to be integrated
  !      rab(i) = r(i) * dr(i)/di * di
  !      For the logarithmic grid not including r=0 :
  !      r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !      For the logarithmic grid including r=0 :
  !      r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !      Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !      where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  */
	//  simpson's rule integrator for function stored on the
	//  radial logarithmic mesh
	//	routine assumes that mesh is an odd number so run check
	if (mesh % 2 == 0)
	{
		cout << "\n error in subroutine simpson ";
		cout << "\n routine assumes mesh is odd but mesh = "
		<< mesh << endl;
		return;
	}

	asum = 0.00;
	const double r12 = 1.00 / 12.00;
	double f3 = func [0] * rab [0] * r12;
	for (int i = 1;i < mesh;i += 2)
	{
		const double f1 = f3;
		const double f2 = func [i] * rab [i] * r12;
		f3 = func [i + 1] * rab [i + 1] * r12;
		asum += 4.00 * f1 + 16.00 * f2 + 4.00 * f3;
	}
	return;
}// end subroutine simpson


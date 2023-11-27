#include "math_sphbes.h"
#include "timer.h"
#include "constants.h"
#include <algorithm>

#include <cassert>

namespace ModuleBase
{

Sphbes::Sphbes(){}
Sphbes::~Sphbes(){}

void Sphbes::BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{
    const int XMIN = 2.0;
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
		std::cout << "Sphbes::BESSJY, bad arguments" << std::endl;
        //ModuleBase::WARNING_QUIT("Sphbes::BESSJY","bad arguments");
		exit(0); // mohan add 2021-05-06
    }


    nl = (x < XMIN ? (int)(xnu + 0.5) : IMAX(0, (int)(xnu - x + 1.5)));
    const double xmu = xnu - nl;
    const double xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    w = xi2 / ModuleBase::PI;
    isign = 1;
    h = xnu * xi;

    if (h < FPMIN)
	{
		h = FPMIN;
	}

    b = xi2 * xnu;

    d = 0.0;

    c = h;

    for (i = 1;i <= MAXIT;i++)
    {
        b += xi2;
        d = b - d;

        if (std::fabs(d) < FPMIN) d = FPMIN;

        c = b - 1.0 / c;

        if (std::fabs(c) < FPMIN) c = FPMIN;

        d = 1.0 / d;

        del = c * d;

        h = del * h;

        if (d < 0.0) isign = -isign;

        if (std::fabs(del - 1.0) < EPS) break;
    }

    if (i > MAXIT)
	{
		std::cout << "x too large in bessjy; try asymptotic expansion" << std::endl;
	}

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

    if (rjl == 0.0)
	{
		rjl = EPS;
	}

    f = rjpl / rjl;

    if (x < XMIN)
    {
        x2 = 0.5 * x;
        pimu = ModuleBase::PI * xmu;
        fact = (std::fabs(pimu) < EPS ? 1.0 : pimu / std::sin(pimu));
        d = -log(x2);
        e = xmu * d;
        fact2 = (std::fabs(e) < EPS ? 1.0 : std::sinh(e) / e);
        // call BESCHB
        BESCHB(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = 2.0 / ModuleBase::PI * fact * (gam1 * std::cosh(e) + gam2 * fact2 * d);
        e = std::exp(e);
        p = e / (gampl * ModuleBase::PI);
        q = 1.0 / (e * ModuleBase::PI * gammi);
        pimu2 = 0.5 * pimu;
        fact3 = (std::fabs(pimu2) < EPS ? 1.0 : std::sin(pimu2) / pimu2);
        r = ModuleBase::PI * pimu2 * fact3 * fact3;
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

            if (std::fabs(del) < (1.0 + std::fabs(sum))*EPS) break;
        }

        if (i > MAXIT) std::cout << "bessy series failed to converge";

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

            if (std::fabs(dr) + std::fabs(di) < FPMIN) dr = FPMIN;

            fact = a / (cr * cr + ci * ci);

            cr = br + cr * fact;

            ci = bi - ci * fact;

            if (std::fabs(cr) + std::fabs(ci) < FPMIN) cr = FPMIN;

            den = dr * dr + di * di;

            dr /= den;

            di /= -den;

            dlr = cr * dr - ci * di;

            dli = cr * di + ci * dr;

            temp = p * dlr - q * dli;

            q = p * dli + q * dlr;

            p = temp;

            if (std::fabs(dlr - 1.0) + std::fabs(dli) < EPS) break;
        }

        if (i > MAXIT) std::cout << "cf2 failed in bessjy";

        gam = (p - f) / q;

        rjmu = std::sqrt(w / ((p - f) * gam + q));

        if (rjl >=0 ) rjmu = std::fabs(rjmu);
        else rjmu = -std::fabs(rjmu);

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


int Sphbes::IMAX(int a, int b)
{
    if (a > b) return a;
    else return b;
}


void Sphbes::BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi)
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

double Sphbes::CHEBEV(double a, double b, double c[], int m, double x)
{
    double d = 0.0;
	double dd = 0.0;
	double sv = 0.0;
	double y = 0.0;
	double y2 = 0.0;
    int j=0;

    if ((x - a)*(x - b) > 0.0)
	{
		std::cout << "x not in range in routine chebev" << std::endl;
	}

    y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));

    for (j = m - 1;j >= 1;j--)
    {
        sv = d;
        d = y2 * d - dd + c[j];
        dd = sv;
    }

    return y*d - dd + 0.5*c[0];
}


double Sphbes::Spherical_Bessel_7(const int n, const double &x)
{
    if (x==0)
    {
        if (n!=0) return 0;
        if (n==0) return 1;
    }
    double order, rj, rjp, ry, ryp;

    if (n < 0 || x <= 0.0)
    {
		std::cout << "Spherical_Bessel_7, bad arguments in sphbes" << std::endl;
        //ModuleBase::WARNING_QUIT("Sphbes::Spherical_Bessel_7","bad arguments in sphbes");
		exit(0);
    }

    order = n + 0.5;

    // call BESSSJY
    BESSJY(x, order, &rj, &ry, &rjp, &ryp);

    const double RTPIO2=1.2533141;

    const double factor = RTPIO2 / std::sqrt(x);

    return factor*rj;
}


void Sphbes::Spherical_Bessel_Roots
(
    const int &num,
    const int &l,
    const double &epsilon,
    double* eigenvalue,
    const double &rcut
)
{
    //ModuleBase::TITLE("Sphbes","Spherical_Bessel_Roots");
    if (num<=0)
	{
		std::cout << "Spherical_Bessel_Roots, num<=0" << std::endl;
		//ModuleBase::WARNING_QUIT("Sphbes::Spherical_Bessel_Roots","num<=0");
		exit(0);
	}
    if (rcut<=0.0)
	{
		std::cout << "Spherical_Bessel_Roots, rcut<=0" << std::endl;
		//ModuleBase::WARNING_QUIT("Sphbes::Spherical_Bessel_Roots","rcut<=0.0");
		exit(0);
	}

    double min = 0.0;
    double max = 2*ModuleBase::PI + (num + (l+0.5)/2 + 0.75)*ModuleBase::PI/2 +
                 std::sqrt((num + (l+0.5)/2+0.75)*(num + (l+0.5)/2+0.75)*ModuleBase::PI*ModuleBase::PI/4-(l+0.5)*(l+0.5)/2);

    // magic number !!
    // guess : only need to > 1
    const int msh = 10 * num;
//	std::cout<<"\n msh = "<<msh;

    // delta don't need to be small,
    // it only needs to make sure can find the eigenstates
    const double delta = (max - min) / static_cast<double>(msh);
//	std::cout<<"\n delta = "<<delta;

    double *r = new double[msh];
    for (int i=0; i<msh; i++)
    {
        r[i] = i*delta;
    }
    double *jl = new double[msh];

    Sphbes::Spherical_Bessel(msh, r, 1, l, jl);

    int n=0;
    for (int i=0; i<msh-1 && n<num; i++)
    {
        if (jl[i]*jl[i+1] < 0.0)
        {
            double y_1 = jl[i];
            double y_2 = jl[i+1];
            double x_1 = r[i];
            double x_2 = r[i+1];
            double acc = std::fabs(y_2 - y_1);

            const int grid=100;
            double *rad = new double[grid];
            double *jl_new = new double[grid];
            while (acc > epsilon)
            {
                // if not enough accurate, divide again.
                const double delta2 = (x_2 - x_1)/(grid-1);
                for (int j=0;j<grid;j++)
                {
                    rad[j] = x_1 + j*delta2;
                }
                Sphbes::Spherical_Bessel(grid,rad,1,l,jl_new);

                int j=0;
                for (;j<grid-1;j++)
                {
                    if (jl_new[j]*jl_new[j+1]<0)break;
                }

                x_1 = rad[j];
                x_2 = rad[j+1];
                y_1 = jl_new[j];
                y_2 = jl_new[j+1];
                acc = std::fabs( y_2 - y_1 );
            }
            delete[] rad;
            delete[] jl_new;
            eigenvalue[n]=(x_2 + x_1)*0.5/rcut;
            n++;
        }
    }
    delete[] r;
    delete[] jl;
}


void Sphbes::Spherical_Bessel
(
    const int &msh,	 // number of grid points
    const double *r, // radial grid
    const double &q, // wave std::vector
    const int &l,	 // angular momentum
    double *jl		 // jl(1:msh) = j_l(q*r(i)),spherical bessel function
)
{
    ModuleBase::timer::tick("Sphbes","Spherical_Bessel");
    double x1=0.0;

    int i=0;
	int ir=0;
	int ir0=0;

    if (l>=7)
    {
        for (int ir=0; ir<msh; ir++)
        {
            x1 = q * r[ir];
            jl[ir] = Spherical_Bessel_7(l, x1);
        }
        return;
    }

    if (std::fabs(q) < 1.0e-8)
    {
        if (l == -1)
        {
            std::cout << "\n sph_bes, j_{-1}(0) ????";
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
        if (std::fabs(q * r [0]) > 1.0e-8)
        {
            ir0 = 0;//mohan modify 2007-10-13
        }
        else
        {
            if (l == -1)
            {
                std::cout << "\n sph_bes, j_{-1}(0) ?//?";
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
                jl [ir] = std::cos(x1) / x1;
            }
        }
        else if (l == 0)
        {
            for (ir = ir0;ir < msh;ir++)
            {
                x1 = q * r[ir];
                jl [ir] = std::sin(x1) / x1;
            }
        }
        else if (l == 1)
        {
            for (ir = ir0;ir < msh;ir++)
            {
                x1 = q * r[ir];
                const double sinx = std::sin(x1);
                const double cosx = std::cos(x1);
                jl [ir] = (sinx / x1 - cosx) / x1;
            }
        }
        else if (l == 2)
        {
            for (ir = ir0;ir < msh;ir++)
            {
                const double x1 = q * r[ir];
                const double sinx = std::sin(x1);
                const double cosx = std::cos(x1);
                jl [ir] = ((3.0 / x1  - x1) * sinx
                           - 3.0 * cosx) / (x1 * x1);
            }
        }
        else if (l == 3)
        {
            for (ir = ir0;ir < msh;ir++)
            {
                x1 = q * r[ir];
                jl [ir] = (std::sin(x1) * (15.0 / x1 - 6.0 * x1) +
                           std::cos(x1) * (x1 * x1 - 15.0)) / std::pow(x1, 3);//mohan modify 2007-10-13
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
                jl [ir] = (std::sin(x1) * (105.0 - 45.0 * x2 + x4) +
                           std::cos(x1)  * (10.0 * x3 - 105.0 * x1)) / x5;   // mohan modify 2007-10-13
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
                    double cx1 = std::cos(x1);
                    double sx1 = std::sin(x1);
                    jl [ir] = (-cx1 -
                               (945.0 * cx1) / std::pow(x1, 4) +
                               (105.0 * cx1) / (x1 * x1)  +
                               (945.0 * sx1) / std::pow(x1, 5) -
                               (420.0 * sx1) / std::pow(x1, 3) +
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
                    double cx1 = std::cos(x1);
                    double sx1 = std::sin(x1);
                    jl [ir] = ((-10395.0 * cx1) / std::pow(x1, 5) +
                               (1260.0 * cx1) / std::pow(x1, 3) -
                               (21.0 * cx1) / x1 - sx1 +
                               (10395.0 * sx1) / std::pow(x1, 6) -
                               (4725.0 * sx1) / std::pow(x1, 4) +
                               (210.0 * sx1) / (x1 * x1)) / x1;
                }
            }
        }//mohan modify 2007-11-20 reduce cos , sin , q*r[ir] times;
        else
        {
            std::cout << "\n error in sph_bes, l out of {-1 ... 6},l = " << l ;
			exit(0);
        }
    }

    ModuleBase::timer::tick("Sphbes","Spherical_Bessel");
    return;
}


void Sphbes::Spherical_Bessel
(
	const int &msh, //number of grid points
	const double *r,//radial grid
	const double &q,    //
	const int &l,   //angular momentum
	double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
	double *sjp
)
{
	ModuleBase::timer::tick("Sphbes","Spherical_Bessel");

	//calculate jlx first
	Spherical_Bessel (msh, r, q, l, sj);

	for (int ir = 0; ir < msh; ir++)
	{
		sjp[ir] = 1.0;
	}
	return;
}

void Sphbes::dSpherical_Bessel_dx
(
    const int &msh,	 // number of grid points
    const double *r, // radial grid
    const double &q, // wave std::vector
    const int &l,	 // angular momentum
    double *djl		 // jl(1:msh) = j_l(q*r(i)),spherical bessel function
)
{
    ModuleBase::timer::tick("Sphbes","dSpherical_Bessel_dq");
    if (l < 0 )
    {
		std::cout << "We temporarily only calculate derivative of l >= 0." << std::endl;
		exit(0);
    }

    double djl0 = 0;
    if(l == 1)
    {
        djl0 = 1.0/3.0;
    }

    if(l == 0 )
    {
        for (int ir = 0;ir < msh; ir++)
        {
            double x1 = q * r[ir];
            if(x1 < 1e-8)
            {
                djl[ir] = djl0;
            }
            else
            {
                djl[ir] = (x1 * std::cos(x1) - std::sin(x1)) / (x1*x1);
            }

        }
    }
    else
    {
        double *jl = new double [msh];
        Spherical_Bessel (msh, r, q, l-1, jl);
        Spherical_Bessel (msh, r, q, l, djl);
        for (int ir = 0;ir < msh; ir++)
        {
            double x1 = q * r[ir];
            if(x1 < 1e-8)
            {
                djl[ir] = djl0;
            }
            else
            {
                djl[ir] = jl[ir] - double(l+1)/x1 * djl[ir];
            }
        }
        delete[] jl;
    }
    ModuleBase::timer::tick("Sphbes","dSpherical_Bessel_dq");
    return;
}

double Sphbes::_sphbesj_ascending_recurrence(int l, double x) {

    // should be used when x > l && l > 0

    double invx = 1.0 / x;
    double j0 = std::sin(x) * invx;
    double j1 = ( j0 - std::cos(x) ) * invx;

    double jl;
    for (int i = 2; i <= l; ++i) {
        jl = (2*i-1) * invx * j1 - j0;
        j0 = j1;
        j1 = jl;
    }

    return j1; // at the end of the loop j1 == jl
}

double Sphbes::_sphbesj_series(int l, double x) {

    // should be used when x < l

    // the absolute ratio between the k-th and (k-1)-th terms is x^2 / k / (4l+4k+2) (k >= 1)
    // terms are guaranteed to be monotonically decreasing from the beginning for x < sqrt(4l+6).
    // terms are guaranteed to be monotonically decreasing from the (l/4)-th term for x < l

    double jl = 0.0;
    constexpr double eps = 1e-17; // series terminate when the k-th term is less than eps * jl

    // zeroth order term: x^l / (2l+1)!!
    int k = 0;
    double kth_term = 1.0;
    for (int i = 1; i <= l; ++i) {
        kth_term *= x / (2 * i + 1);
    }

    double x_sqr_half = 0.5 * x * x;
    do {
        jl += kth_term;
        k += 1;
        kth_term *= -x_sqr_half / ( k * (2*(l+k)+1) );
    } while ( std::abs(kth_term) > std::abs(eps * jl) );

    return jl;
}

double Sphbes::sphbesj(const int l, const double x)
{
    assert( l >= 0 );
    assert( x >= 0 );

    // j_l(0)
    if ( x == 0 )
    {
        return l ? 0.0 : 1.0;
    }

    if ( x < l )
    {
        return _sphbesj_series(l, x);
    }
    else
    {
        double invx = 1.0 / x;
        switch (l)
        {
          case 0:
            return std::sin(x) * invx;
          case 1:
            return ( std::sin(x) * invx - std::cos(x) ) * invx;
            // NOTE: the following explicit expressions are not necessarily faster than ascending recurrence,
            // but we keep them just in case we need them in the future.
          //case 2:
          //  return ( (3.0 * invx  - x) * std::sin(x) - 3.0 * std::cos(x) ) * (invx * invx);
          //case 3:
          //  return ( std::sin(x) * (15.0 * invx - 6.0 * x) + std::cos(x) * (x * x - 15.0) ) * std::pow(invx, 3);
          //case 4:
          //  return ( std::sin(x) * (std::pow(x,3) - 45.0 * x + 105.0 * invx)
          //          + std::cos(x) * (10.0 * x * x - 105.0) ) * std::pow(invx, 4);
          //case 5:
          //  return ( std::sin(x) * (15.0 * std::pow(x,3) - 420.0 * x + 945.0 * invx)
          //          + std::cos(x) * (-std::pow(x, 4) + 105.0 * x * x - 945.0) ) * std::pow(invx, 5);
          //case 6:
          //  return ( std::sin(x) * (-std::pow(x, 5) + 210.0 * std::pow(x, 3) - 4725.0 * x + 10395.0 * invx)
          //          + std::cos(x) * (-21.0 * std::pow(x, 4) + 1260.0 * x * x - 10395.0) ) * std::pow(invx, 6);
          default:
            return _sphbesj_ascending_recurrence(l, x);
        }
    }
}

double Sphbes::dsphbesj(const int l, const double x)
{
    assert( l >= 0 );
    assert( x >= 0 );
    return l == 0 ? -sphbesj(1, x) : ( l * sphbesj(l - 1, x) - (l + 1) * sphbesj(l + 1, x) ) / (2 * l + 1);
}

void Sphbes::sphbesj(const int n,
                     const double* const r,
                     const double q,
                     const int l,
                     double* const jl)
{
    for (int i = 0; i != n; ++i)
    {
        jl[i] = Sphbes::sphbesj(l, q * r[i]);
    }
}

void Sphbes::dsphbesj(const int n,
             const double* const r,
             const double q,
             const int l,
             double* const djl)
{
    for (int i = 0; i != n; ++i)
    {
        djl[i] = Sphbes::dsphbesj(l, q * r[i]);
    }
}

void Sphbes::sphbes_zeros(const int l, const int n, double* const zeros)
{
    assert( n > 0 );
    assert( l >= 0 );

    // The zeros of j_l and j_{l-1} are interlaced;
    // So do the zeros of j_l and j_{l-2}.
    // This property enables us to use bracketing method recursively
    // to find all zeros of j_l from the zeros of j_0.

    // if l is odd , j_0 --> j_1 --> j_3 --> j_5 --> ...
    // if l is even, j_0 --> j_2 --> j_4 --> j_6 --> ...

    int nz = n + (l+1)/2; // number of effective zeros in buffer
    double* buffer = new double[nz];

    // zeros of j_0 = sin(x)/x is just n*pi
    double PI = std::acos(-1.0);
    for (int i = 0; i < nz; i++)
    {
        buffer[i] = (i+1) * PI;
    }

    int ll = 1;
    auto jl = [&ll] (double x) { return sphbesj(ll, x); };

    if (l % 2 == 1)
    {
        for (int i = 0; i < nz-1; i++)
        {
            buffer[i] = illinois(jl, buffer[i], buffer[i+1], 1e-15, 50);
        }
        --nz;
    }

    for (ll = 2 + l%2; ll <= l; ll += 2, --nz)
    {
        for (int i = 0; i < nz-1; i++)
        {
            buffer[i] = illinois(jl, buffer[i], buffer[i+1], 1e-15, 50);
        }
    }

    std::copy(buffer, buffer + n, zeros);
    delete[] buffer;
}

double Sphbes::illinois(std::function<double(double)> func, double x0, double x1, const double tol, const int max_iter)
{
    assert(tol > 0.0 && max_iter > 0);

    double f0 = func(x0);
    double f1 = func(x1);
    assert(f0 * f1 <= 0);

    if (std::abs(f0) < std::abs(f1)) {
        std::swap(x0, x1);
        std::swap(f0, f1);
    }

    int iter = 0;
    double x, f;
    while (++iter <= max_iter && std::abs(f1) > tol)
    {
        // regula falsi
        x = (x0 * f1 - x1 * f0) / (f1 - f0);
        f = func(x);

        // Illinois anti-stalling variant
        if (f * f1 < 0)
        {
            x0 = x1;
            f0 = f1;
        }
        else
        {
            f0 *= 0.5;
        }
        x1 = x;
        f1 = f;
    }

    if (iter > max_iter)
    {
        std::cout << "Maximum number of iterations reached in illinois." << std::endl;
    }

    return x1;
}

}

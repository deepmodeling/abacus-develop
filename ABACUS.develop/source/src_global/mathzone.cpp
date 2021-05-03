#include "../src_pw/global.h"
#include "mathzone.h"
#include "../src_parallel/parallel_reduce.h"

Mathzone::Mathzone()
{
}

Mathzone::~Mathzone()
{
}

void Mathzone::norm_pw(complex<double> *u,const int n)
{
    if (n < 1) WARNING_QUIT("Mathzone::norm","n < 1");
    complex<double> rrc = ZERO;//complex < double> (0.0,0.0);

    for (int i = 0;i < n;i++) rrc += conj(u[i]) * u[i];

    Parallel_Reduce::reduce_complex_double_pool( rrc );

    double rr = rrc.real();

    if (rr <= 1.e-20) WARNING_QUIT("Mathzone::norm","rr <= 1.e-20");

    rr = 1.0/sqrt(rr);

    if (rr==1.0) return;
    else
    {
        for (int i = 0;i < n;i++) u[i] *= rr;
    }
}

void Mathzone::Polynomial_Interpolation
(
    const realArray &table,
    const int &dim1,
    const int &dim2,
    realArray &y,
    const int &dim_y,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
    timer::tick("Mathzone","Poly_Interpo_1");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
	if(iq>=table_length-4)
	{
		cout << "\n x = " << x;
		cout << "\n iq = " << iq << " table_length = " << table_length << endl;
	}	
    assert(iq < table_length-4);

    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    y(dim1, dim2, dim_y)=
        table(dim1, dim2, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, iq+3) * x1 * x2 * x0 / 6.0 ;

    timer::tick("Mathzone","Poly_Interpo_1");
    return;
}

double Mathzone::Polynomial_Interpolation
(
    const realArray &table,
    const int &dim1,
    const int &dim2,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	timer::tick("Mathzone","Poly_Interpo_2");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
    
	if(iq>table_length-4)
	{
		cout << "\n x = " << x;
		cout << "\n table_interval = " << table_interval;
		cout << "\n iq=" << iq << " table_length = " << table_length << endl;
        cout << "\n Not enough space allocated for radial FFT: try restarting with a larger cell_factor" << endl; //LiuXh add 20180619
        cout << "\n Now cell_factor is: " << ppcell.cell_factor << endl; //LiuXh add 20180619
	}
	assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
        table(dim1, dim2, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, iq+3) * x1 * x2 * x0 / 6.0 ;

//	timer::tick("Mathzone","Poly_Interpo_2");
    return y;
}

double Mathzone::Polynomial_Interpolation            // pengfei Li 2018-3-23
(
    const realArray &table,
    const int &dim1,
    const int &dim2,
	const int &dim3,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	timer::tick("Mathzone","Poly_Interpo_3");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
    
	if(iq>table_length-4)
	{
		cout << "\n x = " << x;
		cout << "\n table_interval = " << table_interval;
		cout << "\n iq=" << iq << " table_length = " << table_length << endl;
	}
	assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
        table(dim1, dim2, dim3, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, dim3, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, dim3, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, dim3, iq+3) * x1 * x2 * x0 / 6.0 ;

//	timer::tick("Mathzone","Poly_Interpo_3");
    return y;
}

double Mathzone::Polynomial_Interpolation
(
    const double *table,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	assert(table_interval>0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
//	if(iq >= table_length-4)
//		cout << "\n iq = " << iq << " table_length = " << table_length;
  
   assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;

    /*
    const double y=
    	table[iq]   * x1 * x2 * x3 / 6.0 +
    	table[iq+1] * x0 * x2 * x3 / 2.0 -
    	table[iq+2] * x1 * x0 * x3 / 2.0 +
    	table[iq+3] * x1 * x2 * x0 / 6.0 ;
    	*/

    return x1*x2*(table[iq]*x3+table[iq+3]*x0)/6.0
         + x0*x3*(table[iq+1]*x2-table[iq+2]*x1)/2.0;
}

double Mathzone::Polynomial_Interpolation_xy
(
    const double *xpoint,
    const double *ypoint,
    const int table_length,
    const double &x             // input value
)
{
    int position = -1;

    if (x < xpoint[0])
    {
        return ypoint[0];
    }
    // timer::tick("Mathzone","Poly_Inter_xy");

    for (int ik = 0; ik < table_length; ik++)
    {
        if (x < xpoint[ik])
        {
            break;
        }
        else
        {
            position ++;
        }
    }

    assert(position >= 0);
    assert(position <= table_length-1);

    if (position + 6 < table_length)
    {
        double dx1, dx2, dx3, dx4, dx5, dx6;
        dx1 = x - xpoint[position];
        dx2 = x - xpoint[position+1];
        dx3 = x - xpoint[position+2];
        dx4 = x - xpoint[position+3];
        dx5 = x - xpoint[position+4];
        dx6 = x - xpoint[position+5];


        double x12, x13, x14, x15, x16, x23, x24, x25, x26, x34, x35, x36, x45, x46, x56;
        x12 = xpoint[position] - xpoint[position+1];
        x13 = xpoint[position] - xpoint[position+2];
        x14 = xpoint[position] - xpoint[position+3];
        x15 = xpoint[position] - xpoint[position+4];
        x16 = xpoint[position] - xpoint[position+5];


        x23 = xpoint[position+1] - xpoint[position+2];
        x24 = xpoint[position+1] - xpoint[position+3];
        x25 = xpoint[position+1] - xpoint[position+4];
        x26 = xpoint[position+1] - xpoint[position+5];

        x34 = xpoint[position+2] - xpoint[position+3];
        x35 = xpoint[position+2] - xpoint[position+4];
        x36 = xpoint[position+2] - xpoint[position+5];

        x45 = xpoint[position+3] - xpoint[position+4];
        x46 = xpoint[position+3] - xpoint[position+5];

        x56 = xpoint[position+4] - xpoint[position+5];

        double part1, part2, part3, part4, part5, part6;
        part1 = dx2 * dx3 * dx4 * dx5 * dx6 / x12 / x13 / x14 / x15 / x16 * ypoint[position];
        part2 = dx1 * dx3 * dx4 * dx5 * dx6 / (-x12) / x23 / x24 / x25 / x26 * ypoint[position+1];
        part3 = dx1 * dx2 * dx4 * dx5 * dx6 / (-x13) / (-x23) / x34 / x35 / x36 * ypoint[position+2];
        part4 = dx1 * dx2 * dx3 * dx5 * dx6 / (-x14) / (-x24) / (-x34) / x45 / x46 * ypoint[position+3];
        part5 = dx1 * dx2 * dx3 * dx4 * dx6 / (-x15) / (-x25) / (-x35) / (-x45) / x56 * ypoint[position+4];
        part6 = dx1 * dx2 * dx3 * dx4 * dx5 / (-x16) / (-x26) / (-x36) / (-x46) / (-x56) * ypoint[position+5];

        // 	timer::tick("Mathzone","Poly_Inter_xy");
        return part1 + part2 + part3 + part4 + part5 + part6;
    }
    else
    {
        // 	timer::tick("Mathzone","Poly_Inter_xy");
        return ypoint[position];
    }
}


void Mathzone::BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
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
        WARNING_QUIT("Mathzone::BESSJY","bad arguments");
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

        if (rjl >=0 ) rjmu = fabs(rjmu);
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
    if (x==0)
    {
        if (n!=0) return 0;
        if (n==0) return 1;
    }
    double order, rj, rjp, ry, ryp;

    if (n < 0 || x <= 0.0)
    {
        WARNING_QUIT("Mathzone::Spherical_Bessel_7","bad arguments in sphbes");
    }

    order = n + 0.5;

    // call BESSSJY
    BESSJY(x, order, &rj, &ry, &rjp, &ryp);

    const int RTPIO2=1.2533141;

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
    TITLE("Mathzone","Spherical_Bessel_Roots");
    if (num<=0) WARNING_QUIT("Mathzone::Spherical_Bessel_Roots","num<=0");
    if (rcut<=0.0) WARNING_QUIT("Mathzone::Spherical_Bessel_Roots","rcut<=0.0");

    double min = 0.0;
    double max = 2*PI + (num + (l+0.5)/2 + 0.75)*PI/2 +
                 sqrt((num + (l+0.5)/2+0.75)*(num + (l+0.5)/2+0.75)*PI*PI/4-(l+0.5)*(l+0.5)/2);

    // magic number !!
    // guess : only need to > 1
    const int msh = 10 * num;
//	cout<<"\n msh = "<<msh;

    // delta don't need to be small,
    // it only needs to make sure can find the eigenstates
    const double delta = (max - min) / static_cast<double>(msh);
//	cout<<"\n delta = "<<delta;

    double *r = new double[msh];
    for (int i=0; i<msh; i++)
    {
        r[i] = i*delta;
    }
    double *jl = new double[msh];

    Mathzone::Spherical_Bessel(msh, r, 1, l, jl);

    int n=0;
    for (int i=0; i<msh && n<num; i++)
    {
        if (jl[i]*jl[i+1] < 0.0)
        {
            double y_1 = jl[i];
            double y_2 = jl[i+1];
            double x_1 = r[i];
            double x_2 = r[i+1];
            double acc = abs( y_2 - y_1 );
            while (acc > epsilon)
            {
                double *rad = new double[100];
                double *jl_new = new double[100];

                // if not enough accurate, divide again.
                const double delta2 = (x_2 - x_1)/99.0;
                for (int j=0;j<100;j++)
                {
                    rad[j] = x_1 + j*delta2;
                }
                Mathzone::Spherical_Bessel(100,rad,1,l,jl_new);

                int j=0;
                for (;j<100;j++)
                {
                    if (jl_new[j]*jl_new[j+1]<0)break;
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
            n++;
        }
    }
    delete[] r;
    delete[] jl;
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
    timer::tick("Mathzone","Spherical_Bessel");
    double x1;
    int i, ir, ir0;

    if (l>=7)
    {
        for (int ir=0; ir<msh; ir++)
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
    timer::tick("Mathzone","Spherical_Bessel");
    return;
}


void Mathzone::Spherical_Bessel
(           
	const int &msh, //number of grid points
	const double *r,//radial grid
	const double &q,    //
	const int &l,   //angular momentum
	double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
	double *sjp
)
{
	timer::tick("Mathzone","Spherical_Bessel");

	//calculate jlx first
	Spherical_Bessel (msh, r, q, l, sj);
	
	for (int ir = 0; ir < msh; ir++)
	{
		sjp[ir] = 1.0;
	}
	return;
}

void Mathzone::rlylm
(
    const int lmax, 	
    const double& x,				
    const double& y,
	const double& z, // g_cartesian_vec(x,y,z)
    double* rly 	 // output
)
{
//	TITLE("Ylm","rlylm");
	timer::tick("Mathzone","rlylm");

	assert(lmax >= 0);

	//get xy_dependence
	assert(lmax <= 19);
	
	double Am[20];
	double Bm[20];
	
	ZEROS(Am, 20);
	ZEROS(Bm, 20);
	
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
		ZEROS(zdep[il], 20);
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

	timer::tick("Mathzone", "rlylm");
	return;
}

void Mathzone::Ylm_Real2
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
void Mathzone::Ylm_Real
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
        WARNING("YLM_REAL","ng<1 or lmax2<1");
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
long double Mathzone::Fact(const int n)
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

int Mathzone::Semi_Fact(const int n)
{
    int semif = 1;
    for (int i=n; i>2; i -= 2)
    {
        semif *= i;
    }
    return semif;
}


void Mathzone::To_Polar_Coordinate
(
    const double &x_cartesian,
    const double &y_cartesian,
    const double &z_cartesian,
    double &r,
    double &theta,
    double &phi
)
{
//	OUT("x",x_cartesian);
//	OUT("y",y_cartesian);
//	OUT("z",z_cartesian);
    r = sqrt(x_cartesian*x_cartesian
             + y_cartesian*y_cartesian
             + z_cartesian*z_cartesian);
//----------------------------------------------------------
// Calculate theta
//----------------------------------------------------------
    double small = 1.0e-9;
    if (r < small)
    {
        theta = PI/2.0;
    }
    else
    {
        theta = acos(z_cartesian / r);
    }
//----------------------------------------------------------
// Calculate phi
//----------------------------------------------------------
    if (x_cartesian > small && y_cartesian > small)
    {
        phi = atan( y_cartesian / x_cartesian );
    }
    else if (  x_cartesian < small  )
    {
        phi = atan( y_cartesian / x_cartesian ) + PI;
    }
    else if ( x_cartesian > small && y_cartesian < -small)
    {
        phi = atan( y_cartesian/x_cartesian ) + TWO_PI;
    }
    else
    {
        phi = PI/2.0 * ( (y_cartesian >= 0.0) ? 1.0 : -1.0);
    }
//----------------------------------------------------------
// degress = radians / PI * 180
//----------------------------------------------------------
    theta = theta/PI*180;
    phi = phi/PI*180;

//	OUT("r",r);
//	OUT("theta(degree)",theta);
//	OUT("phi(degree)",phi);
    return;
}


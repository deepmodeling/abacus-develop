#include "bessel_basis.h"
#include "../src_parallel/parallel_common.h"

Bessel_Basis::Bessel_Basis()
{
	Dk = 0.0;
}

Bessel_Basis::~Bessel_Basis()
{
}


// the function is called in numerical_basis.
void Bessel_Basis::init( 
	const double &dk,
	const double &dr)
{
	TITLE("Bessel_Basis", "init");
	this->Dk = dk;

	// to make a table
	this->init_TableOne( SMOOTH, SIGMA, ECUT, RCUT, dr, Dk, LMAXALL, NE, TOLERENCE);

	return;
}

double Bessel_Basis::Polynomial_Interpolation2
	(const int &l, const int &ie, const double &gnorm)const
{
	const double position =  gnorm / this->Dk;
	const int iq = static_cast<int>(position);
//	if(iq<kmesh-4)
//	{
//		cout << "\n iq=" << iq << " kmesh=" << kmesh;
//		QUIT();
//	}
	assert(iq < kmesh-4);
	const double x0 = position - static_cast<double>(iq);
	const double x1 = 1.0 - x0;
	const double x2 = 2.0 - x0;
	const double x3 = 3.0 - x0;
	const double y=
		this->TableOne(l, ie, iq) * x1 * x2 * x3 / 6.0 +
        this->TableOne(l, ie, iq) * x0 * x2 * x3 / 2.0 -
        this->TableOne(l, ie, iq) * x1 * x0 * x3 / 2.0 +
        this->TableOne(l, ie, iq) * x1 * x2 * x0 / 6.0 ;
	return y;
}

// be called in Bessel_Basis::init()
void Bessel_Basis::init_TableOne(
	const bool smooth_in, // mohan add 2009-08-28
	const double &sigma_in, // mohan add 2009-08-28
	const double &ecutwfc,
	const double &rcut,
	const double &dr,
	const double &dk,
	const int &lmax,
	const int &ecut_number,
	const double &tolerence)
{
	TITLE("Bessel_Basis","init_TableOne");
	timer::tick("Spillage","TableOne");
	// check
	assert(ecutwfc > 0.0);
	assert(dr > 0.0);
	assert(dk > 0.0);

	// init kmesh
	this->kmesh = static_cast<int>(sqrt(ecutwfc) / dk) +1 + 4;
	if (kmesh % 2 == 0)++kmesh;
	cout << "\n kmesh = " << kmesh;

	// init Table One
	this->TableOne.create(lmax+1, ecut_number, kmesh);

	// init rmesh
	int rmesh = static_cast<int>( rcut / dr ) + 4;
    if (rmesh % 2 == 0) ++rmesh;
	cout << "\n rmesh = " << rmesh;

	// allocate rmesh and Jlk and eigenvalue of Jlq
	double *r = new double[rmesh];
	double *rab = new double[rmesh];
	double *jle = new double[rmesh];
	double *jlk = new double[rmesh];
	double *g = new double[rmesh]; // smooth function
	double *function = new double[rmesh];
	double *en = new double[ecut_number];
	
	for(int ir=0; ir<rmesh; ir++)
	{
		r[ir] = static_cast<double>(ir) * dr;
		rab[ir] = dr;
		if(smooth_in)
		{
			g[ir] = 1.0 - std::exp(-( (r[ir]-rcut)*(r[ir]-rcut)/2.0/sigma_in/sigma_in ) );
		}
	}
	
	// init eigenvalue of Jl
	for(int l=0; l<lmax+1; l++)
	{
		ZEROS(en, ecut_number);ZEROS(jle, rmesh);ZEROS(jlk, rmesh);

		// calculate eigenvalue for l
		Mathzone::Spherical_Bessel_Roots(ecut_number, l, tolerence, en, rcut);

		// for each eigenvalue
		for (int ie=0; ie<ecut_number; ie++)
		{
			// calculate J_{l}( en[ir]*r) 
			Mathzone::Spherical_Bessel(rmesh, r, en[ie], l, jle);

			for(int ir=0; ir<rmesh; ir++)
			{
				jle[ir] = jle[ir] * r[ir] * r[ir];
			}

			// mohan add 2009-08-28
			if(smooth_in)
			{
				for(int ir=0; ir<rmesh; ir++)
				{
					jle[ir] *= g[ir];
				}
			}
			
			for(int ik=0; ik<kmesh; ik++)
			{
				// calculate J_{l}( ik*dk*r )
				Mathzone::Spherical_Bessel(rmesh, r, ik*dk, l, jlk);

				// calculate the function will be integrated
				for(int ir=0; ir<rmesh; ir++)
				{
					function[ir] = jle[ir] * jlk[ir];
				}
				
				// make table value
				Mathzone::Simpson_Integral(rmesh, function, rab, this->TableOne(l, ie, ik) );
			}
			
		}// end ie
	}// end ;
	
	delete[] en;
	delete[] jle;
	delete[] jlk;
	delete[] rab;
	delete[] g;
	delete[] r;
	delete[] function;
	timer::tick("Spillage","TableOne");
	return;
}

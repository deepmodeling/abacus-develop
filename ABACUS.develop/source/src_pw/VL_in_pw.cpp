#include "global.h"
#include "tools.h"
#include "../src_pw/myfunc.h"
#include "VL_in_pw.h"
#include "../src_global/math_integral.h"

pseudopot_cell_vl::pseudopot_cell_vl()
{
	numeric = new bool[1];
	zp = new double[1]; 
}

pseudopot_cell_vl::~pseudopot_cell_vl()
{
	delete[] numeric;
	delete[] zp;
}


void pseudopot_cell_vl::init_vloc(const int &nggm, matrix &vloc_in)
{
	TITLE("pseudopot_cell_vl","init_vloc");

	// This routine computes the fourier coefficient of the local
	// potential vloc(ig,it) for each type of atom
	timer::tick("ppcell_vl","init_vloc",'C');

	double *vloc1d = new double[nggm];
	ZEROS(vloc1d, nggm);

	this->allocate();
	
	for (int it = 0; it < ucell.ntype; it++) 
	{
		const Atom* atom = &ucell.atoms[it];

		ZEROS(vloc1d, nggm);

		this->zp[it] = atom->zv;

		// compute V_loc(G) for a given type of atom
		if(numeric[it]==true)
		{
			this->vloc_of_g(
					atom->msh, // after cutoff 
					atom->rab,
		          	atom->r, 
					atom->vloc_at, // local potential in real space radial form.  
		          	this->zp[it],
					vloc1d);
		}
		else
		{
			WARNING_QUIT("init_vloc","not available now.");
		}

		dcopy(vloc1d, vloc_in, it);
	} 


	delete[] vloc1d;

	this->print_vloc();

	timer::tick("ppcell_vl","init_vloc",'C');
	return;
}


void pseudopot_cell_vl::allocate(void)
{
	if(test_pp>0) TITLE("pseudopot_cell_vl","allocate");
	this->vloc.create(ucell.ntype, pw.nggm);

	delete[] numeric;
	this->numeric = new bool[ucell.ntype];
	ZEROS(numeric, ucell.ntype);

	for (int it = 0; it < ucell.ntype; it++)
	{ this->numeric[it] = true; }

	// mohan change global variable 'npsx' to local variable,
	// npsx( max number of different PPs)
	// 2021-02-22
	int npsx = 50;
	delete[] zp; 
	this->zp = new double[npsx];
	ZEROS(zp, npsx);

	return;
}

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  This code transfer local pseudopotential in real space
!      radial logarithmic mesh to G space.
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


// Here we always have numeric form, i.e. numeric=ture
void pseudopot_cell_vl::vloc_of_g(
		const int& msh, 
		const double *rab, 
		const double *r, 
		const double *vloc_at, 
		const double &zp_in, 
		double *vloc_1d) const
{
	//----------------------------------------------------------------
	//    This routine computes the Fourier transform of the local
	//    part of the pseudopotential. Two types of local potentials
	//    are allowed:

	//    a) The pseudopotential is in analytic form and its fourier
	//       transform is computed analytically
	//    b) The pseudopotential is in numeric form and its fourier
	//       transform is computed numerically
	//
	//    The local pseudopotential of the US case is always in
	//    numerical form, expressed in Ry units.
	// ----------------------------------------------------------------
	int igl0;// start from |G|=0 or not.
	int ig=0;// counters on g-vectors;
	int ir=0;// counter on mesh points

	// Pseudopotentials in numerical form (Vloc_at) contain the local part
	// in order to perform the Fourier transform, a term erf(r)/r is
	// subtracted in real space and added again in G space

	double *aux = new double[msh];
	double *aux1 = new double[msh];
	ZEROS(aux, msh);
	ZEROS(aux1, msh);

	// for tests
	/*
	for(ir=0; ir<msh; ir++)
	{
		aux[ir] = r[ir] * zp_in * e2 / ucell.omega;
	}
	Integral::Simpson_Integral(msh, aux, rab, vloc_1d[0] );
	vloc_1d[0] *= 4*3.1415926;
	cout << "  vloc_1d[0]=" <<  vloc_1d[0]/pw.ngmc << endl;
	cout << "  vloc_1d[0]=" <<  vloc_1d[0]/pw.ncxyz << endl;
	*/

	// (1)
	if(pw.ggs[0] < 1.0e-8)
	{
		// first the g=0 term
		for (ir=0; ir<msh; ir++) 
		{
			// This is the |G| = 0 component of the local
			// potential giving rise to the so-called
			// "alpha*Z" term in the energy.
			aux[ir] = r [ir] * (r [ir] * vloc_at [ir] + zp_in * e2);
			//aux[ir] = r [ir] * (r [ir] * vloc_at [ir] );
		}
		Integral::Simpson_Integral(msh, aux, rab, vloc_1d[0] );
		igl0 = 1;	
	}
	else
	{
		igl0 = 0;
	}

	// (2) here the |G|>0 terms, we first compute the part of the integrand func
	// indipendent of |G| in real space
	double fac = zp_in * e2;
	for (ir = 0;ir < msh;ir++)  
	{
		aux1 [ir] = r[ir] * vloc_at [ir] + fac * erf(r[ir]);
	} 

	// here we perform the integral, after multiplying for the |G|
	// dependent part
	for (ig = igl0;ig < pw.nggm;ig++) 
	{
		double gx2= pw.ggs [ig] * ucell.tpiba2;
		double gx = std::sqrt(gx2);
		for (ir = 0;ir < msh;ir++) 
		{
			aux [ir] = aux1 [ir] * sin(gx * r [ir]) / gx;
		}
		Integral::Simpson_Integral(msh, aux, rab, vloc_1d[ig] );
		//  here we add the analytic fourier transform of the erf function
		vloc_1d[ig] -= fac * exp(- gx2 * 0.25)/ gx2;
	} // enddo

	const double d_fpi_omega = FOUR_PI/ucell.omega;//mohan add 2008-06-04
	for (ig = 0;ig < pw.nggm; ig++)
	{
		vloc_1d[ig] *= d_fpi_omega;
	}

	delete [] aux;
	delete [] aux1;
	return;
} // end subroutine vloc_of_g


void pseudopot_cell_vl::print_vloc(void)const
{
	if(MY_RANK!=0) return; //mohan fix bug 2011-10-13
	bool check_vl = true;
	if(check_vl)
	{
		for(int it=0; it<ucell.ntype; it++)
		{
			stringstream ss ;
			ss << global_out_dir << ucell.atoms[it].label << "/v_loc_g.dat" ;
			ofstream ofs_vg( ss.str().c_str() );
			for(int ig=0;ig<pw.nggm;ig++)
			{
				ofs_vg << setw(15) << pw.ggs [ig] * ucell.tpiba2 
				   	<< setw(15) << this->vloc(it, ig) << endl;
			}
			ofs_vg.close();
		}
	}
	return;
}

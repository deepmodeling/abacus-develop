#include "global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "myfunc.h"
#include "VL_in_pw.h"
#include "../module_base/math_integral.h"
#include "../module_base/timer.h"

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


void pseudopot_cell_vl::init_vloc(const int &nggm, ModuleBase::matrix &vloc_in)
{
	ModuleBase::TITLE("pseudopot_cell_vl","init_vloc");

	// This routine computes the fourier coefficient of the local
	// potential vloc(ig,it) for each type of atom
	ModuleBase::timer::tick("ppcell_vl","init_vloc");

	double *vloc1d = new double[nggm];
	ModuleBase::GlobalFunc::ZEROS(vloc1d, nggm);

	this->allocate();
	
	for (int it = 0; it < GlobalC::ucell.ntype; it++) 
	{
		const Atom* atom = &GlobalC::ucell.atoms[it];

		ModuleBase::GlobalFunc::ZEROS(vloc1d, nggm);

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
			ModuleBase::WARNING_QUIT("init_vloc","not available now.");
		}

		dcopy(vloc1d, vloc_in, it);
	} 


	delete[] vloc1d;

	this->print_vloc();

	ModuleBase::timer::tick("ppcell_vl","init_vloc");
	return;
}


void pseudopot_cell_vl::allocate(void)
{
	if(GlobalV::test_pp>0) ModuleBase::TITLE("pseudopot_cell_vl","allocate");
	this->vloc.create(GlobalC::ucell.ntype, GlobalC::pw.nggm);

	delete[] numeric;
	this->numeric = new bool[GlobalC::ucell.ntype];
	ModuleBase::GlobalFunc::ZEROS(numeric, GlobalC::ucell.ntype);

	for (int it = 0; it < GlobalC::ucell.ntype; it++)
	{ 
		this->numeric[it] = true; 
	}

	// mohan change global variable 'npsx' to local variable,
	// npsx( max number of different PPs)
	// 2021-02-22
	int npsx = 50;
	delete[] zp; 
	this->zp = new double[npsx];
	ModuleBase::GlobalFunc::ZEROS(zp, npsx);

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
	int igl0=0;// start from |G|=0 or not.
	int ig=0;// counters on g-vectors;
	int ir=0;// counter on mesh points

	// Pseudopotentials in numerical form (Vloc_at) contain the local part
	// in order to perform the Fourier transform, a term erf(r)/r is
	// subtracted in real space and added again in G space

	double *aux = new double[msh];
	double *aux1 = new double[msh];
	ModuleBase::GlobalFunc::ZEROS(aux, msh);
	ModuleBase::GlobalFunc::ZEROS(aux1, msh);

	// for tests
	/*
	for(ir=0; ir<msh; ir++)
	{
		aux[ir] = r[ir] * zp_in * e2 / GlobalC::ucell.omega;
	}
	ModuleBase::Integral::Simpson_Integral(msh, aux, rab, vloc_1d[0] );
	vloc_1d[0] *= 4*3.1415926;
	std::cout << "  vloc_1d[0]=" <<  vloc_1d[0]/GlobalC::pw.ngmc << std::endl;
	std::cout << "  vloc_1d[0]=" <<  vloc_1d[0]/GlobalC::pw.ncxyz << std::endl;
	*/

	// (1)
	if(GlobalC::pw.ggs[0] < 1.0e-8)
	{
		// first the g=0 term
		for (ir=0; ir<msh; ir++) 
		{
			// This is the |G| = 0 component of the local
			// potential giving rise to the so-called
			// "alpha*Z" term in the energy.
			aux[ir] = r [ir] * (r [ir] * vloc_at [ir] + zp_in * ModuleBase::e2);
			//aux[ir] = r [ir] * (r [ir] * vloc_at [ir] );
		}
		ModuleBase::Integral::Simpson_Integral(msh, aux, rab, vloc_1d[0] );
		igl0 = 1;	
	}
	else
	{
		igl0 = 0;
	}

	// (2) here the |G|>0 terms, we first compute the part of the integrand func
	// indipendent of |G| in real space
	double fac = zp_in * ModuleBase::e2;
	for (ir = 0;ir < msh;ir++)  
	{
		aux1 [ir] = r[ir] * vloc_at [ir] + fac * erf(r[ir]);
	} 

	// here we perform the integral, after multiplying for the |G|
	// dependent part
	for (ig = igl0;ig < GlobalC::pw.nggm;ig++) 
	{
		double gx2= GlobalC::pw.ggs [ig] * GlobalC::ucell.tpiba2;
		double gx = std::sqrt(gx2);
		for (ir = 0;ir < msh;ir++) 
		{
			aux [ir] = aux1 [ir] * sin(gx * r [ir]) / gx;
		}
		ModuleBase::Integral::Simpson_Integral(msh, aux, rab, vloc_1d[ig] );
		//  here we add the analytic fourier transform of the erf function
		vloc_1d[ig] -= fac * exp(- gx2 * 0.25)/ gx2;
	} // enddo

	const double d_fpi_omega = ModuleBase::FOUR_PI/GlobalC::ucell.omega;//mohan add 2008-06-04
	for (ig = 0;ig < GlobalC::pw.nggm; ig++)
	{
		vloc_1d[ig] *= d_fpi_omega;
	}

	delete [] aux;
	delete [] aux1;
	return;
} // end subroutine vloc_of_g


void pseudopot_cell_vl::print_vloc(void)const
{
	if(GlobalV::MY_RANK!=0) return; //mohan fix bug 2011-10-13
	bool check_vl = GlobalV::out_element_info;
	if(check_vl)
	{
		for(int it=0; it<GlobalC::ucell.ntype; it++)
		{
			std::stringstream ss ;
			ss << GlobalV::global_out_dir << GlobalC::ucell.atoms[it].label << "/v_loc_g.dat" ;
			std::ofstream ofs_vg( ss.str().c_str() );
			for(int ig=0;ig<GlobalC::pw.nggm;ig++)
			{
				ofs_vg << std::setw(15) << GlobalC::pw.ggs [ig] * GlobalC::ucell.tpiba2 
				   	<< std::setw(15) << this->vloc(it, ig) << std::endl;
			}
			ofs_vg.close();
		}
	}
	return;
}

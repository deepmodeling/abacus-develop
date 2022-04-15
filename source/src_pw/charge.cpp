//----------------------------------------------------------
// EXPLAIN : This routine calculates rhoa as the
// superposition of atomic charges.
//
// nspina is the number of spin components to be calculated
//
// nspina = 1 the total atomic charge density is calculated
// nspina = 2 the spin up and spin down atomic charge
// densities are calculated assuming an uniform atomic
// spin-polarization equal to starting_magnetization(nt)
// nspina = 4 noncollinear case. The total density is
// calculated in the first component and the magnetization
// std::vector in the other three.
//
// NB: nspina may not be equal to nspin because in some cases
// (as in update) the total charge only could be needed,
// even in a LSDA calculation.
//----------------------------------------------------------
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/memory.h"
#include "../src_parallel/parallel_reduce.h"
#include "global.h"
#include "charge.h"
#include "magnetism.h"
#include "../src_parallel/parallel_grid.h"
#include "../module_base/math_integral.h"
#include "../module_base/math_sphbes.h"
#include <vector>
#include "../module_base/timer.h"

Charge::Charge()
{
	allocate_rho = false;
    allocate_rho_final_scf = false; //LiuXh add 20180619
}


Charge::~Charge()
{
	//if(allocate_rho) //LiuXh modify 20180619
	if(allocate_rho || allocate_rho_final_scf) //LiuXh add 20180619
	{
		for(int i=0; i<GlobalV::NSPIN; i++)
		{
			delete[] rho[i];
			delete[] rhog[i];
			delete[] rho_save[i];
			delete[] rhog_save[i];
			if(XC_Functional::get_func_type() == 3)
			{
				delete[] kin_r[i];
				delete[] kin_r_save[i];
			}
		}
		delete[] rho;
		delete[] rhog;
		delete[] rho_save;
		delete[] rhog_save;
    	delete[] rho_core;
		delete[] rhog_core;
		if(XC_Functional::get_func_type() == 3)
		{
			delete[] kin_r;
			delete[] kin_r_save;
		}
	}
}


void Charge::allocate(const int &nspin_in, const int &nrxx_in, const int &ngmc_in)
{
    ModuleBase::TITLE("Charge","allocate");

	assert(allocate_rho == false);

	//  mohan add 2021-02-20
	this->nspin = nspin_in;
	this->nrxx = nrxx_in;
	this->ngmc = ngmc_in;

    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << nspin
             << " real_point_number = " << nrxx;
    }

	// allocate memory
	rho = new double*[nspin];
	rhog = new std::complex<double>*[nspin];
	rho_save = new double*[nspin];
	rhog_save = new std::complex<double>*[nspin];
	if(XC_Functional::get_func_type() == 3)
	{
		kin_r = new double*[nspin];
		kin_r_save = new double*[nspin];
	}

	for(int is=0; is<nspin; is++)
	{
		rho[is] = new double[nrxx];
		rhog[is] = new std::complex<double>[ngmc];
		rho_save[is] = new double[nrxx];
		rhog_save[is] = new std::complex<double>[ngmc];			
		ModuleBase::GlobalFunc::ZEROS(rho[is], nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog[is], ngmc);
		ModuleBase::GlobalFunc::ZEROS(rho_save[is], nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog_save[is], ngmc);
		if(XC_Functional::get_func_type() == 3)
		{
			kin_r[is] = new double[nrxx];
			ModuleBase::GlobalFunc::ZEROS(kin_r[is], nrxx);
			kin_r_save[is] = new double[nrxx];
			ModuleBase::GlobalFunc::ZEROS(kin_r_save[is], nrxx);
		}
	}

    ModuleBase::Memory::record("Charge","rho",nspin*nrxx,"double");
    ModuleBase::Memory::record("Charge","rho_save",nspin*nrxx,"double");
    ModuleBase::Memory::record("Charge","rhog",nspin*ngmc,"double");
    ModuleBase::Memory::record("Charge","rhog_save",nspin*ngmc,"double");
    ModuleBase::Memory::record("Charge","kin_r",nspin*ngmc,"double");
    ModuleBase::Memory::record("Charge","kin_r_save",nspin*ngmc,"double");

    this->rho_core = new double[nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS( rho_core, nrxx);

	this->rhog_core = new std::complex<double>[ngmc]; // reciprocal core charge
	ModuleBase::GlobalFunc::ZEROS( rhog_core, ngmc);

    ModuleBase::Memory::record("Charge","rho_core",nrxx,"double");
    ModuleBase::Memory::record("Charge","rhog_core",ngmc,"double");

	this->allocate_rho = true;
    return;
}


double Charge::sum_rho(void) const
{
	ModuleBase::TITLE("Charge","sum_rho");

    double sum_rho = 0.0;
    int nspin0 = 1;
    if(nspin==2) 
	{
		nspin0 = 2;
	}
	
	for(int is=0; is<nspin0; is++)
	{
		for(int ir=0; ir<nrxx; ir++) 
		{
			sum_rho += this->rho[is][ir];
		}
	}

	// multiply the sum of charge density by a factor
    sum_rho *= GlobalC::ucell.omega / static_cast<double>( GlobalC::pw.ncxyz );
    Parallel_Reduce::reduce_double_pool( sum_rho );

	// mohan fixed bug 2010-01-18, 
	// sum_rho may be smaller than 1, like Na bcc.
    if (sum_rho <= 0.1)
    {
		GlobalV::ofs_warning << " sum_rho=" << sum_rho << std::endl;
        ModuleBase::WARNING_QUIT("Charge::renormalize_rho","Can't find even an electron!");
    }

    return sum_rho;
}


void Charge::renormalize_rho(void)
{
    ModuleBase::TITLE("Charge","renormalize_rho");

    const double sr = this->sum_rho();
	GlobalV::ofs_warning << std::setprecision(15);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge before normalized",sr);
    const double normalize_factor = nelec / sr;

	for(int is=0; is<nspin; is++)
	{
    	for(int ir=0; ir<nrxx; ir++)
		{
			rho[is][ir] *= normalize_factor;
		}
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge after normalized",this->sum_rho());

	GlobalV::ofs_running << std::setprecision(6);
    return;
}

//-------------------------------------------------------
// superposition of atomic charges contained in the array
// rho_at (read from pseudopotential files)
// allocate work space (psic must already be allocated)
//-------------------------------------------------------
void Charge::atomic_rho(const int spin_number_need, double** rho_in)const		// Peize Lin refactor 2021.04.08
{
    ModuleBase::TITLE("Charge","atomic_rho");
    ModuleBase::timer::tick("Charge","atomic_rho");

	const ModuleBase::ComplexMatrix rho_g3d = [&]()->ModuleBase::ComplexMatrix
	{
		// use interpolation to get three dimension charge density.
		ModuleBase::ComplexMatrix rho_g3d( spin_number_need, GlobalC::pw.ngmc);
		
		// check the start magnetization
		const int startmag_type = [&]()->int
		{
			for(int it=0; it<GlobalC::ucell.ntype; it++)
			{
				if( GlobalC::ucell.magnet.start_magnetization[it] != 0.0) return 1;
			}
			return 2;
		}();
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"startmag_type",startmag_type);

		for (int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			const Atom* const atom = &GlobalC::ucell.atoms[it];

			if(!atom->flag_empty_element)		// Peize Lin add for bsse 2021.04.07
			{		
				const std::vector<double> rho_lgl = [&]()->std::vector<double>
				{
					// one dimension of charge in G space.
					std::vector<double> rho_lgl(GlobalC::pw.nggm,0);

					// mesh point of this element.
					const int mesh = atom->msh;

					//----------------------------------------------------------
					// Here we check the electron number 
					//----------------------------------------------------------
					const std::vector<double> rhoatm = [&]()->std::vector<double>
					{
						std::vector<double> rhoatm(mesh);		
						for(int ir=0; ir<mesh; ++ir)
						{
							double r2=atom->r[ir]*atom->r[ir];
							rhoatm[ir]=atom->rho_at[ir]/ModuleBase::FOUR_PI/r2;
						}
						rhoatm[0] = pow( (rhoatm[2]/rhoatm[1]), 1./(atom->r[2]-atom->r[1]) );//zws add
						rhoatm[0] = pow(rhoatm[0], atom->r[1]);
						rhoatm[0] = rhoatm[1] / rhoatm[0];

						double charge = 0.0;
						ModuleBase::Integral::Simpson_Integral(atom->msh,atom->rho_at,atom->rab,charge);
						ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge from rho_at",charge);
						assert(charge!=0.0 || charge==atom->zv);		// Peize Lin add charge==atom->zv for bsse 2021.04.07

						double scale=1.0;
						if(charge!=atom->zv)
						{
							ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"charge should be",atom->zv);
							scale = atom->zv/charge;
						}

						for(int ir=0; ir<mesh; ++ir)
						{
							rhoatm[ir] *= scale;
							rhoatm[ir] *= (ModuleBase::FOUR_PI*atom->r[ir]*atom->r[ir]);
						}
						return rhoatm;
					}();

					assert(GlobalC::ucell.meshx>0);
					std::vector<double> rho1d(GlobalC::ucell.meshx);
					//----------------------------------------------------------
					// Here we compute the G=0 term
					//----------------------------------------------------------
					if (GlobalC::pw.gstart == 1)
					{
						for (int ir = 0;ir < mesh;ir++)
						{
			//              rho1d [ir] = atom->rho_at[ir];
							rho1d[ir] = rhoatm[ir];
						}
						ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->rab, rho_lgl[0]);
					}
					if (GlobalV::test_charge>0) std::cout<<"\n |G|=0 term done." <<std::endl;
					//----------------------------------------------------------
					// Here we compute the G<>0 term
					// But if in parallel case
					// G=0 term only belong to 1 cpu.
					// Other processors start from '0'
					//----------------------------------------------------------
					for (int ig = GlobalC::pw.gstart; ig < GlobalC::pw.nggm ;ig++)
					{
						const double gx = sqrt(GlobalC::pw.ggs [ig]) * GlobalC::ucell.tpiba;
						for (int ir = 0; ir < mesh;ir++)
						{
							if ( atom->r[ir] < 1.0e-8 )
							{
								rho1d[ir] = rhoatm[ir];
								//rho1d[ir] = atom->rho_at[ir];
							}
							else
							{
								const double gxx = gx * atom->r[ir];
								rho1d[ir] = rhoatm[ir] * sin(gxx) / gxx;
								rho1d[ir] = rhoatm[ir] * sin(gxx) / gxx;
							}
						}
						ModuleBase::Integral::Simpson_Integral(mesh, rho1d.data(), atom->rab, rho_lgl[ig]);
					}
					
					if (GlobalV::test_charge>0) std::cout<<" |G|>0 term done." <<std::endl;
					//----------------------------------------------------------
					// EXPLAIN : Complete the transfer of rho from real space to
					// reciprocal space
					//----------------------------------------------------------
					for (int ig=0; ig< GlobalC::pw.nggm ; ig++)
						rho_lgl[ig] /= GlobalC::ucell.omega;
					return rho_lgl;
				}();
				//----------------------------------------------------------
				// EXPLAIN : compute the 3D atomic charge in reciprocal space
				//----------------------------------------------------------
				if(spin_number_need==1)
				{
					for (int ig=0; ig< GlobalC::pw.ngmc ;ig++)
					{
						rho_g3d(0, ig) += GlobalC::pw.strucFac(it, ig) * rho_lgl[ GlobalC::pw.ig2ngg[ig] ];
					}
				}
				// mohan add 2011-06-14, initialize the charge density according to each atom 
				else if(spin_number_need==2)
				{
					if(startmag_type==1)
					{
						for (int ig = 0; ig < GlobalC::pw.ngmc ; ig++)
						{
							const std::complex<double> swap = GlobalC::pw.strucFac(it, ig)* rho_lgl[GlobalC::pw.ig2ngg[ig]];
							//rho_g3d(0, ig) += swap * GlobalC::ucell.magnet.nelup_percent(it);
							//rho_g3d(1, ig) += swap * GlobalC::ucell.magnet.neldw_percent(it);
							const double up = 0.5 * ( 1 + GlobalC::ucell.magnet.start_magnetization[it] / atom->zv );
							const double dw = 0.5 * ( 1 - GlobalC::ucell.magnet.start_magnetization[it] / atom->zv );
							rho_g3d(0, ig) += swap * up;
							rho_g3d(1, ig) += swap * dw;
						}
					}
					// mohan add 2011-06-14
					else if(startmag_type==2)
					{
						std::complex<double> swap = ModuleBase::ZERO;
						std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
						for (int ia = 0; ia < atom->na; ia++)
						{
							//const double up = 0.5 * ( 1 + atom->mag[ia] );
							//const double dw = 0.5 * ( 1 - atom->mag[ia] );
							const double up = 0.5 * ( 1 + atom->mag[ia] / atom->zv );
							const double dw = 0.5 * ( 1 - atom->mag[ia] / atom->zv );
							//std::cout << " atom " << ia << " up=" << up << " dw=" << dw << std::endl;

							for (int ig = 0; ig < GlobalC::pw.ngmc ; ig++)
							{
								const double Gtau =
									GlobalC::pw.get_G_cartesian_projection(ig, 0) * atom->tau[ia].x + 
									GlobalC::pw.get_G_cartesian_projection(ig, 1) * atom->tau[ia].y + 
									GlobalC::pw.get_G_cartesian_projection(ig, 2) * atom->tau[ia].z;

								swap = exp(ci_tpi * Gtau) * rho_lgl[GlobalC::pw.ig2ngg[ig]];

								rho_g3d(0, ig) += swap * up;
								rho_g3d(1, ig) += swap * dw;
							}
						}
					}
				}
				else if(spin_number_need==4)
				{
					//noncolinear case
					if(startmag_type == 1)
					{
						for (int ig = 0; ig < GlobalC::pw.ngmc ; ig++)
						{
							const std::complex<double> swap = GlobalC::pw.strucFac(it, ig)* rho_lgl[GlobalC::pw.ig2ngg[ig]];
							rho_g3d(0, ig) += swap ;
							if(GlobalV::DOMAG)
							{
								rho_g3d(1, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->zv) 
								* sin(atom->angle1[0]) * cos(atom->angle2[0]);
								rho_g3d(2, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->zv) 
								* sin(atom->angle1[0]) * sin(atom->angle2[0]);
								rho_g3d(3, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->zv) 
								* cos(atom->angle1[0]);
							}
							else if(GlobalV::DOMAG_Z)
							{
								//rho_g3d(3, ig) += swap * GlobalC::ucell.magnet.start_magnetization[it];
								rho_g3d(3, ig) += swap * (GlobalC::ucell.magnet.start_magnetization[it] / atom->zv);
							}
						}
					}
					else if(startmag_type == 2)
					{//zdy-warning-not-available
						std::complex<double> swap = ModuleBase::ZERO;
						std::complex<double> ci_tpi = ModuleBase::NEG_IMAG_UNIT * ModuleBase::TWO_PI;
						for(int ia = 0;ia<atom->na;ia++)
						{
							for (int ig = 0; ig < GlobalC::pw.ngmc ; ig++)
							{
								const double Gtau =
									GlobalC::pw.get_G_cartesian_projection(ig, 0) * atom->tau[ia].x + 
									GlobalC::pw.get_G_cartesian_projection(ig, 1) * atom->tau[ia].y + 
									GlobalC::pw.get_G_cartesian_projection(ig, 2) * atom->tau[ia].z;

								swap = exp(ci_tpi * Gtau) * rho_lgl[GlobalC::pw.ig2ngg[ig]];

								rho_g3d(0, ig) += swap;
								if(GlobalV::DOMAG)
								{
									rho_g3d(1, ig) += swap * (atom->mag[ia] / atom->zv) 
										* sin(atom->angle1[ia]) * cos(atom->angle2[ia]);
									rho_g3d(2, ig) += swap * (atom->mag[ia] / atom->zv) 
										* sin(atom->angle1[ia]) * sin(atom->angle2[ia]);
									rho_g3d(3, ig) += swap * (atom->mag[ia] / atom->zv) 
										* cos(atom->angle1[ia]);
								}
								else if(GlobalV::DOMAG_Z)
								{
									rho_g3d(3, ig) += swap * (atom->mag[ia] / atom->zv);
								}
							}
						}
					}
				}
				else
				{
					ModuleBase::WARNING_QUIT("Charge::spin_number_need"," Either 1 or 2 or 4, check SPIN number !");
				}
			}
		}
		return rho_g3d;
	}();

	assert( spin_number_need > 0 );
	std::vector<double> ne(spin_number_need);
    for (int is = 0; is < spin_number_need;is++)
    {
        GlobalC::UFFT.ToRealSpace( is, rho_g3d, rho_in[is]);

		for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
			ne[is] += rho_in[is][ir];
		ne[is] *= GlobalC::ucell.omega/(double)GlobalC::pw.ncxyz; 
		Parallel_Reduce::reduce_double_pool( ne[is] );

        // we check that everything is correct
        double neg = 0.0;
        double rea = 0.0;
        double ima = 0.0;
		double sumrea = 0.0;
        for (int ir=0;ir < GlobalC::pw.nrxx; ir++)
        {
            rea = GlobalC::UFFT.porter[ir].real();
			sumrea += rea;
            neg += std::min(0.0, rea);
            ima += abs(GlobalC::UFFT.porter[ir].imag());
        }

		Parallel_Reduce::reduce_double_pool( neg );	
		Parallel_Reduce::reduce_double_pool( ima );	
		Parallel_Reduce::reduce_double_pool( sumrea );	

		// mohan fix bug 2011-04-03
        neg = neg / (double)GlobalC::pw.ncxyz * GlobalC::ucell.omega;
        ima = ima / (double)GlobalC::pw.ncxyz * GlobalC::ucell.omega;
		sumrea = sumrea / (double)GlobalC::pw.ncxyz * GlobalC::ucell.omega;

        if( ((neg<-1.0e-4) && (is==0||GlobalV::NSPIN==2)) || ima>1.0e-4)
        {
            GlobalV::ofs_warning << " Warning: negative or imaginary starting charge : " ;
            GlobalV::ofs_warning << " neg = " << neg
                 << " ima = " << ima
                 << " SPIN = " << is << std::endl;
        }

//		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;
//		std::cout << " sum rho for spin " << is << " = " << sumrea << std::endl;

    }//end is

//	for(int it=0; it<GlobalC::ucell.ntype; it++)
//	{
		//std::cout << " nelup_percent = " << GlobalC::ucell.magnet.nelup_percent(it) << std::endl;
		//std::cout << " neldw_percent = " << GlobalC::ucell.magnet.neldw_percent(it) << std::endl;
//	}


	double ne_tot = 0.0;
	int spin0=1;
	if(spin_number_need == 2) spin0 = spin_number_need;
	for(int is=0; is<spin0; ++is)
	{
		GlobalV::ofs_warning << "\n SETUP ATOMIC RHO FOR SPIN " << is+1 << std::endl;
		ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"Electron number from rho",ne[is]);
		ne_tot += ne[is];
	}
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"total electron number from rho",ne_tot);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_warning,"should be",nelec);
	for(int is=0; is<spin_number_need; ++is)
		for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
			rho_in[is][ir] = rho_in[is][ir] / ne_tot * nelec;

	//wenfei 2021-7-29 : initial tau = 3/5 rho^2/3, Thomas-Fermi
	if(XC_Functional::get_func_type() == 3)
	{
		const double pi = 3.141592653589790;
		double fact = (3.0/5.0)*pow(3.0*pi*pi,2.0/3.0);
		int nspin = spin_number_need;
		//ofstream test_tau0("tau0");
		for(int is=0; is<spin_number_need; ++is)
			for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
			{
				kin_r[is][ir] = fact * pow(abs(rho_in[is][ir])*nspin,5.0/3.0)/nspin;
				//test_tau0 << rho_in[is][ir] << " " << kin_r[is][ir] << endl;
			}
	}

	// if TWO_EFEMI, 
	// the total magnetism will affect the calculation of
	// occupations.
	// GlobalC::ucell.magnet.compute_magnetization();

	//GlobalV::ofs_running << " Superposition of atomic wave function as First-Charge done." << std::endl;
	//2014-06-22

    ModuleBase::timer::tick("Charge","atomic_rho");
    return;
}


//==========================================================
// computes the core charge on the real space 3D mesh.
//==========================================================
void Charge::set_rho_core(
    const ModuleBase::ComplexMatrix &structure_factor
)
{
    ModuleBase::TITLE("Charge","set_rho_core");
    ModuleBase::timer::tick("Charge","set_rho_core");

    //double eps = 1.e-10;
    GlobalC::en.etxcc = 0.0;
//----------------------------------------------------------
// LOCAL VARIABLES :
// counter on mesh points
// counter on atomic types
// counter on g vectors
//----------------------------------------------------------
    //int ir = 0;
    //int it = 0;
    //int ig = 0;

    bool bl = false;
    for (int it = 0; it<GlobalC::ucell.ntype; it++)
    {
        if (GlobalC::ucell.atoms[it].nlcc)
        {
            bl = true;
            break;
        }
    }

    if (!bl)
    {
        ModuleBase::GlobalFunc::ZEROS( this->rho_core, GlobalC::pw.nrxx);
    	ModuleBase::timer::tick("Charge","set_rho_core");
        return;
    }

    double *rhocg = new double[GlobalC::pw.nggm];
    ModuleBase::GlobalFunc::ZEROS(rhocg, GlobalC::pw.nggm );

	// three dimension.
    std::complex<double> *vg = new std::complex<double>[GlobalC::pw.ngmc];	

    for (int it = 0; it < GlobalC::ucell.ntype;it++)
    {
        if (GlobalC::ucell.atoms[it].nlcc)
        {
//----------------------------------------------------------
// EXPLAIN : drhoc compute the radial fourier transform for
// each shell of g vec
//----------------------------------------------------------
            this->non_linear_core_correction(
                GlobalC::ppcell.numeric,
                GlobalC::ucell.atoms[it].msh,
                GlobalC::ucell.atoms[it].r,
                GlobalC::ucell.atoms[it].rab,
                GlobalC::ucell.atoms[it].rho_atc,
                rhocg);
//----------------------------------------------------------
// EXPLAIN : multiply by the structure factor and sum
//----------------------------------------------------------
            for (int ig = 0; ig < GlobalC::pw.ngmc ; ig++)
            {
                vg[ig] += structure_factor(it, ig) * rhocg[GlobalC::pw.ig2ngg[ig]];
            }
        }
    }

	// for tmp use.
	for(int ig=0; ig< GlobalC::pw.ngmc; ig++)
	{
		this->rhog_core[ig] = vg[ig];
	}

    GlobalC::UFFT.ToRealSpace(vg, this->rho_core);

    // test on the charge and computation of the core energy
    double rhoima = 0.0;
    double rhoneg = 0.0;
    for (int ir = 0; ir < GlobalC::pw.nrxx; ir++)
    {
        rhoneg += min(0.0, GlobalC::UFFT.porter[ir].real());
        rhoima += abs(GlobalC::UFFT.porter[ir].imag());
        // NOTE: Core charge is computed in reciprocal space and brought to real
        // space by FFT. For non smooth core charges (or insufficient cut-off)
        // this may result in negative values in some grid points.
        // Up to October 1999 the core charge was forced to be positive definite.
        // This induces an error in the force, and probably stress, calculation if
        // the number of grid points where the core charge would be otherwise neg
        // is large. The error disappears for sufficiently high cut-off, but may be
        // rather large and it is better to leave the core charge as it is.
        // If you insist to have it positive definite (with the possible problems
        // mentioned above) uncomment the following lines.  SdG, Oct 15 1999
    }

	// mohan fix bug 2011-04-03
	Parallel_Reduce::reduce_double_pool( rhoneg );
	Parallel_Reduce::reduce_double_pool( rhoima );

	// mohan changed 2010-2-2, make this same as in atomic_rho.
	// still lack something......
    rhoneg /= GlobalC::pw.ncxyz * GlobalC::ucell.omega;
    rhoima /= GlobalC::pw.ncxyz * GlobalC::ucell.omega;

    // calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
    // The term was present in previous versions of the code but it shouldn't
    delete [] rhocg;
    delete [] vg;
    ModuleBase::timer::tick("Charge","set_rho_core");
    return;
} // end subroutine set_rhoc



void Charge::non_linear_core_correction
(
    const bool &numeric,
    const int mesh,
    const double *r,
    const double *rab,
    const double *rhoc,
    double *rhocg)
{
    ModuleBase::TITLE("charge","drhoc");
    double gx = 0.0;
    double rhocg1 = 0.0;
    double *aux;

    // here we compute the fourier transform is the charge in numeric form
    if (numeric)
    {
        aux = new double [mesh];
        // G=0 term

        int igl0 = 0;
        if (GlobalC::pw.ggs [0] < 1.0e-8)
        {
            for (int ir = 0;ir < mesh; ir++)
            {
                aux [ir] = r [ir] * r [ir] * rhoc [ir];
            }
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            //rhocg [1] = fpi * rhocg1 / omega;
            rhocg [0] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;//mohan modify 2008-01-19
            igl0 = 1;
        }

        // G <> 0 term
        for (int igl = igl0; igl < GlobalC::pw.nggm;igl++) 
        {
            gx = sqrt(GlobalC::pw.ggs [igl] * GlobalC::ucell.tpiba2);
            ModuleBase::Sphbes::Spherical_Bessel(mesh, r, gx, 0, aux);
            for (int ir = 0;ir < mesh; ir++) 
            {
                aux [ir] = r[ir] * r[ir] * rhoc [ir] * aux [ir];
            } //  enddo
            ModuleBase::Integral::Simpson_Integral(mesh, aux, rab, rhocg1);
            rhocg [igl] = ModuleBase::FOUR_PI * rhocg1 / GlobalC::ucell.omega;
        } //  enddo
        delete [] aux;
    }
    else
    {
        // here the case where the charge is in analytic form,
        // check old version before 2008-12-9
    }
    return;
}


//----------------------------------------------------------
// NAME : sum_band
//----------------------------------------------------------
void Charge::sum_band(void)
{
    ModuleBase::TITLE("Charge","sum_band");
    ModuleBase::timer::tick("Charge","sum_band");
//----------------------------------------------------------
// Calculates the symmetrized charge density and sum of
// occupied eigenvalues.
//----------------------------------------------------------

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		ModuleBase::GlobalFunc::ZEROS(rho[is], GlobalC::pw.nrxx);
		if (XC_Functional::get_func_type() == 3)
		{
			ModuleBase::GlobalFunc::ZEROS(kin_r[is], GlobalC::pw.nrxx);	
		}
	}
	
    sum_band_k();

    // Symmetrization of the charge density (and local magnetization)
    ModuleBase::timer::tick("Charge","sum_band");
    return;
}

void Charge::sum_band_k(void)
{
	ModuleBase::TITLE("Charge","sum_band_k");
	GlobalC::en.eband = 0.0;

	std::complex<double>* porter = GlobalC::UFFT.porter;
	std::complex<double>* porter1 = nullptr;
	if(GlobalV::NSPIN==4) porter1 = new std::complex<double>[GlobalC::pw.nrxx];//added by zhengdy-soc

	for (int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		//std::cout << "\n ik=" << ik;
		if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		
		//  here we compute the band energy: the sum of the eigenvalues
		if(GlobalV::NSPIN==4)
		{
			for (int ibnd = 0;ibnd < GlobalV::NBANDS;ibnd++)
			{
				///
				///only occupied band should be calculated.
				///
				if(GlobalC::wf.wg(ik, ibnd)<ModuleBase::threshold_wg) continue;
				GlobalC::en.eband += GlobalC::wf.ekb[ik][ibnd] * GlobalC::wf.wg(ik, ibnd);
				ModuleBase::GlobalFunc::ZEROS( porter, GlobalC::pw.nrxx );
				for (int ig = 0;ig < GlobalC::kv.ngk[ik] ; ig++)
 				{
					porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik, ig)] ] = GlobalC::wf.evc[ik](ibnd, ig);
				}
				GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter, 1);
				if(GlobalV::NPOL ==2)
				{
					ModuleBase::GlobalFunc::ZEROS( porter1, GlobalC::pw.nrxx );
					for (int ig = 0;ig < GlobalC::kv.ngk[ik] ; ig++)
					{
						porter1[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik, ig)] ] = GlobalC::wf.evc[ik](ibnd, ig + GlobalC::wf.npwx);
					}
					GlobalC::pw.FFT_wfc.FFT3D(porter1, 1);
				}
				const double w1 = GlobalC::wf.wg(ik, ibnd) / GlobalC::ucell.omega;

				// Increment the charge density in chr.rho for real space
				if (w1 != 0.0)
				{
					for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
					{
						rho[0][ir]+=w1* (norm( porter[ir])+ norm(porter1[ir]));
					}
				}
				// In this case, calculate the three components of the magnetization
				if(GlobalV::DOMAG){
					if(w1 != 0.0)
						for(int ir= 0;ir<GlobalC::pw.nrxx;ir++)
						{
							rho[1][ir] += w1 * 2.0 * (porter[ir].real()* porter1[ir].real()
								+ porter[ir].imag()* porter1[ir].imag());
							rho[2][ir] += w1 * 2.0 * (porter[ir].real()* porter1[ir].imag()
								- porter1[ir].real()* porter[ir].imag());
							rho[3][ir] += w1 * (norm(porter[ir]) - norm(porter1[ir]));
						}
				}
				else if(GlobalV::DOMAG_Z){
					if(w1 != 0.0)
						for(int ir= 0;ir<GlobalC::pw.nrxx;ir++)
						{
							rho[1][ir] = 0;
							rho[2][ir] = 0;
							rho[3][ir] += w1 * (norm(porter[ir]) - norm(porter1[ir]));
						}
				}
				else for(int is= 1;is<4;is++)
					for(int ir = 0;ir<GlobalC::pw.nrxx;ir++) rho[is][ir] = 0;
			}
		}
		else
		for (int ibnd = 0;ibnd < GlobalV::NBANDS;ibnd++)
		{
			///
			///only occupied band should be calculated.
			///
			if(GlobalC::wf.wg(ik, ibnd)<ModuleBase::threshold_wg) continue;
			GlobalC::en.eband += GlobalC::wf.ekb[ik][ibnd] * GlobalC::wf.wg(ik, ibnd);
			//std::cout << "\n ekb = " << GlobalC::wf.ekb[ik][ibnd] << " wg = " << GlobalC::wf.wg(ik, ibnd);

			ModuleBase::GlobalFunc::ZEROS( porter, GlobalC::pw.nrxx );
			for (int ig = 0;ig < GlobalC::kv.ngk[ik] ; ig++)
			{
				porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik, ig)] ] = GlobalC::wf.evc[ik](ibnd, ig);
			}
			GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter, 1);

			const double w1 = GlobalC::wf.wg(ik, ibnd) / GlobalC::ucell.omega;

			if (w1 != 0.0)
			{
				for (int ir=0; ir<GlobalC::pw.nrxx; ir++) 
				{
					rho[GlobalV::CURRENT_SPIN][ir]+=w1* norm( porter[ir] );
				}
			}

			//kinetic energy density
			if (XC_Functional::get_func_type() == 3)
			{
				for (int j=0; j<3; j++)
				{
					ModuleBase::GlobalFunc::ZEROS( porter, GlobalC::pw.nrxx );
					for (int ig = 0;ig < GlobalC::kv.ngk[ik] ; ig++)
					{
						double fact = GlobalC::pw.get_GPlusK_cartesian_projection(ik,GlobalC::wf.igk(ik,ig),j) * GlobalC::ucell.tpiba;
						porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik, ig)] ] = GlobalC::wf.evc[ik](ibnd, ig) * complex<double>(0.0,fact);
					}
					GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter, 1);
					for (int ir=0; ir<GlobalC::pw.nrxx; ir++) 
					{
						kin_r[GlobalV::CURRENT_SPIN][ir]+=w1* norm( porter[ir] );
					}
				}	
			}
		}
	} // END DO k_loop
	if(GlobalV::NSPIN==4) delete[] porter1;

#ifdef __MPI
	this->rho_mpi();
	if(GlobalV::CALCULATION!="scf-sto" && GlobalV::CALCULATION!="relax-sto" && GlobalV::CALCULATION!="md-sto") //qinarui add it temporarily.
	{
    //==================================
    // Reduce all the Energy in each cpu
    //==================================
	GlobalC::en.eband /= GlobalV::NPROC_IN_POOL;
	Parallel_Reduce::reduce_double_all( GlobalC::en.eband );
	}
#endif
	// check how many electrons on this grid.
	/*
	double sum = 0.0;
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		sum += rho1[ir];
	}
	std::cout << "\n sum=" << sum * GlobalC::ucell.omega / GlobalC::pw.nrxx << std::endl;
	*/

    return;
}


#ifdef __MPI
void Charge::rho_mpi(void)
{
	ModuleBase::TITLE("Charge","rho_mpi");
    if (GlobalV::NPROC==1) return;
	if((GlobalV::CALCULATION=="scf-sto" || GlobalV::CALCULATION=="relax-sto" || GlobalV::CALCULATION=="md-sto")&&GlobalV::NPROC_IN_POOL==1) 
		return;//qinarui add it temporarily.
    ModuleBase::timer::tick("Charge","rho_mpi");
    int ir;//counters on real space mesh point.
    int iz;//counters on z direction of fft grid.
    int ip;//counters on processors

    //=========================================
    // There are two steps to do before getting
    // the final charge:
    // (1) sum up the plane rhos in each pool.
    // (2) sum up all the rhos from all pools.
    //=========================================

    //=================================================
    // Searching in all planes, and calculated each
    // plane belong to which processor.
    // Count number of planes for each cpu in this pool
	// num_z: how many planes on processor 'ip'
    //=================================================
    int *num_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(num_z, GlobalV::NPROC_IN_POOL);
    for (iz=0;iz<GlobalC::pw.nbz;iz++)
    {
        ip = iz % GlobalV::NPROC_IN_POOL;
        num_z[ip]++;
    }

	// mohan update 2011-04-26
	for(int ip=0; ip<GlobalV::NPROC_IN_POOL; ip++)
	{
		num_z[ip]*=GlobalC::pw.bz;
	}

    //=======================================
    // Find current number of planes (nz)
	// start_z: start position of z in 
	// processor ip.
    //=======================================
    int *start_z = new int[GlobalV::NPROC_IN_POOL];
    ModuleBase::GlobalFunc::ZEROS(start_z, GlobalV::NPROC_IN_POOL);
    for (ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        start_z[ip] = start_z[ip-1]+num_z[ip-1];
    }

    //====================================================
    // Find "number of data" in each processor in each pool
    //====================================================
    int *rec = new int[GlobalV::NPROC_IN_POOL];
	ModuleBase::GlobalFunc::ZEROS(rec, GlobalV::NPROC_IN_POOL);
    const int ncxy = GlobalC::pw.ncx * GlobalC::pw.ncy;
    for (ip=0;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        rec[ip] = num_z[ip]*ncxy;
    }

    //======================================================
    // Find current "index of data" in each cpu in this pool
	// also, we mean start position of data.
    //======================================================
    int *dis = new int[GlobalV::NPROC_IN_POOL];
	ModuleBase::GlobalFunc::ZEROS(dis, GlobalV::NPROC_IN_POOL);
    for (ip=1;ip<GlobalV::NPROC_IN_POOL;ip++)
    {
        dis[ip]=dis[ip-1]+rec[ip-1];
    }

    //==========================================
    // Collection of rho in each pool
    // ( according to different k distribution,
    // so the rho in each pool is different
    //==========================================
    double *rho_tmp = new double[GlobalC::pw.nrxx];
    double *rho_tot = new double[GlobalC::pw.ncxyz];
    double *rho_tot_aux = new double[GlobalC::pw.ncxyz];
	ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, GlobalC::pw.ncxyz);

	double *tau_tmp;
	double *tau_tot;
	double *tau_tot_aux;

	if(XC_Functional::get_func_type() == 3)
	{
    	tau_tmp = new double[GlobalC::pw.nrxx];
	    tau_tot = new double[GlobalC::pw.ncxyz];
    	tau_tot_aux = new double[GlobalC::pw.ncxyz];
		ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, GlobalC::pw.ncxyz);
	}

    for (int is=0; is< GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(rho_tot, GlobalC::pw.ncxyz);
		if(XC_Functional::get_func_type() == 3) ModuleBase::GlobalFunc::ZEROS(tau_tot, GlobalC::pw.ncxyz);

		for (ir=0;ir<GlobalC::pw.nrxx;ir++)
		{
			rho_tmp[ir] = this->rho[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
			if(XC_Functional::get_func_type() == 3)
			{
				tau_tmp[ir] = this->kin_r[is][ir] / static_cast<double>(GlobalV::NPROC_IN_POOL);
			}
		}

        MPI_Allgatherv(rho_tmp, GlobalC::pw.nrxx, MPI_DOUBLE, rho_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
		if(XC_Functional::get_func_type() == 3)
		{
        	MPI_Allgatherv(tau_tmp, GlobalC::pw.nrxx, MPI_DOUBLE, tau_tot, rec, dis, MPI_DOUBLE, POOL_WORLD);
		}
        //=================================================================
        // Change the order of rho_tot in each pool , make them consistent
        // this is the most complicated part !!
        //=================================================================
        ModuleBase::GlobalFunc::ZEROS(rho_tot_aux, GlobalC::pw.ncxyz);
		if(XC_Functional::get_func_type() == 3)
		{
        	ModuleBase::GlobalFunc::ZEROS(tau_tot_aux, GlobalC::pw.ncxyz);
		}

        for (ip=0;ip<GlobalV::NPROC_IN_POOL;ip++)
        {
            for (ir=0;ir<ncxy;ir++)
            {
                for (iz=0;iz<num_z[ip];iz++)
                {
					// -------------------------------------------------
					// very carefully with the order of charge density.
					// the data (ir,iz) is now in processor 'ip'.
					// different POOL has different ordering.
					// we want to collect them in each processor
					// in a unit format,
					// and then reduce among all POOLS to yield
					// the correct charge density.
					// we know the division of 'z' is indipendent
					// in each processor, so the 'unit format'
					// must have no relationship with 'z' divide method.
					// -------------------------------------------------
					// rot_tot_aux : suitable among all pools.
					// (1) the data save along z direction.
					// (2) and each element number of group 'z data' 
					// is 'GlobalC::pw.ncz'
					// (3) however, the data rearrange is occured
					// between [ start_z[ip], start_z[ip]+num_z[ip] )
					// (4) start_z[ip] + iz yields correct z coordiante.
					// -------------------------------------------------
					// rot_tot: suitable for local pool.
					// (1) the data save along z direction, only 
					// in a small distance.
					// (2) however, the number of z in each processor 
					// 'ip' is num_z[ip]
					// (3) the index of data increases with the ip,
					// so the data on large 'ip' processor must 
					// have large 'start position', which we label
					// start_z[ip] * ncxy.
					// -------------------------------------------------
                    rho_tot_aux[GlobalC::pw.ncz*ir    + start_z[ip]      + iz]
                      = rho_tot[num_z[ip]*ir + start_z[ip]*ncxy + iz];
					if(XC_Functional::get_func_type() == 3)
					{
                    	tau_tot_aux[GlobalC::pw.ncz*ir    + start_z[ip]      + iz]
                      	  = tau_tot[num_z[ip]*ir + start_z[ip]*ncxy + iz];
					}	
                }
            }
        }
        //==================================
        // Reduce all the rho in each cpu
        //==================================
		if(GlobalV::CALCULATION=="scf-sto" || GlobalV::CALCULATION=="relax-sto" || GlobalV::CALCULATION=="md-sto") //qinarui add it temporarily.
		{
			MPI_Allreduce(rho_tot_aux,rho_tot,GlobalC::pw.ncxyz,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
			if(XC_Functional::get_func_type() == 3)
			{
				MPI_Allreduce(tau_tot_aux,tau_tot,GlobalC::pw.ncxyz,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
			}
		}
		else
        MPI_Allreduce(rho_tot_aux,rho_tot,GlobalC::pw.ncxyz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		if(XC_Functional::get_func_type() == 3)
		{
   	    	MPI_Allreduce(tau_tot_aux,tau_tot,GlobalC::pw.ncxyz,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		}
        
		//=====================================
        // Change the order of rho in each cpu
        //=====================================
        for (ir =0;ir<ncxy;ir++)
        {
            for (iz=0;iz<num_z[GlobalV::RANK_IN_POOL];iz++)
            {
                this->rho[is][num_z[GlobalV::RANK_IN_POOL]*ir+iz] = rho_tot[GlobalC::pw.ncz*ir + start_z[GlobalV::RANK_IN_POOL] + iz ];
				if(XC_Functional::get_func_type() == 3)
				{
                	this->kin_r[is][num_z[GlobalV::RANK_IN_POOL]*ir+iz] = tau_tot[GlobalC::pw.ncz*ir + start_z[GlobalV::RANK_IN_POOL] + iz ];
				}	
            }
        }
    }
    delete[] rho_tot_aux;
    delete[] rho_tot;
    delete[] rho_tmp;
    
	if(XC_Functional::get_func_type() == 3)
	{
		delete[] tau_tot_aux;
    	delete[] tau_tot;
	    delete[] tau_tmp;
	}
    delete[] rec;
    delete[] dis;

    delete[] num_z;
    delete[] start_z;
    ModuleBase::timer::tick("Charge","rho_mpi");
    return;
}
#endif


void Charge::save_rho_before_sum_band(void)
{
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
    	ModuleBase::GlobalFunc::DCOPY( rho[is], rho_save[is], GlobalC::pw.nrxx);
    	if(XC_Functional::get_func_type() == 3) ModuleBase::GlobalFunc::DCOPY( kin_r[is], kin_r_save[is], GlobalC::pw.nrxx);
    }
    return;
}


double Charge::check_ne(const double *rho_in) const 
{
	double ne= 0.0;
	for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
	{
		ne += rho_in[ir];
	}
	Parallel_Reduce::reduce_double_pool( ne );
	ne = ne * GlobalC::ucell.omega / (double)GlobalC::pw.ncxyz;
	std::cout << std::setprecision(10);
	std::cout << " check the electrons number from rho, ne =" << ne << std::endl;
	std::cout << std::setprecision(6);
	return ne;
}


//LiuXh add 20180619
void Charge::init_final_scf()
{
    ModuleBase::TITLE("Charge","init_after_scf");

	assert(allocate_rho_final_scf == false);

    if (GlobalV::test_charge > 1)
    {
        std::cout << "\n spin_number = " << GlobalV::NSPIN
             << " real_point_number = " << GlobalC::pw.nrxx;
    }

	// allocate memory
	rho = new double*[GlobalV::NSPIN];
	rhog = new std::complex<double>*[GlobalV::NSPIN];
	rho_save = new double*[GlobalV::NSPIN];
	rhog_save = new std::complex<double>*[GlobalV::NSPIN];

	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		rho[is] = new double[GlobalC::pw.nrxx];
		rhog[is] = new std::complex<double>[GlobalC::pw.ngmc];
		rho_save[is] = new double[GlobalC::pw.nrxx];
		rhog_save[is] = new std::complex<double>[GlobalC::pw.ngmc];			
		ModuleBase::GlobalFunc::ZEROS(rho[is], GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog[is], GlobalC::pw.ngmc);
		ModuleBase::GlobalFunc::ZEROS(rho_save[is], GlobalC::pw.nrxx);
		ModuleBase::GlobalFunc::ZEROS(rhog_save[is], GlobalC::pw.ngmc);
	}

    ModuleBase::Memory::record("Charge","rho",GlobalV::NSPIN*GlobalC::pw.nrxx,"double");
    ModuleBase::Memory::record("Charge","rho_save",GlobalV::NSPIN*GlobalC::pw.nrxx,"double");
    ModuleBase::Memory::record("Charge","rhog",GlobalV::NSPIN*GlobalC::pw.ngmc,"double");
    ModuleBase::Memory::record("Charge","rhog_save",GlobalV::NSPIN*GlobalC::pw.ngmc,"double");

    this->rho_core = new double[GlobalC::pw.nrxx]; // core charge in real space
    ModuleBase::GlobalFunc::ZEROS( rho_core, GlobalC::pw.nrxx);

	this->rhog_core = new std::complex<double>[GlobalC::pw.ngmc]; // reciprocal core charge
	ModuleBase::GlobalFunc::ZEROS( rhog_core, GlobalC::pw.ngmc);

    ModuleBase::Memory::record("Charge","rho_core",GlobalC::pw.nrxx,"double");
    ModuleBase::Memory::record("Charge","rhog_core",GlobalC::pw.ngmc,"double");

	this->allocate_rho_final_scf = true;
    return;
}

//=========================================================
// calculate total number of electrons (nelec) and default
// number of bands (GlobalV::NBANDS).
//=========================================================
#include "occupy.h"
void Charge::cal_nelec(void)
{
	ModuleBase::TITLE("UnitCell_pseudo","cal_nelec");
	//=======================================================
	// calculate the total number of electrons in the system
	// if nelec <>0; use input number (setup.f90)
	//=======================================================

	GlobalV::ofs_running << "\n SETUP THE ELECTRONS NUMBER" << std::endl;

	if (nelec == 0)
	{
		for (int it = 0; it < GlobalC::ucell.ntype;it++)
		{
			std::stringstream ss1, ss2;
			ss1 << "electron number of element " << GlobalC::ucell.atoms[it].label;
			const int nelec_it = GlobalC::ucell.atoms[it].zv * GlobalC::ucell.atoms[it].na;
			nelec += nelec_it;
			ss2 << "total electron number of element " << GlobalC::ucell.atoms[it].label; 
			
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss1.str(),GlobalC::ucell.atoms[it].zv);
			ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,ss2.str(),nelec_it);
		}
	}

	//ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total nelec",nelec);

	//=======================================
	// calculate number of bands (setup.f90)
	//=======================================
	double occupied_bands = static_cast<double>(nelec/ModuleBase::DEGSPIN);	
	if(GlobalV::LSPINORB==1) occupied_bands = static_cast<double>(nelec);

	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0; //mohan fix 2012-04-16
	}

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"occupied bands",occupied_bands);
	
	// mohan add 2010-09-04
    //std::cout << "nbands(GlobalC::ucell) = " <<GlobalV::NBANDS <<std::endl;
	if(GlobalV::NBANDS==occupied_bands)
	{
		if( Occupy::gauss() || Occupy::tetra() )
		{
			ModuleBase::WARNING_QUIT("UnitCell_pseudo::cal_nelec","for smearing, num. of bands > num. of occupied bands");
		}
	}
	
	if ( GlobalV::CALCULATION!="scf-sto" && GlobalV::CALCULATION!="relax-sto" && GlobalV::CALCULATION!="md-sto" ) //qianrui 2021-2-20
	{
		if(GlobalV::NBANDS == 0)
		{
			if(GlobalV::NSPIN == 1)
			{
				const int nbands1 = static_cast<int>(occupied_bands) + 10;
				const int nbands2 = static_cast<int>(1.2 * occupied_bands);
				GlobalV::NBANDS = std::max(nbands1, nbands2);
				if(GlobalV::BASIS_TYPE!="pw") GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
			}
			else if (GlobalV::NSPIN ==2 || GlobalV::NSPIN == 4)
			{
				const int nbands3 = nelec + 20;
				const int nbands4 = 1.2 * nelec;
				GlobalV::NBANDS = std::max(nbands3, nbands4);
				if(GlobalV::BASIS_TYPE!="pw") GlobalV::NBANDS = std::min(GlobalV::NBANDS, GlobalV::NLOCAL);
			}
			ModuleBase::GlobalFunc::AUTO_SET("NBANDS",GlobalV::NBANDS);
		}
		//else if ( GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax") //pengfei 2014-10-13
		else
		{
			if(GlobalV::NBANDS < occupied_bands) ModuleBase::WARNING_QUIT("unitcell","Too few bands!");
			if(GlobalV::NBANDS < GlobalC::ucell.magnet.get_nelup() ) 
			{
				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nelup",GlobalC::ucell.magnet.get_nelup());
				ModuleBase::WARNING_QUIT("unitcell","Too few spin up bands!");
			}
			if(GlobalV::NBANDS < GlobalC::ucell.magnet.get_neldw() )
			{
				ModuleBase::WARNING_QUIT("unitcell","Too few spin down bands!");
			}
		}
	}

	// mohan update 2021-02-19
    // mohan add 2011-01-5
    if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
    {
        if( GlobalV::NBANDS > GlobalV::NLOCAL )
        {
            ModuleBase::WARNING_QUIT("UnitCell_pseudo::cal_nwfc","NLOCAL < NBANDS");
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NLOCAL",GlobalV::NLOCAL);
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NBANDS",GlobalV::NBANDS);
        }
    }

	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"NBANDS",GlobalV::NBANDS);
	return;
}

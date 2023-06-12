#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "bessel_basis.h"
#include "module_base/parallel_common.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"
#include "module_base/timer.h"

Bessel_Basis::Bessel_Basis()
{
	Ecut_number = 0;
	Dk = 0.0;
}

Bessel_Basis::~Bessel_Basis()
{
}


// the function is called in numerical_basis.
void Bessel_Basis::init(
	const bool start_from_file,
	const double &ecutwfc,
	const int &ntype,
	const int &lmax_in,
	const bool &smooth,
	const double &sigma,
	const double &rcut_in,
	const double &tol_in,
	const UnitCell& ucell,
	const double &dk,
	const double &dr
	)
{
	ModuleBase::TITLE("Bessel_Basis", "init");
	this->Dk = dk;
	this->ecut = ecutwfc;
	this->rcut = rcut_in;
	this->tolerence = tol_in;
    this->smooth = smooth;
    this->sigma = sigma;

    //----------------------------------------------
    // setup Ecut_number
    // ne * pi / rcut = sqrt(ecut) (Rydberg)
    //----------------------------------------------
    //	this->Ecut_number = static_cast<int>( sqrt( 2.0 * ecut )* rcut/ModuleBase::PI );// hartree
    this->Ecut_number = static_cast<int>(sqrt(ecut) * rcut / ModuleBase::PI); // Rydberg Unit.
    assert(this->Ecut_number > 0);

    //------------------
	// Making a table
	//------------------

	this->init_TableOne( smooth, sigma, ecutwfc, rcut, dr, Dk, lmax_in, Ecut_number, tolerence);

//-----------------------------------------------
// for test.
//-----------------------------------------------
//	GlobalV::ofs_running << "\n TableOne:";
//	for(int i=0; i<TableOne.getSize(); i++)
//	{
//		GlobalV::ofs_running << "\n" << TableOne.ptr[i];
//	}

	if( start_from_file )
	{
		// setup C4
		this->allocate_C4(ntype, lmax_in, ucell.nmax, Ecut_number, ucell);
		// check tolerence
		this->readin_C4("INPUTs", ntype, ecut, rcut, Ecut_number, tolerence, ucell);
#ifdef __MPI
		Parallel_Common::bcast_double( C4.ptr, C4.getSize() );
#endif
		this->init_Faln(ntype, lmax_in, ucell.nmax, Ecut_number, ucell);
	}

	return;
}

double Bessel_Basis::Polynomial_Interpolation2
	(const int &l, const int &ie, const double &gnorm)const
{
	const double position =  gnorm / this->Dk;
	const int iq = static_cast<int>(position);
	/*
	if(iq >= kmesh-4)
	{
		std::cout << "\n iq = " << iq;
		std::cout << "\n kmesh = " << kmesh;
		ModuleBase::QUIT();
	}
	*/
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

double Bessel_Basis::Polynomial_Interpolation(
	const int &it, const int &l, const int &ic, const double &gnorm)const
{
	const double position =  gnorm / this->Dk;
	const int iq = static_cast<int>(position);
	assert(iq < kmesh-4);
	const double x0 = position - static_cast<double>(iq);
	const double x1 = 1.0 - x0;
	const double x2 = 2.0 - x0;
	const double x3 = 3.0 - x0;
	const double y=
		this->Faln(it, l, ic, iq) * x1 * x2 * x3 / 6.0 +
        this->Faln(it, l, ic, iq) * x0 * x2 * x3 / 2.0 -
        this->Faln(it, l, ic, iq) * x1 * x0 * x3 / 2.0 +
        this->Faln(it, l, ic, iq) * x1 * x2 * x0 / 6.0 ;
	return y;
}

void Bessel_Basis::init_Faln(
	const int &ntype,
	const int &lmax,
	const int &nmax,
	const int &ecut_number,
	const UnitCell& ucell)
{
	ModuleBase::TITLE("Bessel_Basis","init_Faln");
	ModuleBase::timer::tick("Spillage","init_Faln");
	assert( this->kmesh > 0);

	this->Faln.create(ntype, lmax+1, nmax, this->kmesh);

	this->nwfc = 0;
	for(int it=0; it<ntype; it++)
	{
		for(int il=0; il<ucell.atoms[it].nwl+1; il++)
		{
			for(int in=0; in<ucell.atoms[it].l_nchi[il]; in++)
			{
				for(int ie=0; ie<ecut_number; ie++)
				{
					for(int ik=0; ik< this->kmesh; ik++)
					{
						this->Faln(it, il, in, ik) += this->C4(it, il, in, ie) * this->TableOne(il, ie, ik);
					}
				}
				nwfc+=2*il+1;
			}
		}
	}
	ModuleBase::GlobalFunc::OUT("nwfc = ",nwfc);

	ModuleBase::timer::tick("Spillage","init_Faln");
	return;
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
	ModuleBase::TITLE("Bessel_Basis","init_TableOne");
	ModuleBase::timer::tick("Spillage","TableONe");
	// check
	assert(ecutwfc > 0.0);
	assert(dr > 0.0);
	assert(dk > 0.0);

	// init kmesh
	this->kmesh = static_cast<int>(sqrt(ecutwfc) / dk) +1 + 4;
	if (kmesh % 2 == 0)++kmesh;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "kmesh",kmesh);
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "dk",dk);

	// init Table One
	this->TableOne.create(lmax+1, ecut_number, kmesh);

	// init rmesh
	int rmesh = static_cast<int>( rcut / dr ) + 4;
    if (rmesh % 2 == 0) ++rmesh;
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "rmesh",rmesh);
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "dr",dr);

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

	//caoyu add 2021-3-10
	//=========output .orb format=============
	std::stringstream ss;
	ss << GlobalV::global_out_dir << "jle.orb";
	std::ofstream ofs(ss.str().c_str());
	ofs << "---------------------------------------------------------------------------"<< std::endl;
	ofs << std::setiosflags(ios::left) << std::setw(28) << "Energy Cutoff(Ry)" << ecut << std::endl;
	ofs << std::setiosflags(ios::left) << std::setw(28) << "Radius Cutoff(a.u.)" << rcut << std::endl;
	ofs << std::setiosflags(ios::left) << std::setw(28) << "Lmax" << lmax << std::endl;
	for (int l = 0; l < lmax + 1; l++)
	{
		switch (l)
		{
			case 0:
			ofs << std::setiosflags(ios::left) << std::setw(28) << "Number of Sorbitals-->" << ecut_number << std::endl;
			break;
			case 1:
			ofs << std::setiosflags(ios::left) << std::setw(28) << "Number of Porbitals-->" << ecut_number << std::endl;
			break;
			case 2:
			ofs << std::setiosflags(ios::left) << std::setw(28) << "Number of Dorbitals-->" << ecut_number << std::endl;
			break;
			case 3:
			ofs << std::setiosflags(ios::left) << std::setw(28) << "Number of Forbitals-->" << ecut_number << std::endl;
			break;
			default:
			ofs << std::setiosflags(ios::left) << std::setw(28) << "Number of Gorbitals-->" << ecut_number << std::endl;
		}
	}
	ofs << "---------------------------------------------------------------------------"<< std::endl;
	ofs << "SUMMARY END" << std::endl << std::endl;
	ofs << std::setiosflags(ios::left) << std::setw(28) << "Mesh" << rmesh << std::endl;
	ofs << std::setiosflags(ios::left) << std::setw(28) << "dr" << dr << std::endl ;
	//=========output	 .orb format=============

	// init eigenvalue of Jl
	for(int l=0; l<lmax+1; l++)
	{
		ModuleBase::GlobalFunc::ZEROS(en, ecut_number);
		ModuleBase::GlobalFunc::ZEROS(jle, rmesh);
		ModuleBase::GlobalFunc::ZEROS(jlk, rmesh);

		// calculate eigenvalue for l
		ModuleBase::Sphbes::Spherical_Bessel_Roots(ecut_number, l, tolerence, en, rcut);
//		for (int ie=0; ie<ecut_number; ie++)
//		{
//			std::cout << "\n en[" << ie << "]=" << en[ie];
//		}

		// for each eigenvalue
		for (int ie=0; ie<ecut_number; ie++)
		{
			// calculate J_{l}( en[ir]*r)
			ModuleBase::Sphbes::Spherical_Bessel(rmesh, r, en[ie], l, jle);

			//caoyu add 2021-3-10
			//=========output .orb format=============
			ofs << std::setiosflags(ios::right) << std::setw(20) << "Type"<< std::setw(20) << "L" << std::setw(20) << "N" << std::endl;
			ofs << std::setiosflags(ios::right) << std::setw(20) << "0"<< std::setw(20) << l << std::setw(20) << ie << std::endl;
			for (int ir = 0; ir < rmesh; ir++)
			{
				ofs << std::setiosflags(ios::scientific)
				<< std::setprecision(12) << jle[ir]<< " "; if ((ir+1) % 4 == 0) ofs << std::endl;
			}
			ofs << std::endl;
			//=========output .orb format=============

			for(int ir=0; ir<rmesh; ir++)
			{
				jle[ir] = jle[ir] * r[ir] * r[ir];
			}

			//====== output ========
//			std::stringstream ss;
//			ss << GlobalV::global_out_dir << l << "." << ie << ".txt";
//			std::ofstream ofs(ss.str().c_str());

//			for(int ir=0; ir<rmesh; ir++) ofs << r[ir] << " " << jle[ir] << " " << jle[ir]*g[ir] << std::endl;

//			ofs.close();
			//====== output ========

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
				ModuleBase::Sphbes::Spherical_Bessel(rmesh, r, ik*dk, l, jlk);

				// calculate the function will be integrated
				for(int ir=0; ir<rmesh; ir++)
				{
					function[ir] = jle[ir] * jlk[ir];
				}

				// make table value
				ModuleBase::Integral::Simpson_Integral(rmesh, function, rab, this->TableOne(l, ie, ik) );
			}

		}// end ie
	}// end ;

	if (ofs)
	{
		ofs.close();	//caoyu add 2020-3-10
	}

	delete[] en;
	delete[] jle;
	delete[] jlk;
	delete[] rab;
	delete[] g;
	delete[] r;
	delete[] function;
	ModuleBase::timer::tick("Spillage","TableONe");
	return;
}

void Bessel_Basis::readin_C4(
	const std::string &name,
	const int &ntype,
	const int &ecut,
	const int &rcut,
	const int &ecut_number,
	const double &tolerence,
	const UnitCell& ucell)
{
	ModuleBase::TITLE("Bessel_Basis","readin_C4");

	if(GlobalV::MY_RANK != 0) return;

	std::ifstream ifs( name.c_str() );

	if(!ifs)
	{
		GlobalV::ofs_warning << " File name : " << name << std::endl;
		ModuleBase::WARNING_QUIT("Bessel_Basis::readin_C4","Can not find file.");
	}

	if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "<FILE>"))
	{
		// mohan modify 2009-11-29
		for (int it = 0; it < ntype; it++)
		{
			std::string filec4;
			ifs >> filec4;
			for(int il=0; il< ucell.atoms[it].nwl+1; il++)
			{
				for(int in=0; in< ucell.atoms[it].l_nchi[il]; in++)
				{
					//for tests
					//std::cout << "\n" << std::setw(5) << it << std::setw(5) << il << std::setw(5) << in;
					//std::cout << "\n file=" << filec4;
					std::ifstream inc4( filec4.c_str() );

					if(!inc4)
					{
						GlobalV::ofs_warning << " File name : " << filec4 << std::endl;
						ModuleBase::WARNING_QUIT("Bessel_Basis::readin_C4","Can not find file.");
					}

					if(ModuleBase::GlobalFunc::SCAN_BEGIN(inc4, "<INPUTS>"))
					{
						double tmp_ecut;
						double tmp_rcut;
						double tmp_enumber;
						double tmp_tolerence;
						ModuleBase::GlobalFunc::READ_VALUE( inc4, tmp_ecut);
						ModuleBase::GlobalFunc::READ_VALUE( inc4, tmp_rcut);
						ModuleBase::GlobalFunc::READ_VALUE( inc4, tmp_enumber);
						ModuleBase::GlobalFunc::READ_VALUE( inc4, tmp_tolerence);
						assert( tmp_ecut == this->ecut );
						assert( tmp_rcut == this->rcut );
						assert( tmp_enumber == this->Ecut_number);
						assert( tmp_tolerence == this->tolerence );
					}

					bool find = false;
					if(ModuleBase::GlobalFunc::SCAN_BEGIN(inc4, "<C4>"))
					{
						int total_nchi = 0;
						ModuleBase::GlobalFunc::READ_VALUE(inc4, total_nchi);

						for(int ichi=0; ichi<total_nchi; ichi++)
						{
							std::string title1, title2, title3;
							inc4 >> title1 >> title2 >> title3;

							int tmp_type, tmp_l, tmp_n;
							inc4 >> tmp_type >> tmp_l >> tmp_n;
							//std::cout << "\n Find T=" << tmp_type << " L=" << tmp_l << " N=" << tmp_n;

							if(tmp_l == il && tmp_n == in)
							//if(tmp_type == it && tmp_l == il && tmp_n == in) // mohan modify 2009-11-29
							{
								find = true;
								for(int ie=0; ie<ecut_number; ie++)
								{
									inc4 >> this->C4(it, il, in, ie);
									// for tests
									//std::cout << "\n" << std::setw(5) << ie << std::setw(25) << this->C4(it, il, in, ie);
								}
							}
							else
							{
								double no_use_c4;
								for(int ie=0; ie<ecut_number; ie++)
								{
									inc4 >> no_use_c4;
								}
							}
							if(find) break;
						}
					}
					if(!find)
					{
						std::cout << "\n T=" << it << " L=" << il << " N=" << in;
						ModuleBase::WARNING_QUIT("Bessel_Basis::readin_C4","Can't find needed c4!");
					}
					inc4.close();
				}
			}
		}
		ModuleBase::GlobalFunc::SCAN_END(ifs, "</FILE>");
	}
	ifs.close();
	return;
}

void Bessel_Basis::allocate_C4(
	const int &ntype,
	const int &lmax,
	const int &nmax,
	const int &ecut_number,
	const UnitCell& ucell)
{
	ModuleBase::TITLE("Bessel_Basis","allocate_C4");

	this->C4.create(ntype, lmax+1, nmax, ecut_number);

	for(int it=0; it<ntype; it++)
	{
		for(int il=0; il<ucell.atoms[it].nwl+1; il++)
		{
			for(int in=0; in<ucell.atoms[it].l_nchi[il]; in++)
			{
				for(int ie=0; ie<ecut_number; ie++)
				{
					this->C4(it, il, in, ie) = 1.0;
				}
			}
		}
	}
	return;
}
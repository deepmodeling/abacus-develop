#include "forces.h"
#include "global.h"
#include "vdwd2.h"
#include "vdwd3.h"				  
#include "../module_symmetry/symmetry.h"
// new
#include "../module_xc/xc_functional.h"
#include "../module_base/math_integral.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

double Forces::output_acc = 1.0e-8; // (Ryd/angstrom).	

Forces::Forces()
{
}

Forces::~Forces() {}

#include "../module_base/mathzone.h"
#include "efield.h"
void Forces::init(ModuleBase::matrix& force)
{
	ModuleBase::TITLE("Forces", "init");
	this->nat = GlobalC::ucell.nat;
	force.create(nat, 3);
	
	ModuleBase::matrix forcelc(nat, 3);
	ModuleBase::matrix forceion(nat, 3);
	ModuleBase::matrix forcecc(nat, 3);
	ModuleBase::matrix forcenl(nat, 3);
	ModuleBase::matrix forcescc(nat, 3);
    this->cal_force_loc(forcelc);
    this->cal_force_ew(forceion);
    this->cal_force_nl(forcenl);
	this->cal_force_cc(forcecc);
	this->cal_force_scc(forcescc);

	ModuleBase::matrix stress_vdw_pw;//.create(3,3);
    ModuleBase::matrix force_vdw;
    force_vdw.create(nat, 3);
	if(GlobalC::vdwd2_para.flag_vdwd2)													//Peize Lin add 2014.04.03, update 2021.03.09
	{
        Vdwd2 vdwd2(GlobalC::ucell,GlobalC::vdwd2_para);
		vdwd2.cal_force();
		for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		{
			force_vdw(iat,0) = vdwd2.get_force()[iat].x;
			force_vdw(iat,1) = vdwd2.get_force()[iat].y;
			force_vdw(iat,2) = vdwd2.get_force()[iat].z;
		}
		if(GlobalV::TEST_FORCE)
		{
			Forces::print("VDW      FORCE (Ry/Bohr)", force_vdw);
		}
	}
	else if(GlobalC::vdwd3_para.flag_vdwd3)													//jiyy add 2019-05-18, update 2021-05-02
	{
        Vdwd3 vdwd3(GlobalC::ucell,GlobalC::vdwd3_para);
		vdwd3.cal_force();
		for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		{
			force_vdw(iat,0) = vdwd3.get_force()[iat].x;
			force_vdw(iat,1) = vdwd3.get_force()[iat].y;
			force_vdw(iat,2) = vdwd3.get_force()[iat].z;
		}
		if(GlobalV::TEST_FORCE)
		{
			Forces::print("VDW      FORCE (Ry/Bohr)", force_vdw);
		}
	}
    //impose total force = 0
    int iat = 0;

	ModuleBase::matrix force_e;
	if(GlobalV::EFIELD)
	{
		force_e.create(GlobalC::ucell.nat, 3);
		Efield::compute_force(force_e);
	}
	
	for (int ipol = 0; ipol < 3; ipol++)
	{
		double sum = 0.0;
		iat = 0;

		for (int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
			{
				force(iat, ipol) =
					forcelc(iat, ipol)
					+ forceion(iat, ipol)
					+ forcenl(iat, ipol)
					+ forcecc(iat, ipol)
					+ forcescc(iat, ipol);

				if(GlobalC::vdwd2_para.flag_vdwd2 || GlobalC::vdwd3_para.flag_vdwd3)		//linpz and jiyy added vdw force, modified by zhengdy
				{
                    force(iat, ipol) += force_vdw(iat, ipol);
                }																										   
					
				if(GlobalV::EFIELD)
				{
					force(iat,ipol) = force(iat, ipol) + force_e(iat, ipol);
				}

				sum += force(iat, ipol);

				iat++;
			}
		}

		double compen = sum / GlobalC::ucell.nat;
		for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		{
			force(iat, ipol) = force(iat, ipol) - compen;
		}	
	}
	
	if(ModuleSymmetry::Symmetry::symm_flag)
	{
		double *pos;
		double d1,d2,d3;
		pos = new double[GlobalC::ucell.nat*3];
		ModuleBase::GlobalFunc::ZEROS(pos, GlobalC::ucell.nat*3);
		int iat = 0;
		for(int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			//Atom* atom = &GlobalC::ucell.atoms[it];
			for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
			{
				pos[3*iat  ] = GlobalC::ucell.atoms[it].taud[ia].x ;
				pos[3*iat+1] = GlobalC::ucell.atoms[it].taud[ia].y ;
				pos[3*iat+2] = GlobalC::ucell.atoms[it].taud[ia].z;
				for(int k=0; k<3; ++k)
				{
					GlobalC::symm.check_translation( pos[iat*3+k], -floor(pos[iat*3+k]));
					GlobalC::symm.check_boundary( pos[iat*3+k] );
				}
				iat++;				
			}
		}
		
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			ModuleBase::Mathzone::Cartesian_to_Direct(force(iat,0),force(iat,1),force(iat,2),
                                        GlobalC::ucell.a1.x, GlobalC::ucell.a1.y, GlobalC::ucell.a1.z,
                                        GlobalC::ucell.a2.x, GlobalC::ucell.a2.y, GlobalC::ucell.a2.z,
                                        GlobalC::ucell.a3.x, GlobalC::ucell.a3.y, GlobalC::ucell.a3.z,
                                        d1,d2,d3);
			
			force(iat,0) = d1;force(iat,1) = d2;force(iat,2) = d3;
		}
		GlobalC::symm.force_symmetry(force , pos, GlobalC::ucell);
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			ModuleBase::Mathzone::Direct_to_Cartesian(force(iat,0),force(iat,1),force(iat,2),
                                        GlobalC::ucell.a1.x, GlobalC::ucell.a1.y, GlobalC::ucell.a1.z,
                                        GlobalC::ucell.a2.x, GlobalC::ucell.a2.y, GlobalC::ucell.a2.z,
                                        GlobalC::ucell.a3.x, GlobalC::ucell.a3.y, GlobalC::ucell.a3.z,
                                        d1,d2,d3);
			force(iat,0) = d1;force(iat,1) = d2;force(iat,2) = d3;
		}
		// std::cout << "nrotk =" << GlobalC::symm.nrotk << std::endl;
		delete[] pos;
		
	}

 	GlobalV::ofs_running << std::setiosflags(ios::fixed) << std::setprecision(6) << std::endl;
	/*if(GlobalV::TEST_FORCE)
	{
		Forces::print("LOCAL    FORCE (Ry/Bohr)", forcelc);
		Forces::print("NONLOCAL FORCE (Ry/Bohr)", forcenl);
		Forces::print("NLCC     FORCE (Ry/Bohr)", forcecc);
		Forces::print("ION      FORCE (Ry/Bohr)", forceion);
		Forces::print("SCC      FORCE (Ry/Bohr)", forcescc);
		if(GlobalV::EFIELD) Forces::print("EFIELD   FORCE (Ry/Bohr)", force_e);
	}*/
	
/*
	Forces::print("   TOTAL-FORCE (Ry/Bohr)", force);
	
	if(INPUT.out_force)                                                   // pengfei 2016-12-20
	{
		std::ofstream ofs("FORCE.dat");
		if(!ofs)
		{
			std::cout << "open FORCE.dat error !" <<std::endl;
		}
		for(int iat=0; iat<GlobalC::ucell.nat; iat++)
		{
			ofs << "   " << force(iat,0)*ModuleBase::Ry_to_eV / 0.529177 
				<< "   " << force(iat,1)*ModuleBase::Ry_to_eV / 0.529177 
				<< "   " << force(iat,2)*ModuleBase::Ry_to_eV / 0.529177 << std::endl;
		}
		ofs.close();
	}
*/
		
	// output force in unit eV/Angstrom
	GlobalV::ofs_running << std::endl;

	if(GlobalV::TEST_FORCE)
	{
		Forces::print("LOCAL    FORCE (eV/Angstrom)", forcelc,0);
		Forces::print("NONLOCAL FORCE (eV/Angstrom)", forcenl,0);
		Forces::print("NLCC     FORCE (eV/Angstrom)", forcecc,0);
		Forces::print("ION      FORCE (eV/Angstrom)", forceion,0);
		Forces::print("SCC      FORCE (eV/Angstrom)", forcescc,0);
		if(GlobalV::EFIELD) Forces::print("EFIELD   FORCE (eV/Angstrom)", force_e,0);
	}
	Forces::print("   TOTAL-FORCE (eV/Angstrom)", force,0);

    return;
}

void Forces::print_to_files(std::ofstream &ofs, const std::string &name, const ModuleBase::matrix &f)
{
    int iat = 0;
    ofs << " " << name;
    ofs << std::setprecision(8);
	//ofs << std::setiosflags(ios::showpos);
   
	double fac = ModuleBase::Ry_to_eV / 0.529177;// (eV/A)

	if(GlobalV::TEST_FORCE)
	{
		std::cout << std::setiosflags(ios::showpos);
		std::cout << " " << name;
		std::cout << std::setprecision(8);
	}

    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
            ofs << " " << std::setw(5) << it
            << std::setw(8) << ia+1
            << std::setw(20) << f(iat, 0)*fac
            << std::setw(20) << f(iat, 1)*fac
            << std::setw(20) << f(iat, 2)*fac << std::endl;
			
			if(GlobalV::TEST_FORCE)
			{
            	std::cout << " " << std::setw(5) << it
            	<< std::setw(8) << ia+1
            	<< std::setw(20) << f(iat, 0)*fac
            	<< std::setw(20) << f(iat, 1)*fac
            	<< std::setw(20) << f(iat, 2)*fac << std::endl;
			}
            iat++;
        }
    }

	GlobalV::ofs_running << std::resetiosflags(ios::showpos);
	std::cout << std::resetiosflags(ios::showpos);
    return;
}



void Forces::print(const std::string &name, const ModuleBase::matrix &f, bool ry)
{
	ModuleBase::GlobalFunc::NEW_PART(name);

	GlobalV::ofs_running << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
	GlobalV::ofs_running << std::setiosflags(ios::showpos);
    GlobalV::ofs_running << std::setprecision(8);

	const double fac = ModuleBase::Ry_to_eV / 0.529177;
	
	if(GlobalV::TEST_FORCE)
	{
		std::cout << " --------------- " << name << " ---------------" << std::endl;
		std::cout << " " << std::setw(8) << "atom" << std::setw(15) << "x" << std::setw(15) << "y" << std::setw(15) << "z" << std::endl;
		std::cout << std::setiosflags(ios::showpos);
		std::cout << std::setprecision(6);
		std::cout << std::setiosflags(ios::fixed);
	}

    int iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
			std::stringstream ss;
			ss << GlobalC::ucell.atoms[it].label << ia+1;

			if(ry) // output Rydberg Unit
			{
				GlobalV::ofs_running << " " << std::setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,0);
				else GlobalV::ofs_running << std::setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,1);
				else GlobalV::ofs_running << std::setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,2);
				else GlobalV::ofs_running << std::setw(15) << "0";
				GlobalV::ofs_running << std::endl;
			}
			else
			{
				GlobalV::ofs_running << " " << std::setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,0)*fac;
				else GlobalV::ofs_running << std::setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,1)*fac;
				else GlobalV::ofs_running << std::setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) GlobalV::ofs_running << std::setw(15) << f(iat,2)*fac;
				else GlobalV::ofs_running << std::setw(15) << "0";
				GlobalV::ofs_running << std::endl;
			}

			if(GlobalV::TEST_FORCE && ry)
			{
				std::cout << " " << std::setw(8) << ss.str();
                std::cout<<fixed;
				if( abs(f(iat,0)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,0);
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,1);
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,2);
				else std::cout << std::setw(15) << "0";
				std::cout << std::endl;
			}
			else if (GlobalV::TEST_FORCE)
			{
				std::cout << " " << std::setw(8) << ss.str();
                std::cout<<fixed;
				if( abs(f(iat,0)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,0)*fac;
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,1)*fac;
				else std::cout << std::setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) std::cout << std::setw(15) << f(iat,2)*fac;
				else std::cout << std::setw(15) << "0";
				std::cout << std::endl;
			}	
				
            iat++;
        }
    }

	GlobalV::ofs_running << std::resetiosflags(ios::showpos);
	std::cout << std::resetiosflags(ios::showpos);
    return;
}


void Forces::cal_force_loc(ModuleBase::matrix& forcelc)
{
	ModuleBase::timer::tick("Forces","cal_force_loc");

    std::complex<double> *aux = new std::complex<double>[GlobalC::pw.nrxx];
    ModuleBase::GlobalFunc::ZEROS(aux, GlobalC::pw.nrxx);

    // now, in all pools , the charge are the same,
    // so, the force calculated by each pool is equal.
    
	for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
        	aux[ir] += std::complex<double>( GlobalC::CHR.rho[is][ir], 0.0 );
		}
	}

	// to G space.
    GlobalC::pw.FFT_chg.FFT3D(aux, -1);

    int gstart_here = GlobalC::pw.gstart;
    if (GlobalC::pw.ggs[0] != 0) gstart_here = 0;

//  GlobalV::ofs_running << "\n ggs = " << GlobalC::pw.ggs[0];
//  GlobalV::ofs_running << "\n gstart_here = " << gstart_here;
    int iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
            for (int ig = gstart_here; ig < GlobalC::pw.ngmc; ig++)
            {
                const double phase = ModuleBase::TWO_PI * (GlobalC::pw.get_G_cartesian(ig) * GlobalC::ucell.atoms[it].tau[ia]);
                const double factor = GlobalC::ppcell.vloc(it, GlobalC::pw.ig2ngg[ig]) *
									  ( cos(phase) * aux[ GlobalC::pw.ig2fftc[ig] ].imag()
                                      + sin(phase) * aux[ GlobalC::pw.ig2fftc[ig] ].real()); 
                forcelc(iat, 0) += GlobalC::pw.get_G_cartesian_projection(ig, 0) * factor;
                forcelc(iat, 1) += GlobalC::pw.get_G_cartesian_projection(ig, 1) * factor;
                forcelc(iat, 2) += GlobalC::pw.get_G_cartesian_projection(ig, 2) * factor;
            }
            for (int ipol = 0;ipol < 3;ipol++)
            {
                forcelc(iat, ipol) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
            }
            ++iat;
        }
    }
    //this->print(GlobalV::ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_double_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
	ModuleBase::timer::tick("Forces","cal_force_loc");
    return;
}

#include "H_Ewald_pw.h"
void Forces::cal_force_ew(ModuleBase::matrix& forceion)
{
	ModuleBase::timer::tick("Forces","cal_force_ew");

    double fact = 2.0;
    std::complex<double> *aux = new std::complex<double> [GlobalC::pw.ngmc];
    ModuleBase::GlobalFunc::ZEROS(aux, GlobalC::pw.ngmc);

    int gstart = GlobalC::pw.gstart;

    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ig = gstart; ig < GlobalC::pw.ngmc; ig++)
        {
            aux[ig] += static_cast<double>(GlobalC::ucell.atoms[it].zv) * conj(GlobalC::pw.strucFac(it, ig));
        }
    }

	// calculate total ionic charge
    double charge = 0.0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        charge += GlobalC::ucell.atoms[it].na * GlobalC::ucell.atoms[it].zv;//mohan modify 2007-11-7
    }
	
	double alpha = 1.1;
	double upperbound ;
    do
    {
        alpha -= 0.10;
        // choose alpha in order to have convergence in the sum over G
        // upperbound is a safe upper bound for the error in the sum over G

        if (alpha <= 0.0)
        {
            ModuleBase::WARNING_QUIT("ewald","Can't find optimal alpha.");
        }
        upperbound = 2.0 * charge * charge * sqrt(2.0 * alpha / ModuleBase::TWO_PI) *
                     erfc(sqrt(GlobalC::ucell.tpiba2 * GlobalC::pw.ggchg / 4.0 / alpha));
    }
    while (upperbound > 1.0e-6);
//	std::cout << " GlobalC::en.alpha = " << alpha << std::endl;
//	std::cout << " upperbound = " << upperbound << std::endl;
	


    for (int ig = gstart; ig < GlobalC::pw.ngmc; ig++)
    {
        if(GlobalC::pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            aux[ig] *= exp(-1.0 * GlobalC::pw.gg[ig] * GlobalC::ucell.tpiba2 / alpha / 4.0) / (GlobalC::pw.gg[ig] * GlobalC::ucell.tpiba2);
        }
    }

    int iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
        {
            for (int ig = gstart; ig < GlobalC::pw.ngmc; ig++)
            {
                const double arg = ModuleBase::TWO_PI * (GlobalC::pw.get_G_cartesian(ig) * GlobalC::ucell.atoms[it].tau[ia]);
                double sumnb =  -cos(arg) * aux[ig].imag() + sin(arg) * aux[ig].real();
                forceion(iat, 0) += GlobalC::pw.get_G_cartesian_projection(ig, 0) * sumnb;
                forceion(iat, 1) += GlobalC::pw.get_G_cartesian_projection(ig, 1) * sumnb;
                forceion(iat, 2) += GlobalC::pw.get_G_cartesian_projection(ig, 2) * sumnb;
            }
            for (int ipol = 0;ipol < 3;ipol++)
            {
                forceion(iat, ipol) *= GlobalC::ucell.atoms[it].zv * ModuleBase::e2 * GlobalC::ucell.tpiba * ModuleBase::TWO_PI / GlobalC::ucell.omega * fact;
            }

	//		std::cout << " atom" << iat << std::endl;
	//		std::cout << std::setw(15) << forceion(iat, 0) << std::setw(15) << forceion(iat,1) << std::setw(15) << forceion(iat,2) << std::endl; 
            iat++;
        }
    }
    delete [] aux;


	// means that the processor contains G=0 term.
    if (gstart == 1)
    {
        double rmax = 5.0 / (sqrt(alpha) * GlobalC::ucell.lat0);
        int nrm = 0;
		
        //output of rgen: the number of vectors in the sphere
        const int mxr = 50;
        // the maximum number of R vectors included in r
        ModuleBase::Vector3<double> *r  = new ModuleBase::Vector3<double>[mxr];
        double *r2 = new double[mxr];
		ModuleBase::GlobalFunc::ZEROS(r2, mxr);
        int *irr = new int[mxr];
		ModuleBase::GlobalFunc::ZEROS(irr, mxr);
        // the square modulus of R_j-tau_s-tau_s'

		int iat1 = 0;
        for (int T1 = 0; T1 < GlobalC::ucell.ntype; T1++)
        {
			Atom* atom1 = &GlobalC::ucell.atoms[T1]; 
            for (int I1 = 0; I1 < atom1->na; I1++)
            {
				int iat2 = 0; // mohan fix bug 2011-06-07
                for (int T2 = 0; T2 < GlobalC::ucell.ntype; T2++)
                {
                    for (int I2 = 0; I2 < GlobalC::ucell.atoms[T2].na; I2++)
                    {
                        if (iat1 != iat2)
                        {
                            ModuleBase::Vector3<double> d_tau = GlobalC::ucell.atoms[T1].tau[I1] - GlobalC::ucell.atoms[T2].tau[I2];
                            H_Ewald_pw::rgen(d_tau, rmax, irr, GlobalC::ucell.latvec, GlobalC::ucell.G, r, r2, nrm);

                            for (int n = 0;n < nrm;n++)
                            {
								const double rr = sqrt(r2[n]) * GlobalC::ucell.lat0;

                                double factor = GlobalC::ucell.atoms[T1].zv * GlobalC::ucell.atoms[T2].zv * ModuleBase::e2 / (rr * rr)
                                                * (erfc(sqrt(alpha) * rr) / rr
                                    + sqrt(8.0 * alpha / ModuleBase::TWO_PI) * exp(-1.0 * alpha * rr * rr)) * GlobalC::ucell.lat0;

								forceion(iat1, 0) -= factor * r[n].x;
                                forceion(iat1, 1) -= factor * r[n].y;
                                forceion(iat1, 2) -= factor * r[n].z;

//								std::cout << " r.z=" << r[n].z << " r2=" << r2[n] << std::endl;
						//		std::cout << " " << iat1 << " " << iat2 << " n=" << n
						//		 << " rn.z=" << r[n].z 
						//		 << " r2=" << r2[n] << " rr=" << rr << " fac=" << factor << " force=" << forceion(iat1,2) 
						//		 << " new_part=" << factor*r[n].z <<  std::endl;
                            }
                        }

                        ++iat2;
                    }
                }//atom b

//				std::cout << " atom" << iat1 << std::endl;
//				std::cout << std::setw(15) << forceion(iat1, 0) << std::setw(15) << forceion(iat1,1) << std::setw(15) << forceion(iat1,2) << std::endl; 

                ++iat1;
            }
        }//atom a
    }

    Parallel_Reduce::reduce_double_pool(forceion.c, forceion.nr * forceion.nc);

    //this->print(GlobalV::ofs_running, "ewald forces", forceion);

	ModuleBase::timer::tick("Forces","cal_force_ew");

    return;
}

void Forces::cal_force_cc(ModuleBase::matrix& forcecc)
{
	// recalculate the exchange-correlation potential.
	
    ModuleBase::matrix v(GlobalV::NSPIN,GlobalC::pw.nrxx);

	if(XC_Functional::get_func_type() == 3)
	{
#ifdef USE_LIBXC
    	const auto etxc_vtxc_v = XC_Functional::v_xc_meta(
            GlobalC::pw.nrxx, GlobalC::pw.ncxyz, GlobalC::ucell.omega,
            GlobalC::CHR.rho, GlobalC::CHR.rho_core, GlobalC::CHR.kin_r);
        
        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
        v = std::get<2>(etxc_vtxc_v);
#else
        ModuleBase::WARNING_QUIT("cal_force_cc","to use mGGA, compile with LIBXC");
#endif
	}
	else
	{	
    	const auto etxc_vtxc_v = XC_Functional::v_xc(
            GlobalC::pw.nrxx, GlobalC::pw.ncxyz, GlobalC::ucell.omega,
            GlobalC::CHR.rho, GlobalC::CHR.rho_core);
        
        GlobalC::en.etxc = std::get<0>(etxc_vtxc_v);
        GlobalC::en.vtxc = std::get<1>(etxc_vtxc_v);
	    v = std::get<2>(etxc_vtxc_v);
	}

	const ModuleBase::matrix vxc = v;
    std::complex<double> * psiv = new std::complex<double> [GlobalC::pw.nrxx];
    ModuleBase::GlobalFunc::ZEROS(psiv, GlobalC::pw.nrxx);
    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
        for (int ir = 0;ir < GlobalC::pw.nrxx;ir++)
        {
            psiv[ir] = std::complex<double>(vxc(0, ir),  0.0);
        }
    }
    else
    {
        for (int ir = 0;ir < GlobalC::pw.nrxx;ir++)
        {
            psiv[ir] = 0.5 * (vxc(0 ,ir) + vxc(1, ir));
        }
    }

	// to G space
    GlobalC::pw.FFT_chg.FFT3D(psiv, -1);

    //psiv contains now Vxc(G)
    double * rhocg = new double [GlobalC::pw.nggm];
    ModuleBase::GlobalFunc::ZEROS(rhocg, GlobalC::pw.nggm);
    int iat = 0;
    for (int T1 = 0;T1 < GlobalC::ucell.ntype;T1++)
    {
        if (GlobalC::ucell.atoms[T1].nlcc)
        {
            //call drhoc
            GlobalC::CHR.non_linear_core_correction(
                GlobalC::ppcell.numeric,
                GlobalC::ucell.atoms[T1].msh,
                GlobalC::ucell.atoms[T1].r,
                GlobalC::ucell.atoms[T1].rab,
                GlobalC::ucell.atoms[T1].rho_atc,
                rhocg);


			std::complex<double> ipol0, ipol1, ipol2;
            for (int I1 = 0;I1 < GlobalC::ucell.atoms[T1].na;I1++)
            {
                for (int ig = GlobalC::pw.gstart; ig < GlobalC::pw.ngmc; ig++)
                {
                    const double arg = ModuleBase::TWO_PI * (GlobalC::pw.get_G_cartesian_projection(ig, 0) * GlobalC::ucell.atoms[T1].tau[I1].x + 
                        GlobalC::pw.get_G_cartesian_projection(ig, 1) * GlobalC::ucell.atoms[T1].tau[I1].y + 
                        GlobalC::pw.get_G_cartesian_projection(ig, 2) * GlobalC::ucell.atoms[T1].tau[I1].z);

                    ipol0 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocg[GlobalC::pw.ig2ngg[ig]] * 
                        GlobalC::pw.get_G_cartesian_projection(ig, 0) * conj(psiv[GlobalC::pw.ig2fftc[ig]]) * std::complex<double>(sin(arg), cos(arg));
                    forcecc(iat, 0) +=  ipol0.real();

                    ipol1 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocg[GlobalC::pw.ig2ngg[ig]] * 
                        GlobalC::pw.get_G_cartesian_projection(ig, 1) * conj(psiv[GlobalC::pw.ig2fftc[ig]]) * std::complex<double>(sin(arg), cos(arg));
                    forcecc(iat, 1) += ipol1.real();

                    ipol2 = GlobalC::ucell.tpiba * GlobalC::ucell.omega * rhocg[GlobalC::pw.ig2ngg[ig]] * 
                        GlobalC::pw.get_G_cartesian_projection(ig, 2) * conj(psiv[GlobalC::pw.ig2fftc[ig]]) * std::complex<double>(sin(arg), cos(arg));

                    forcecc(iat, 2) += ipol2.real();
                }
                ++iat;
            }
        }
        else{
            iat += GlobalC::ucell.atoms[T1].na;
        }
    }
    assert(iat == GlobalC::ucell.nat);
    delete [] rhocg;
	delete [] psiv; // mohan fix bug 2012-03-22
    Parallel_Reduce::reduce_double_pool(forcecc.c, forcecc.nr * forcecc.nc); //qianrui fix a bug for kpar > 1
	return;
}

#include "../module_base/complexarray.h"
#include "../module_base/complexmatrix.h"
void Forces::cal_force_nl(ModuleBase::matrix& forcenl)
{
	ModuleBase::TITLE("Forces","cal_force_nl");
	ModuleBase::timer::tick("Forces","cal_force_nl");

    const int nkb = GlobalC::ppcell.nkb;
	if(nkb == 0) return; // mohan add 2010-07-25
	
	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	ModuleBase::ComplexArray dbecp( nkb, GlobalV::NBANDS, 3);
    ModuleBase::ComplexMatrix becp( nkb, GlobalV::NBANDS);
    
	
	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ModuleBase::ComplexMatrix vkb1( nkb, GlobalC::wf.npwx );

    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
        if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
        GlobalC::wf.npw = GlobalC::kv.ngk[ik];
        // generate vkb
        if (GlobalC::ppcell.nkb > 0)
        {
            GlobalC::ppcell.getvnl(ik);
        }

        // get becp according to wave functions and vkb
        // important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        becp.zero_out();
        for (int ib=0; ib<GlobalV::NBANDS; ib++)
        {
            for (int i=0;i<nkb;i++)
            {
                for (int ig=0; ig<GlobalC::wf.npw; ig++)
                {
                    becp(i,ib) += GlobalC::wf.evc[ik](ib,ig)* conj( GlobalC::ppcell.vkb(i,ig) );
                }
            }
        }
        Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);

        //out.printcm_real("becp",becp,1.0e-4);
        // Calculate the derivative of beta,
        // |dbeta> =  -ig * |beta>
        dbecp.zero_out();
        for (int ipol = 0; ipol<3; ipol++)
        {
			for (int i = 0;i < nkb;i++)
			{
				if (ipol==0)
				{
					for (int ig=0; ig<GlobalC::wf.npw; ig++)
                        vkb1(i, ig) = GlobalC::ppcell.vkb(i, ig) * ModuleBase::NEG_IMAG_UNIT * GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(ik, ig), 0);
                }
				if (ipol==1)
				{
					for (int ig=0; ig<GlobalC::wf.npw; ig++)
                        vkb1(i, ig) = GlobalC::ppcell.vkb(i, ig) * ModuleBase::NEG_IMAG_UNIT * GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(ik, ig), 1);
                }
				if (ipol==2)
				{
					for (int ig=0; ig<GlobalC::wf.npw; ig++)
                        vkb1(i, ig) = GlobalC::ppcell.vkb(i, ig) * ModuleBase::NEG_IMAG_UNIT * GlobalC::pw.get_G_cartesian_projection(GlobalC::wf.igk(ik, ig), 2);
                }
			}
            for (int ib=0; ib<GlobalV::NBANDS; ib++)
            {
                ///
                ///only occupied band should be calculated.
                ///
                if(GlobalC::wf.wg(ik, ib) < ModuleBase::threshold_wg) continue;
                for (int i=0; i<nkb; i++)
                {
                    for (int ig=0; ig<GlobalC::wf.npw; ig++)
                    {
                        dbecp(i,ib, ipol) += conj( vkb1(i,ig) ) * GlobalC::wf.evc[ik](ib,ig) ;
                    }
                }
            }
        }// end ipol

//		don't need to reduce here, keep dbecp different in each processor,
//		and at last sum up all the forces.
//		Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//		double *cf = new double[GlobalC::ucell.nat*3];
//		ModuleBase::GlobalFunc::ZEROS(cf, GlobalC::ucell.nat);
		for (int ib=0; ib<GlobalV::NBANDS; ib++)
		{
            ///
			///only occupied band should be calculated.
			///
            if(GlobalC::wf.wg(ik, ib) < ModuleBase::threshold_wg) continue;
			double fac = GlobalC::wf.wg(ik, ib) * 2.0 * GlobalC::ucell.tpiba;
        	int iat = 0;
        	int sum = 0;
			for (int it=0; it<GlobalC::ucell.ntype; it++)
			{
				const int Nprojs = GlobalC::ucell.atoms[it].nh;
				for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
				{
					for (int ip=0; ip<Nprojs; ip++)
					{
						double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip) ;
						const int inkb = sum + ip; 
						//out<<"\n ps = "<<ps;

						for (int ipol=0; ipol<3; ipol++)
						{
							const double dbb = ( conj( dbecp( inkb, ib, ipol) ) * becp( inkb, ib) ).real();
							forcenl(iat, ipol) = forcenl(iat, ipol) - ps * fac * dbb;
							//cf[iat*3+ipol] += ps * fac * dbb;
						}
					}

					//if ( GlobalC::ucell.atoms[it].nbeta > GlobalC::ucell.atoms[it].lmax+1 )    //{zws add 20160110
					//{
					//std::cout << " \n multi-projector force calculation ... " << std::endl;
					for (int ip=0; ip<Nprojs; ip++)
					{
						const int inkb = sum + ip;
						//for (int ip2=0; ip2<Nprojs; ip2++)
						for (int ip2=ip+1; ip2<Nprojs; ip2++)
						{
						//if ( ip != ip2 )
						//{
							const int jnkb = sum + ip2;
							double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip2, ip) ;

							for (int ipol=0; ipol<3; ipol++)
							{
								const double dbb = ( conj( dbecp( inkb, ib, ipol) ) * becp( jnkb, ib)
										+ dbecp( jnkb, ib, ipol) * conj(becp( inkb, ib) ) ).real();
								//const double dbb = ( conj( dbecp( inkb, ib, ipol) ) * becp( jnkb, ib) ).real();
								forcenl(iat, ipol) = forcenl(iat, ipol) - ps * fac * dbb;
								//cf[iat*3+ipol] += ps * fac * dbb;
							}
						//}
						}
					}
					//}    //}zws add 20160110

					++iat;
					sum+=Nprojs;
				}
			} //end it
		} //end band
    }// end ik

    // sum up forcenl from all processors
    Parallel_Reduce::reduce_double_all(forcenl.c, forcenl.nr * forcenl.nc);
//  this->print(GlobalV::ofs_running, "nonlocal forces", forcenl);
	ModuleBase::timer::tick("Forces","cal_force_nl");
    return;
}

void Forces::cal_force_scc(ModuleBase::matrix& forcescc)
{
    std::complex<double>* psic = new std::complex<double> [GlobalC::pw.nrxx];
    ModuleBase::GlobalFunc::ZEROS(psic, GlobalC::pw.nrxx);

    if (GlobalV::NSPIN == 1 || GlobalV::NSPIN == 4)
    {
        for (int i = 0;i < GlobalC::pw.nrxx;i++)
        {
            psic[i] = GlobalC::pot.vnew(0,i);
        }
    }
    else
    {
        int isup = 0;
        int isdw = 1;
        for (int i = 0;i < GlobalC::pw.nrxx;i++)
        {
            psic[i] = (GlobalC::pot.vnew(isup, i) + GlobalC::pot.vnew(isdw, i)) * 0.5;
        }
    }

    int ndm = 0;

    for (int it = 0;it < GlobalC::ucell.ntype;it++)
    {
        if (ndm < GlobalC::ucell.atoms[it].msh)
        {
            ndm = GlobalC::ucell.atoms[it].msh;
        }
    }

    //work space
    double* aux = new double[ndm];
    ModuleBase::GlobalFunc::ZEROS(aux, ndm);

    double* rhocgnt = new double[GlobalC::pw.nggm];
    ModuleBase::GlobalFunc::ZEROS(rhocgnt, GlobalC::pw.nggm);

    GlobalC::pw.FFT_chg.FFT3D(psic, -1);

    double fact = 2.0;
    for (int nt = 0;nt < GlobalC::ucell.ntype;nt++)
    {
//		Here we compute the G.ne.0 term
        const int mesh = GlobalC::ucell.atoms[nt].msh;
        for (int ig = GlobalC::pw.gstart;ig < GlobalC::pw.nggm;ig++)
        {
            const double gx = sqrt(GlobalC::pw.ggs[ig]) * GlobalC::ucell.tpiba;
            for (int ir = 0;ir < mesh;ir++)
            {
                if (GlobalC::ucell.atoms[nt].r[ir] < 1.0e-8)
                {
                    aux[ir] = GlobalC::ucell.atoms[nt].rho_at[ir];
                }
                else
                {
                    const double gxx = gx * GlobalC::ucell.atoms[nt].r[ir];
                    aux[ir] = GlobalC::ucell.atoms[nt].rho_at[ir] * sin(gxx) / gxx;
                }
            }
            ModuleBase::Integral::Simpson_Integral(mesh , aux, GlobalC::ucell.atoms[nt].rab , rhocgnt [ig]);
        }

        int iat = 0;
        for (int it = 0;it < GlobalC::ucell.ntype;it++)
        {
            for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
            {
                if (nt == it)
                {
                    for (int ig = GlobalC::pw.gstart;ig < GlobalC::pw.ngmc;ig++)
                    {
                        const double arg = ModuleBase::TWO_PI * (GlobalC::pw.get_G_cartesian_projection(ig, 0) * GlobalC::ucell.atoms[it].tau[ia].x + 
                            GlobalC::pw.get_G_cartesian_projection(ig, 1) * GlobalC::ucell.atoms[it].tau[ia].y + 
                            GlobalC::pw.get_G_cartesian_projection(ig, 2) * GlobalC::ucell.atoms[it].tau[ia].z);

                        const std::complex<double> cpm = std::complex<double>(sin(arg), cos(arg)) * conj(psic[GlobalC::pw.ig2fftc[ig] ]);

                        forcescc(iat, 0) += fact * rhocgnt[GlobalC::pw.ig2ngg[ig]] * GlobalC::ucell.tpiba * GlobalC::pw.get_G_cartesian_projection(ig, 0) * cpm.real();
                        forcescc(iat, 1) += fact * rhocgnt[GlobalC::pw.ig2ngg[ig]] * GlobalC::ucell.tpiba * GlobalC::pw.get_G_cartesian_projection(ig, 1) * cpm.real();
                        forcescc(iat, 2) += fact * rhocgnt[GlobalC::pw.ig2ngg[ig]] * GlobalC::ucell.tpiba * GlobalC::pw.get_G_cartesian_projection(ig, 2) * cpm.real();
                    }
					//std::cout << " forcescc = " << forcescc(iat,0) << " " << forcescc(iat,1) << " " << forcescc(iat,2) << std::endl;
                }
                iat++;
            }
        }
    }
    
	Parallel_Reduce::reduce_double_pool(forcescc.c, forcescc.nr * forcescc.nc);

	delete[] psic; //mohan fix bug 2012-03-22
	delete[] aux; //mohan fix bug 2012-03-22
	delete[] rhocgnt;  //mohan fix bug 2012-03-22

    return;
}




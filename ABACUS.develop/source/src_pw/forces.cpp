#include "forces.h"
#include "global.h"
#include "vdwd2.h"
#include "vdwd3.h"				  
#include "symmetry.h"
// new
#include "H_XC_pw.h"
#include "../src_global/math_integral.h"

double Forces::output_acc = 1.0e-8; // (Ryd/angstrom).	

Forces::Forces()
{
}

Forces::~Forces() {}

#include "efield.h"
void Forces::init(matrix& force)
{
	TITLE("Forces", "init");
	this->nat = ucell.nat;
	force.create(nat, 3);
	
	matrix forcelc(nat, 3);
	matrix forceion(nat, 3);
	matrix forcecc(nat, 3);
	matrix forcenl(nat, 3);
	matrix forcescc(nat, 3);
    this->cal_force_loc(forcelc);
    this->cal_force_ew(forceion);
    this->cal_force_nl(forcenl);
	this->cal_force_cc(forcecc);
	this->cal_force_scc(forcescc);

	matrix stress_vdw_pw;//.create(3,3);
    matrix force_vdw;
    force_vdw.create(nat, 3);
	if(vdwd2_para.flag_vdwd2)													//Peize Lin add 2014.04.03, update 2021.03.09
	{
        Vdwd2 vdwd2(ucell,vdwd2_para);
		vdwd2.cal_force();
		for(int iat=0; iat<ucell.nat; ++iat)
		{
			force_vdw(iat,0) = vdwd2.get_force()[iat].x;
			force_vdw(iat,1) = vdwd2.get_force()[iat].y;
			force_vdw(iat,2) = vdwd2.get_force()[iat].z;
		}
		if(TEST_FORCE)
		{
			Forces::print("VDW      FORCE (Ry/Bohr)", force_vdw);
		}
	}
	else if(vdwd3_para.flag_vdwd3)													//jiyy add 2019-05-18, update 2021-05-02
	{
        Vdwd3 vdwd3(ucell,vdwd3_para);
		vdwd3.cal_force();
		for(int iat=0; iat<ucell.nat; ++iat)
		{
			force_vdw(iat,0) = vdwd3.get_force()[iat].x;
			force_vdw(iat,1) = vdwd3.get_force()[iat].y;
			force_vdw(iat,2) = vdwd3.get_force()[iat].z;
		}
		if(TEST_FORCE)
		{
			Forces::print("VDW      FORCE (Ry/Bohr)", force_vdw);
		}
	}
    //impose total force = 0
    int iat = 0;

	matrix force_e;
	if(EFIELD)
	{
		force_e.create(ucell.nat, 3);
		Efield::compute_force(force_e);
	}
	
	for (int ipol = 0; ipol < 3; ipol++)
	{
		double sum = 0.0;
		iat = 0;

		for (int it = 0;it < ucell.ntype;it++)
		{
			for (int ia = 0;ia < ucell.atoms[it].na;ia++)
			{
				force(iat, ipol) =
					forcelc(iat, ipol)
					+ forceion(iat, ipol)
					+ forcenl(iat, ipol)
					+ forcecc(iat, ipol)
					+ forcescc(iat, ipol);

				if(vdwd2_para.flag_vdwd2 || vdwd3_para.flag_vdwd3)		//linpz and jiyy added vdw force, modified by zhengdy
				{
                    force(iat, ipol) += force_vdw(iat, ipol);
                }																										   
					
				if(EFIELD)
				{
					force(iat,ipol) = force(iat, ipol) + force_e(iat, ipol);
				}

				sum += force(iat, ipol);

				iat++;
			}
		}

		double compen = sum / ucell.nat;
		for(int iat=0; iat<ucell.nat; ++iat)
		{
			force(iat, ipol) = force(iat, ipol) - compen;
		}	
	}
	
	if(Symmetry::symm_flag)
	{
		double *pos;
		double d1,d2,d3;
		pos = new double[ucell.nat*3];
		ZEROS(pos, ucell.nat*3);
		int iat = 0;
		for(int it = 0;it < ucell.ntype;it++)
		{
			//Atom* atom = &ucell.atoms[it];
			for(int ia =0;ia< ucell.atoms[it].na;ia++)
			{
				pos[3*iat  ] = ucell.atoms[it].taud[ia].x ;
				pos[3*iat+1] = ucell.atoms[it].taud[ia].y ;
				pos[3*iat+2] = ucell.atoms[it].taud[ia].z;
				for(int k=0; k<3; ++k)
				{
					symm.check_translation( pos[iat*3+k], -floor(pos[iat*3+k]));
					symm.check_boundary( pos[iat*3+k] );
				}
				iat++;				
			}
		}
		
		for(int iat=0; iat<ucell.nat; iat++)
		{
			Mathzone::Cartesian_to_Direct(force(iat,0),force(iat,1),force(iat,2),
                                        ucell.a1.x, ucell.a1.y, ucell.a1.z,
                                        ucell.a2.x, ucell.a2.y, ucell.a2.z,
                                        ucell.a3.x, ucell.a3.y, ucell.a3.z,
                                        d1,d2,d3);
			
			force(iat,0) = d1;force(iat,1) = d2;force(iat,2) = d3;
		}
		symm.force_symmetry(force , pos);
		for(int iat=0; iat<ucell.nat; iat++)
		{
			Mathzone::Direct_to_Cartesian(force(iat,0),force(iat,1),force(iat,2),
                                        ucell.a1.x, ucell.a1.y, ucell.a1.z,
                                        ucell.a2.x, ucell.a2.y, ucell.a2.z,
                                        ucell.a3.x, ucell.a3.y, ucell.a3.z,
                                        d1,d2,d3);
			force(iat,0) = d1;force(iat,1) = d2;force(iat,2) = d3;
		}
		// cout << "nrotk =" << symm.nrotk << endl;
		delete[] pos;
		
	}

 	ofs_running << setiosflags(ios::fixed) << setprecision(6) << endl;
	if(TEST_FORCE)
	{
		Forces::print("LOCAL    FORCE (Ry/Bohr)", forcelc);
		Forces::print("NONLOCAL FORCE (Ry/Bohr)", forcenl);
		Forces::print("NLCC     FORCE (Ry/Bohr)", forcecc);
		Forces::print("ION      FORCE (Ry/Bohr)", forceion);
		Forces::print("SCC      FORCE (Ry/Bohr)", forcescc);
		if(EFIELD) Forces::print("EFIELD   FORCE (Ry/Bohr)", force_e);
	}
	
/*
	Forces::print("   TOTAL-FORCE (Ry/Bohr)", force);
	
	if(INPUT.force_set)                                                   // pengfei 2016-12-20
	{
		ofstream ofs("FORCE.dat");
		if(!ofs)
		{
			cout << "open FORCE.dat error !" <<endl;
		}
		for(int iat=0; iat<ucell.nat; iat++)
		{
			ofs << "   " << force(iat,0)*Ry_to_eV / 0.529177 
				<< "   " << force(iat,1)*Ry_to_eV / 0.529177 
				<< "   " << force(iat,2)*Ry_to_eV / 0.529177 << endl;
		}
		ofs.close();
	}
*/
		
	// output force in unit eV/Angstrom
	ofs_running << endl;

	if(TEST_FORCE)
	{
		Forces::print("LOCAL    FORCE (eV/Angstrom)", forcelc,0);
		Forces::print("NONLOCAL FORCE (eV/Angstrom)", forcenl,0);
		Forces::print("NLCC     FORCE (eV/Angstrom)", forcecc,0);
		Forces::print("ION      FORCE (eV/Angstrom)", forceion,0);
		Forces::print("SCC      FORCE (eV/Angstrom)", forcescc,0);
		if(EFIELD) Forces::print("EFIELD   FORCE (eV/Angstrom)", force_e,0);
	}
	Forces::print("   TOTAL-FORCE (eV/Angstrom)", force,0);

    return;
}

void Forces::print_to_files(ofstream &ofs, const string &name, const matrix &f)
{
    int iat = 0;
    ofs << " " << name;
    ofs << setprecision(8);
	//ofs << setiosflags(ios::showpos);
   
	double fac = Ry_to_eV / 0.529177;// (eV/A)

	if(TEST_FORCE)
	{
		cout << setiosflags(ios::showpos);
		cout << " " << name;
		cout << setprecision(8);
	}

    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
            ofs << " " << setw(5) << it
            << setw(8) << ia+1
            << setw(20) << f(iat, 0)*fac
            << setw(20) << f(iat, 1)*fac
            << setw(20) << f(iat, 2)*fac << endl;
			
			if(TEST_FORCE)
			{
            	cout << " " << setw(5) << it
            	<< setw(8) << ia+1
            	<< setw(20) << f(iat, 0)*fac
            	<< setw(20) << f(iat, 1)*fac
            	<< setw(20) << f(iat, 2)*fac << endl;
			}
            iat++;
        }
    }

	ofs_running << resetiosflags(ios::showpos);
	cout << resetiosflags(ios::showpos);
    return;
}



void Forces::print(const string &name, const matrix &f, bool ry)
{
	NEW_PART(name);

	ofs_running << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
	ofs_running << setiosflags(ios::showpos);

	const double fac = Ry_to_eV / 0.529177;
	
	if(TEST_FORCE)
	{
		cout << " --------------- " << name << " ---------------" << endl;
		cout << " " << setw(8) << "atom" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl;
		cout << setiosflags(ios::showpos);
		cout << setprecision(6);
		cout << setiosflags(ios::fixed);
	}

    int iat = 0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
			stringstream ss;
			ss << ucell.atoms[it].label << ia+1;

			if(ry) // output Rydberg Unit
			{
				ofs_running << " " << setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) ofs_running << setw(15) << f(iat,0);
				else ofs_running << setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) ofs_running << setw(15) << f(iat,1);
				else ofs_running << setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) ofs_running << setw(15) << f(iat,2);
				else ofs_running << setw(15) << "0";
				ofs_running << endl;
			}
			else
			{
				ofs_running << " " << setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) ofs_running << setw(15) << f(iat,0)*fac;
				else ofs_running << setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) ofs_running << setw(15) << f(iat,1)*fac;
				else ofs_running << setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) ofs_running << setw(15) << f(iat,2)*fac;
				else ofs_running << setw(15) << "0";
				ofs_running << endl;
			}

			if(TEST_FORCE && ry)
			{
				cout << " " << setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) cout << setw(15) << f(iat,0);
				else cout << setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) cout << setw(15) << f(iat,1);
				else cout << setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) cout << setw(15) << f(iat,2);
				else cout << setw(15) << "0";
				cout << endl;
			}
			else if (TEST_FORCE)
			{
				cout << " " << setw(8) << ss.str();
				if( abs(f(iat,0)) > Forces::output_acc) cout << setw(15) << f(iat,0)*fac;
				else cout << setw(15) << "0";
				if( abs(f(iat,1)) > Forces::output_acc) cout << setw(15) << f(iat,1)*fac;
				else cout << setw(15) << "0";
				if( abs(f(iat,2)) > Forces::output_acc) cout << setw(15) << f(iat,2)*fac;
				else cout << setw(15) << "0";
				cout << endl;
			}	
				
            iat++;
        }
    }

	ofs_running << resetiosflags(ios::showpos);
	cout << resetiosflags(ios::showpos);
    return;
}


void Forces::cal_force_loc(matrix& forcelc)
{
	timer::tick("Forces","cal_force_loc");

    complex<double> *aux = new complex<double>[pw.nrxx];
    ZEROS(aux, pw.nrxx);

    // now, in all pools , the charge are the same,
    // so, the force calculated by each pool is equal.
    
	for(int is=0; is<NSPIN; is++)
	{
		for (int ir=0; ir<pw.nrxx; ir++)
		{
        	aux[ir] += complex<double>( CHR.rho[is][ir], 0.0 );
		}
	}

	// to G space.
    pw.FFT_chg.FFT3D(aux, -1);

    int gstart_here = pw.gstart;
    if (pw.ggs[0] != 0) gstart_here = 0;

//  ofs_running << "\n ggs = " << pw.ggs[0];
//  ofs_running << "\n gstart_here = " << gstart_here;
    int iat = 0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
            for (int ig = gstart_here; ig < pw.ngmc; ig++)
            {
                const double phase = TWO_PI * (pw.gcar[ig] * ucell.atoms[it].tau[ia]);
                const double factor = ppcell.vloc(it, pw.ig2ngg[ig]) *
									  ( cos(phase) * aux[ pw.ig2fftc[ig] ].imag()
                                      + sin(phase) * aux[ pw.ig2fftc[ig] ].real()); 
                forcelc(iat, 0) += pw.gcar[ig].x * factor;
                forcelc(iat, 1) += pw.gcar[ig].y * factor;
                forcelc(iat, 2) += pw.gcar[ig].z * factor;
            }
            for (int ipol = 0;ipol < 3;ipol++)
            {
                forcelc(iat, ipol) *= (ucell.tpiba * ucell.omega);
            }
            ++iat;
        }
    }
    //this->print(ofs_running, "local forces", forcelc);
    Parallel_Reduce::reduce_double_pool(forcelc.c, forcelc.nr * forcelc.nc);
    delete[] aux;
	timer::tick("Forces","cal_force_loc");
    return;
}

#include "H_Ewald_pw.h"
void Forces::cal_force_ew(matrix& forceion)
{
	timer::tick("Forces","cal_force_ew");

    double fact = 2.0;
    complex<double> *aux = new complex<double> [pw.ngmc];
    ZEROS(aux, pw.ngmc);

    int gstart = pw.gstart;

    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ig = gstart; ig < pw.ngmc; ig++)
        {
            aux[ig] += static_cast<double>(ucell.atoms[it].zv) * conj(pw.strucFac(it, ig));
        }
    }

	// calculate total ionic charge
    double charge = 0.0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        charge += ucell.atoms[it].na * ucell.atoms[it].zv;//mohan modify 2007-11-7
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
            WARNING_QUIT("ewald","Can't find optimal alpha.");
        }
        upperbound = 2.0 * charge * charge * sqrt(2.0 * alpha / TWO_PI) *
                     erfc(sqrt(ucell.tpiba2 * pw.ggchg / 4.0 / alpha));
    }
    while (upperbound > 1.0e-6);
//	cout << " en.alpha = " << alpha << endl;
//	cout << " upperbound = " << upperbound << endl;
	


    for (int ig = gstart; ig < pw.ngmc; ig++)
    {
        if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            aux[ig] *= exp(-1.0 * pw.gg[ig] * ucell.tpiba2 / alpha / 4.0) / (pw.gg[ig] * ucell.tpiba2);
        }
    }

    int iat = 0;
    for (int it = 0;it < ucell.ntype;it++)
    {
        for (int ia = 0;ia < ucell.atoms[it].na;ia++)
        {
            for (int ig = gstart; ig < pw.ngmc; ig++)
            {
                const double arg = TWO_PI * (pw.gcar[ig] * ucell.atoms[it].tau[ia]);
                double sumnb =  -cos(arg) * aux[ig].imag() + sin(arg) * aux[ig].real();
                forceion(iat, 0) += pw.gcar[ig].x * sumnb;
                forceion(iat, 1) += pw.gcar[ig].y * sumnb;
                forceion(iat, 2) += pw.gcar[ig].z * sumnb;
            }
            for (int ipol = 0;ipol < 3;ipol++)
            {
                forceion(iat, ipol) *= ucell.atoms[it].zv * e2 * ucell.tpiba * TWO_PI / ucell.omega * fact;
            }

	//		cout << " atom" << iat << endl;
	//		cout << setw(15) << forceion(iat, 0) << setw(15) << forceion(iat,1) << setw(15) << forceion(iat,2) << endl; 
            iat++;
        }
    }
    delete [] aux;


	// means that the processor contains G=0 term.
    if (gstart == 1)
    {
        double rmax = 5.0 / (sqrt(alpha) * ucell.lat0);
        int nrm = 0;
		
        //output of rgen: the number of vectors in the sphere
        const int mxr = 50;
        // the maximum number of R vectors included in r
        Vector3<double> *r  = new Vector3<double>[mxr];
        double *r2 = new double[mxr];
		ZEROS(r2, mxr);
        int *irr = new int[mxr];
		ZEROS(irr, mxr);
        // the square modulus of R_j-tau_s-tau_s'

		int iat1 = 0;
        for (int T1 = 0; T1 < ucell.ntype; T1++)
        {
			Atom* atom1 = &ucell.atoms[T1]; 
            for (int I1 = 0; I1 < atom1->na; I1++)
            {
				int iat2 = 0; // mohan fix bug 2011-06-07
                for (int T2 = 0; T2 < ucell.ntype; T2++)
                {
                    for (int I2 = 0; I2 < ucell.atoms[T2].na; I2++)
                    {
                        if (iat1 != iat2)
                        {
                            Vector3<double> d_tau = ucell.atoms[T1].tau[I1] - ucell.atoms[T2].tau[I2];
                            H_Ewald_pw::rgen(d_tau, rmax, irr, ucell.latvec, ucell.G, r, r2, nrm);

                            for (int n = 0;n < nrm;n++)
                            {
								const double rr = sqrt(r2[n]) * ucell.lat0;

                                double factor = ucell.atoms[T1].zv * ucell.atoms[T2].zv * e2 / (rr * rr)
                                                * (erfc(sqrt(alpha) * rr) / rr
                                    + sqrt(8.0 * alpha / TWO_PI) * exp(-1.0 * alpha * rr * rr)) * ucell.lat0;

								forceion(iat1, 0) -= factor * r[n].x;
                                forceion(iat1, 1) -= factor * r[n].y;
                                forceion(iat1, 2) -= factor * r[n].z;

//								cout << " r.z=" << r[n].z << " r2=" << r2[n] << endl;
						//		cout << " " << iat1 << " " << iat2 << " n=" << n
						//		 << " rn.z=" << r[n].z 
						//		 << " r2=" << r2[n] << " rr=" << rr << " fac=" << factor << " force=" << forceion(iat1,2) 
						//		 << " new_part=" << factor*r[n].z <<  endl;
                            }
                        }

                        ++iat2;
                    }
                }//atom b

//				cout << " atom" << iat1 << endl;
//				cout << setw(15) << forceion(iat1, 0) << setw(15) << forceion(iat1,1) << setw(15) << forceion(iat1,2) << endl; 

                ++iat1;
            }
        }//atom a
    }

    Parallel_Reduce::reduce_double_pool(forceion.c, forceion.nr * forceion.nc);

    //this->print(ofs_running, "ewald forces", forceion);

	timer::tick("Forces","cal_force_ew");

    return;
}

void Forces::cal_force_cc(matrix& forcecc)
{
	// recalculate the exchange-correlation potential.
    matrix vxc(NSPIN, pw.nrxx);
    H_XC_pw::v_xc(pw.nrxx, pw.ncxyz, ucell.omega, CHR.rho, CHR.rho_core, vxc);

    complex<double> * psiv = new complex<double> [pw.nrxx];
    ZEROS(psiv, pw.nrxx);
    if (NSPIN == 1 || NSPIN == 4)
    {
        for (int ir = 0;ir < pw.nrxx;ir++)
        {
            psiv[ir] = complex<double>(vxc(0, ir),  0.0);
        }
    }
    else
    {
        for (int ir = 0;ir < pw.nrxx;ir++)
        {
            psiv[ir] = 0.5 * (vxc(0 ,ir) + vxc(1, ir));
        }
    }

	// to G space
    pw.FFT_chg.FFT3D(psiv, -1);

    //psiv contains now Vxc(G)
    double * rhocg = new double [pw.nggm];
    ZEROS(rhocg, pw.nggm);

    for (int T1 = 0;T1 < ucell.ntype;T1++)
    {
        if (ucell.atoms[T1].nlcc)
        {
            //call drhoc
            CHR.non_linear_core_correction(
                ppcell.numeric,
                ucell.atoms[T1].msh,
                ucell.atoms[T1].r,
                ucell.atoms[T1].rab,
                ucell.atoms[T1].rho_atc,
                rhocg);
            int iat = 0;


			complex<double> ipol0, ipol1, ipol2;
            for (int T2 = 0;T2 < ucell.ntype;T2++)
            {
                for (int I2 = 0;I2 < ucell.atoms[T2].na;I2++)
                {
                    if (T2 == T1)
                    {
                        for (int ig = pw.gstart; ig < pw.ngmc; ig++)
                        {
                            const double arg = TWO_PI * (pw.gcar[ig].x * ucell.atoms[T2].tau[I2].x
                                                      + pw.gcar[ig].y * ucell.atoms[T2].tau[I2].y
                                                      + pw.gcar[ig].z * ucell.atoms[T2].tau[I2].z);

                            ipol0 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].x * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg))  ;
                            forcecc(iat, 0) +=  ipol0.real();

                            ipol1 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].y * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg)) ;
                            forcecc(iat, 1) += ipol1.real();

                            ipol2 = ucell.tpiba * ucell.omega * rhocg[pw.ig2ngg[ig]]
                                                    * pw.gcar[ig].z * conj(psiv[pw.ig2fftc[ig]]) 
													* complex<double>(sin(arg), cos(arg)) ;

							forcecc(iat, 2) += ipol2.real();
                        }
                        ++iat;
                    }
                }
            }
        }
    }
    delete [] rhocg;
	delete [] psiv; // mohan fix bug 2012-03-22
    Parallel_Reduce::reduce_double_all(forcecc.c, forcecc.nr * forcecc.nc);
	return;
}

void Forces::cal_force_nl(matrix& forcenl)
{
	TITLE("Forces","cal_force_nl");
	timer::tick("Forces","cal_force_nl");

    const int nkb = ppcell.nkb;
	if(nkb == 0) return; // mohan add 2010-07-25
	
	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	ComplexArray dbecp( nkb, NBANDS, 3);
    ComplexMatrix becp( nkb, NBANDS);
    
	
	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ComplexMatrix vkb1( nkb, wf.npwx );

    for (int ik = 0;ik < kv.nks;ik++)
    {
        if (NSPIN==2) CURRENT_SPIN = kv.isk[ik];
        wf.npw = kv.ngk[ik];
        // generate vkb
        if (ppcell.nkb > 0)
        {
            ppcell.getvnl(ik);
        }

        // get becp according to wave functions and vkb
        // important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        becp.zero_out();
        for (int ib=0; ib<NBANDS; ib++)
        {
            for (int i=0;i<nkb;i++)
            {
                for (int ig=0; ig<wf.npw; ig++)
                {
                    becp(i,ib) += wf.evc[ik](ib,ig)* conj( ppcell.vkb(i,ig) );
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
					for (int ig=0; ig<wf.npw; ig++)
						vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * pw.gcar[ wf.igk(ik, ig) ].x;
				}
				if (ipol==1)
				{
					for (int ig=0; ig<wf.npw; ig++)
						vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * pw.gcar[ wf.igk(ik, ig) ].y;
				}
				if (ipol==2)
				{
					for (int ig=0; ig<wf.npw; ig++)
						vkb1(i, ig) = ppcell.vkb(i, ig) * NEG_IMAG_UNIT * pw.gcar[ wf.igk(ik, ig) ].z;
				}
			}
            for (int ib=0; ib<NBANDS; ib++)
            {
                for (int i=0; i<nkb; i++)
                {
                    for (int ig=0; ig<wf.npw; ig++)
                    {
                        dbecp(i,ib, ipol) += conj( vkb1(i,ig) ) * wf.evc[ik](ib,ig) ;
                    }
                }
            }
        }// end ipol

//		don't need to reduce here, keep dbecp different in each processor,
//		and at last sum up all the forces.
//		Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//		double *cf = new double[ucell.nat*3];
//		ZEROS(cf, ucell.nat);
		for (int ib=0; ib<NBANDS; ib++)
		{
			double fac = wf.wg(ik, ib) * 2.0 * ucell.tpiba;
        	int iat = 0;
        	int sum = 0;
			for (int it=0; it<ucell.ntype; it++)
			{
				const int Nprojs = ucell.atoms[it].nh;
				for (int ia=0; ia<ucell.atoms[it].na; ia++)
				{
					for (int ip=0; ip<Nprojs; ip++)
					{
						double ps = ppcell.deeq(CURRENT_SPIN, iat, ip, ip) ;
						const int inkb = sum + ip; 
						//out<<"\n ps = "<<ps;

						for (int ipol=0; ipol<3; ipol++)
						{
							const double dbb = ( conj( dbecp( inkb, ib, ipol) ) * becp( inkb, ib) ).real();
							forcenl(iat, ipol) = forcenl(iat, ipol) - ps * fac * dbb;
							//cf[iat*3+ipol] += ps * fac * dbb;
						}
					}

					//if ( ucell.atoms[it].nbeta > ucell.atoms[it].lmax+1 )    //{zws add 20160110
					//{
					//cout << " \n multi-projector force calculation ... " << endl;
					for (int ip=0; ip<Nprojs; ip++)
					{
						const int inkb = sum + ip;
						//for (int ip2=0; ip2<Nprojs; ip2++)
						for (int ip2=ip+1; ip2<Nprojs; ip2++)
						{
						//if ( ip != ip2 )
						//{
							const int jnkb = sum + ip2;
							double ps = ppcell.deeq(CURRENT_SPIN, iat, ip2, ip) ;

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
//  this->print(ofs_running, "nonlocal forces", forcenl);
	timer::tick("Forces","cal_force_nl");
    return;
}

void Forces::cal_force_scc(matrix& forcescc)
{
    complex<double>* psic = new complex<double> [pw.nrxx];
    ZEROS(psic, pw.nrxx);

    if (NSPIN == 1 || NSPIN == 4)
    {
        for (int i = 0;i < pw.nrxx;i++)
        {
            psic[i] = pot.vnew(0,i);
        }
    }
    else
    {
        int isup = 0;
        int isdw = 1;
        for (int i = 0;i < pw.nrxx;i++)
        {
            psic[i] = (pot.vnew(isup, i) + pot.vnew(isdw, i)) * 0.5;
        }
    }

    int ndm = 0;

    for (int it = 0;it < ucell.ntype;it++)
    {
        if (ndm < ucell.atoms[it].msh)
        {
            ndm = ucell.atoms[it].msh;
        }
    }

    //work space
    double* aux = new double[ndm];
    ZEROS(aux, ndm);

    double* rhocgnt = new double[pw.nggm];
    ZEROS(rhocgnt, pw.nggm);

    pw.FFT_chg.FFT3D(psic, -1);

    double fact = 2.0;
    for (int nt = 0;nt < ucell.ntype;nt++)
    {
//		Here we compute the G.ne.0 term
        const int mesh = ucell.atoms[nt].msh;
        for (int ig = pw.gstart;ig < pw.nggm;ig++)
        {
            const double gx = sqrt(pw.ggs[ig]) * ucell.tpiba;
            for (int ir = 0;ir < mesh;ir++)
            {
                if (ucell.atoms[nt].r[ir] < 1.0e-8)
                {
                    aux[ir] = ucell.atoms[nt].rho_at[ir];
                }
                else
                {
                    const double gxx = gx * ucell.atoms[nt].r[ir];
                    aux[ir] = ucell.atoms[nt].rho_at[ir] * sin(gxx) / gxx;
                }
            }
            Integral::Simpson_Integral(mesh , aux, ucell.atoms[nt].rab , rhocgnt [ig]);
        }

        int iat = 0;
        for (int it = 0;it < ucell.ntype;it++)
        {
            for (int ia = 0;ia < ucell.atoms[it].na;ia++)
            {
                if (nt == it)
                {
                    for (int ig = pw.gstart;ig < pw.ngmc;ig++)
                    {
                        const double arg = TWO_PI * (pw.gcar[ig].x * ucell.atoms[it].tau[ia].x
                                                  + pw.gcar[ig].y * ucell.atoms[it].tau[ia].y
                                                  + pw.gcar[ig].z * ucell.atoms[it].tau[ia].z);
                        
						const complex<double> cpm = complex<double>(sin(arg), cos(arg)) * conj(psic[pw.ig2fftc[ig] ]);
						
                        forcescc(iat, 0) += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].x * cpm.real();
                        forcescc(iat, 1) += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].y * cpm.real();
                        forcescc(iat, 2) += fact * rhocgnt [pw.ig2ngg[ig] ] * ucell.tpiba * pw.gcar[ig].z * cpm.real();
                    }
					//cout << " forcescc = " << forcescc(iat,0) << " " << forcescc(iat,1) << " " << forcescc(iat,2) << endl;
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




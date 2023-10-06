#include "sto_forces.h"

#include "module_base/mathzone.h"
#include "module_cell/module_symmetry/symmetry.h"
#include "module_elecstate/elecstate.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

// new
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_base/math_integral.h"
#include "module_base/parallel_reduce.h"
#include "module_base/timer.h"

void Sto_Forces::cal_stoforce(ModuleBase::matrix& force,
                              const elecstate::ElecState& elec,
                              ModulePW::PW_Basis* rho_basis,
                              ModuleSymmetry::Symmetry* p_symm,
                              Structure_Factor* p_sf,
                              K_Vectors* pkv,
                              ModulePW::PW_Basis_K* wfc_basis,
                              const psi::Psi<std::complex<double>>* psi_in,
                              Stochastic_WF& stowf)
{
    ModuleBase::timer::tick("Sto_Force","cal_force");
	ModuleBase::TITLE("Sto_Forces", "init");
	this->nat =  GlobalC::ucell.nat;
    const ModuleBase::matrix& wg = elec.wg;
    const Charge* chr = elec.charge;
    force.create(nat, 3);

    ModuleBase::matrix forcelc(nat, 3);
    ModuleBase::matrix forceion(nat, 3);
    ModuleBase::matrix forcecc(nat, 3);
    ModuleBase::matrix forcenl(nat, 3);
    ModuleBase::matrix forcescc(nat, 3);
    this->cal_force_loc(forcelc, rho_basis, chr);
    this->cal_force_ew(forceion, rho_basis, p_sf);
    this->cal_sto_force_nl(forcenl, wg, pkv, wfc_basis, psi_in, stowf);
    this->cal_force_cc(forcecc, rho_basis, chr);
    this->cal_force_scc(forcescc, rho_basis, elec.vnew, elec.vnew_exist);

    //impose total force = 0
    int iat = 0;

	ModuleBase::matrix force_e;
	if(GlobalV::EFIELD_FLAG)
	{
		force_e.create( GlobalC::ucell.nat, 3);
		elecstate::Efield::compute_force(GlobalC::ucell, force_e);
	}

    ModuleBase::matrix force_gate;
    if(GlobalV::GATE_FLAG)
    {
        force_gate.create( GlobalC::ucell.nat, 3);
        elecstate::Gatefield::compute_force(GlobalC::ucell, force_gate);
    }

	for (int ipol = 0; ipol < 3; ipol++)
	{
		double sum = 0.0;
		iat = 0;

		for (int it = 0;it <  GlobalC::ucell.ntype;it++)
		{
			for (int ia = 0;ia <  GlobalC::ucell.atoms[it].na;ia++)
			{
				force(iat, ipol) =
					forcelc(iat, ipol)
					+ forceion(iat, ipol)
					+ forcenl(iat, ipol)
					+ forcecc(iat, ipol)
					+ forcescc(iat, ipol);																									   
					
				if(GlobalV::EFIELD_FLAG)
				{
					force(iat,ipol) = force(iat, ipol) + force_e(iat, ipol);
				}

                if(GlobalV::GATE_FLAG)
                {
                    force(iat,ipol) = force(iat, ipol) + force_gate(iat, ipol);
                }

				sum += force(iat, ipol);

				iat++;
			}
		}

        if(!(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG))
        {
            double compen = sum / GlobalC::ucell.nat;
            for (int iat = 0; iat < GlobalC::ucell.nat; ++iat)
            {
                force(iat, ipol) = force(iat, ipol) - compen;
            }
        }
	}

    if(GlobalV::GATE_FLAG || GlobalV::EFIELD_FLAG)
    {
        GlobalV::ofs_running << "Atomic forces are not shifted if gate_flag or efield_flag == true!" << std::endl;
    }
	
	if(ModuleSymmetry::Symmetry::symm_flag == 1)
	{
		double *pos;
		double d1,d2,d3;
		pos = new double[GlobalC::ucell.nat*3];
		ModuleBase::GlobalFunc::ZEROS(pos, GlobalC::ucell.nat*3);
		int iat = 0;
		for(int it = 0;it < GlobalC::ucell.ntype;it++)
		{
			//Atom* atom = &ucell.atoms[it];
			for(int ia =0;ia< GlobalC::ucell.atoms[it].na;ia++)
			{
				pos[3*iat  ] = GlobalC::ucell.atoms[it].taud[ia].x ;
				pos[3*iat+1] = GlobalC::ucell.atoms[it].taud[ia].y ;
				pos[3*iat+2] = GlobalC::ucell.atoms[it].taud[ia].z;
				for(int k=0; k<3; ++k)
				{
                    p_symm->check_translation(pos[iat * 3 + k], -floor(pos[iat * 3 + k]));
                    p_symm->check_boundary(pos[iat * 3 + k]);
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
        p_symm->force_symmetry(force, pos, GlobalC::ucell);
        for (int iat = 0; iat < GlobalC::ucell.nat; iat++)
        {
			ModuleBase::Mathzone::Direct_to_Cartesian(force(iat,0),force(iat,1),force(iat,2),
                                        GlobalC::ucell.a1.x, GlobalC::ucell.a1.y, GlobalC::ucell.a1.z,
                                        GlobalC::ucell.a2.x, GlobalC::ucell.a2.y, GlobalC::ucell.a2.z,
                                        GlobalC::ucell.a3.x, GlobalC::ucell.a3.y, GlobalC::ucell.a3.z,
                                        d1,d2,d3);
			force(iat,0) = d1;force(iat,1) = d2;force(iat,2) = d3;
		}
        // cout << "nrotk =" << p_symm->nrotk << endl;
        delete[] pos;
    }

 	GlobalV::ofs_running << setiosflags(std::ios::fixed) << std::setprecision(6) << std::endl;
	if(GlobalV::TEST_FORCE)
	{
		Sto_Forces::print("LOCAL    FORCE (Ry/Bohr)", forcelc);
		Sto_Forces::print("NONLOCAL FORCE (Ry/Bohr)", forcenl);
		Sto_Forces::print("NLCC     FORCE (Ry/Bohr)", forcecc);
		Sto_Forces::print("ION      FORCE (Ry/Bohr)", forceion);
		Sto_Forces::print("SCC      FORCE (Ry/Bohr)", forcescc);
		if(GlobalV::EFIELD_FLAG) Sto_Forces::print("EFIELD   FORCE (Ry/Bohr)", force_e);
        if(GlobalV::GATE_FLAG) Sto_Forces::print("GATEFIELD   FORCE (Ry/Bohr)", force_gate);
	}
	
		
	// output force in unit eV/Angstrom
	GlobalV::ofs_running << std::endl;
    
	if(GlobalV::TEST_FORCE)
	{
		Sto_Forces::print("LOCAL    FORCE (eV/Angstrom)", forcelc,0);
		Sto_Forces::print("NONLOCAL FORCE (eV/Angstrom)", forcenl,0);
		Sto_Forces::print("NLCC     FORCE (eV/Angstrom)", forcecc,0);
		Sto_Forces::print("ION      FORCE (eV/Angstrom)", forceion,0);
		Sto_Forces::print("SCC      FORCE (eV/Angstrom)", forcescc,0);
		if(GlobalV::EFIELD_FLAG) Sto_Forces::print("EFIELD   FORCE (eV/Angstrom)", force_e,0);
        if(GlobalV::GATE_FLAG) Sto_Forces::print("GATEFIELD   FORCE (eV/Angstrom)", force_gate,0);
	}
	Sto_Forces::print("   TOTAL-FORCE (eV/Angstrom)", force,0);
	ModuleBase::timer::tick("Sto_Force","cal_force");
    return;
}

void Sto_Forces::cal_sto_force_nl(ModuleBase::matrix& forcenl,
                                  const ModuleBase::matrix& wg,
                                  K_Vectors* p_kv,
                                  ModulePW::PW_Basis_K* wfc_basis,
                                  const psi::Psi<std::complex<double>>* psi_in,
                                  Stochastic_WF& stowf)
{
	ModuleBase::TITLE("Sto_Forces","cal_force_nl");
	ModuleBase::timer::tick("Sto_Forces","cal_force_nl");

    const int nkb = GlobalC::ppcell.nkb;
    int* nchip = stowf.nchip;
	if(nkb == 0) return; // mohan add 2010-07-25
	
	const int npwx = wfc_basis->npwk_max;
	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ModuleBase::ComplexMatrix vkb1( nkb, npwx );
	int nksbands = psi_in->get_nbands();
	if(GlobalV::MY_STOGROUP != 0) nksbands = 0;
    

    for (int ik = 0;ik < wfc_basis->nks;ik++)
    {
		const int nstobands = nchip[ik];
		const int nbandstot = nstobands + nksbands;
		const int npw = wfc_basis->npwk[ik];

		// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
		ModuleBase::ComplexArray dbecp( 3, nbandstot, nkb);
    	ModuleBase::ComplexMatrix becp( nbandstot, nkb);
        if (GlobalV::NSPIN == 2)
            GlobalV::CURRENT_SPIN = p_kv->isk[ik];

        // generate vkb
        if (GlobalC::ppcell.nkb > 0)
        {
            GlobalC::ppcell.getvnl(ik, GlobalC::ppcell.vkb);
        }

        // get becp according to wave functions and vkb
        // important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        becp.zero_out();
		char transa = 'C';
        char transb = 'N';
		psi_in[0].fix_k(ik);
		//KS orbitals
		int npmks = GlobalV::NPOL * nksbands;
		zgemm_(&transa,&transb,&nkb,&npmks,&npw,&ModuleBase::ONE,
				GlobalC::ppcell.vkb.c,&npwx,
            	psi_in->get_pointer(),&npwx,
            	&ModuleBase::ZERO,becp.c,&nkb);
		//stochastic orbitals
		int npmsto = GlobalV::NPOL * nstobands;
		zgemm_(&transa,&transb,&nkb,&npmsto,&npw,&ModuleBase::ONE,
				GlobalC::ppcell.vkb.c,&npwx,
            	stowf.shchi[ik].c,&npwx,
            	&ModuleBase::ZERO,&becp(nksbands,0),&nkb);
        
        Parallel_Reduce::reduce_pool(becp.c, becp.size);

        //out.printcm_real("becp",becp,1.0e-4);
        // Calculate the derivative of beta,
        // |dbeta> =  -ig * |beta>
        dbecp.zero_out();
        for (int ipol = 0; ipol<3; ipol++)
        {
			for (int i = 0;i < nkb;i++)
			{
				std::complex<double>* pvkb1 = &vkb1(i,0);
                std::complex<double>* pvkb = &GlobalC::ppcell.vkb(i,0);
				if (ipol==0)
				{
					for (int ig=0; ig<npw; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * wfc_basis->getgcar(ik, ig)[0];
                }
				if (ipol==1)
				{
					for (int ig=0; ig<npw; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * wfc_basis->getgcar(ik,ig)[1];
                }
				if (ipol==2)
				{
					for (int ig=0; ig<npw; ig++)
                        pvkb1[ig] = pvkb[ig] * ModuleBase::NEG_IMAG_UNIT * wfc_basis->getgcar(ik,ig)[2];
                }
			}
            //KS orbitals
			zgemm_(&transa,&transb,&nkb,&npmks,&npw,&ModuleBase::ONE,
					vkb1.c,&npwx,
        	    	psi_in->get_pointer(),&npwx,
        	    	&ModuleBase::ZERO,&dbecp(ipol, 0, 0),&nkb);
			//stochastic orbitals
			zgemm_(&transa,&transb,&nkb,&npmsto,&npw,&ModuleBase::ONE,
					vkb1.c,&npwx,
        	    	stowf.shchi[ik].c,&npwx,
        	    	&ModuleBase::ZERO,&dbecp(ipol, nksbands, 0),&nkb);
        }// end ipol

//		don't need to reduce here, keep dbecp different in each processor,
//		and at last sum up all the forces.
//		Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//		double *cf = new double[ucell.nat*3];
//		ZEROS(cf, ucell.nat);
		for (int ib=0; ib<nbandstot; ib++)
		{
			double fac;
			if(ib < nksbands)
				fac = wg(ik, ib) * 2.0 * GlobalC::ucell.tpiba;
			else
                fac = p_kv->wk[ik] * 2.0 * GlobalC::ucell.tpiba;
            int iat = 0;
        	int sum = 0;
			for (int it=0; it< GlobalC::ucell.ntype; it++)
			{
				const int Nprojs =  GlobalC::ucell.atoms[it].ncpp.nh;
				for (int ia=0; ia< GlobalC::ucell.atoms[it].na; ia++)
				{
					for (int ip=0; ip<Nprojs; ip++)
					{
						double ps =  GlobalC::ppcell.deeq( GlobalV::CURRENT_SPIN, iat, ip, ip) ;
						const int inkb = sum + ip; 
						//out<<"\n ps = "<<ps;

						for (int ipol=0; ipol<3; ipol++)
						{
							const double dbb = ( conj( dbecp( ipol, ib, inkb) ) * becp( ib, inkb) ).real();
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
							double ps =  GlobalC::ppcell.deeq( GlobalV::CURRENT_SPIN, iat, ip2, ip) ;
							for (int ipol=0; ipol<3; ipol++)
							{
								const double dbb = 2.0 * ( conj( dbecp( ipol, ib, inkb) ) * becp( ib, jnkb) ).real();
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
    Parallel_Reduce::reduce_all(forcenl.c, forcenl.nr* forcenl.nc);

	ModuleBase::timer::tick("Sto_Forces","cal_force_nl");
    return;
}



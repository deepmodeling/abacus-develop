#include "surchem.h"
#include "module_base/timer.h"

void force_cor_one(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, ModuleBase::matrix& forcesol)
{
   
   
    //delta phi multiply by the derivative of nuclear charge density with respect to the positions
    std::complex<double> *N = new std::complex<double>[rho_basis->npw];
    std::complex<double> *vloc_at = new std::complex<double>[rho_basis->npw];
    std::complex<double> *delta_phi_g = new complex<double>[rho_basis->npw];
    //ModuleBase::GlobalFunc::ZEROS(delta_phi_g, rho_basis->npw);

    rho_basis->real2recip(GlobalC::solvent_model.delta_phi, delta_phi_g);
    //GlobalC::UFFT.ToReciSpace(GlobalC::solvent_model.delta_phi, delta_phi_g,rho_basis);
    double Ael=0;double Ael1 = 0;
    //ModuleBase::GlobalFunc::ZEROS(vg, ngmc);
    int iat = 0;

    for (int it = 0;it < cell.ntype;it++)
    {
        for (int ia = 0;ia < cell.atoms[it].na ; ia++)
        {
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {   
                complex<double> phase = exp( ModuleBase::NEG_IMAG_UNIT *ModuleBase::TWO_PI * ( rho_basis->gcar[ig] * cell.atoms[it].tau[ia]));
                //vloc for each atom
                vloc_at[ig] = GlobalC::ppcell.vloc(it, rho_basis->ig2igg[ig]) * phase;
                if(rho_basis->ig_gge0 == ig)
                {
                    N[ig] = GlobalC::ucell.atoms[it].ncpp.zv / GlobalC::ucell.omega;
                }
                else
                {
                    const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                               (cell.tpiba2 * rho_basis->gg[ig]);

                    N[ig] = -vloc_at[ig] / fac;
                }
                
                //force for each atom
                forcesol(iat, 0) += rho_basis->gcar[ig][0] * imag(conj(delta_phi_g[ig]) * N[ig]);
                forcesol(iat, 1) += rho_basis->gcar[ig][1] * imag(conj(delta_phi_g[ig]) * N[ig]);
                forcesol(iat, 2) += rho_basis->gcar[ig][2] * imag(conj(delta_phi_g[ig]) * N[ig]);
            }
                
                forcesol(iat, 0) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
                forcesol(iat, 1) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
                forcesol(iat, 2) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
            //unit Ry/Bohr
                forcesol(iat, 0) *= 2 ;
                forcesol(iat, 1) *= 2 ;
                forcesol(iat, 2) *= 2 ;

                //cout<<"Force1"<<iat<<":"<<" "<<forcesol(iat, 0)<<" "<<forcesol(iat, 1)<<" "<<forcesol(iat, 2)<<endl;
               
            ++iat;
        }
    }

    delete[] vloc_at;
    delete[] N;
    delete[] delta_phi_g;

}

void force_cor_two(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, ModuleBase::matrix& forcesol)
{
   
    complex<double> *n_pseudo = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(n_pseudo,rho_basis->npw);

    //GlobalC::solvent_model.gauss_charge(cell, pwb, n_pseudo);
    
    double *Vcav_sum =  new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vcav_sum, rho_basis->nrxx);
    std::complex<double> *Vcav_g = new complex<double>[rho_basis->npw];
    std::complex<double> *Vel_g = new complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(Vcav_g, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(Vel_g, rho_basis->npw);
    for(int is=0; is<GlobalV::NSPIN; is++)
	{
		for (int ir=0; ir<rho_basis->nrxx; ir++)
		{
        	Vcav_sum[ir] += GlobalC::solvent_model.Vcav(is, ir);
		}
	}

    rho_basis->real2recip(Vcav_sum, Vcav_g);
    rho_basis->real2recip(GlobalC::solvent_model.epspot, Vel_g);

    int iat = 0;
    double Ael1 = 0;
    for (int it = 0;it < cell.ntype;it++)
    {
        double RCS = GlobalC::solvent_model.GetAtom.atom_RCS[cell.atoms[it].ncpp.psd];
        double sigma_rc_k = RCS / 2.5;
        for (int ia = 0;ia < cell.atoms[it].na;ia++)
        {
            //cell.atoms[0].tau[0].z = 3.302;
            //cout<<cell.atoms[it].tau[ia]<<endl;
             ModuleBase::GlobalFunc::ZEROS(n_pseudo, rho_basis->npw);
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {
                // G^2
                double gg = rho_basis->gg[ig];
                gg = gg * cell.tpiba2;
                complex<double> phase = exp( ModuleBase::NEG_IMAG_UNIT *ModuleBase::TWO_PI * ( rho_basis->gcar[ig] * cell.atoms[it].tau[ia]));

                n_pseudo[ig].real((GlobalC::solvent_model.GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv) * phase.real()
                             * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
                n_pseudo[ig].imag((GlobalC::solvent_model.GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv) * phase.imag()
                             * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
            }
            
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {   
                n_pseudo[ig] /= cell.omega;
            }
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {
                forcesol(iat, 0) -= rho_basis->gcar[ig][0] * imag(conj(Vcav_g[ig]+Vel_g[ig]) * n_pseudo[ig]);
                forcesol(iat, 1) -= rho_basis->gcar[ig][1] * imag(conj(Vcav_g[ig]+Vel_g[ig]) * n_pseudo[ig]);
                forcesol(iat, 2) -= rho_basis->gcar[ig][2] * imag(conj(Vcav_g[ig]+Vel_g[ig]) * n_pseudo[ig]);
            }

                forcesol(iat, 0) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
                forcesol(iat, 1) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
                forcesol(iat, 2) *= (GlobalC::ucell.tpiba * GlobalC::ucell.omega);
            //eV/Ang
                forcesol(iat, 0) *= 2 ;
                forcesol(iat, 1) *= 2 ;
                forcesol(iat, 2) *= 2 ;

                //cout<<"Force2"<<iat<<":"<<" "<<forcesol(iat, 0)<<" "<<forcesol(iat, 1)<<" "<<forcesol(iat, 2)<<endl;

            ++iat;
        }
    }
    
    delete[] n_pseudo;
    delete[] Vcav_sum;
    delete[] Vcav_g;
    delete[] Vel_g;

}

void surchem::cal_force_sol(const UnitCell& cell, const ModulePW::PW_Basis* rho_basis, ModuleBase::matrix& forcesol)
{
    ModuleBase::TITLE("surchem", "cal_force_sol");
    ModuleBase::timer::tick("surchem", "cal_force_sol");

    int nat = GlobalC::ucell.nat;
	ModuleBase::matrix force1(nat, 3);
    ModuleBase::matrix force2(nat, 3);
    
    force_cor_one(cell, rho_basis,force1);
    force_cor_two(cell, rho_basis,force2);
    
    int iat = 0;
    for (int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;ia++)
		{
            for(int ipol = 0; ipol < 3; ipol++)
            {
                forcesol(iat, ipol) = 0.5*force1(iat, ipol) + force2 (iat, ipol);
            }
				
		    ++iat;
        }
    }
    
    Parallel_Reduce::reduce_pool(forcesol.c, forcesol.nr * forcesol.nc);
    ModuleBase::timer::tick("surchem", "cal_force_sol");
    return;
}
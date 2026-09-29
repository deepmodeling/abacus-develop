#include "surchem.h"
#include "sccs/experimental_gaussian.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"

#include <stdexcept>

void surchem::force_cor_one(const UnitCell& cell,
                            const ModulePW::PW_Basis* rho_basis,
                            const ModuleBase::matrix& vloc,
                            ModuleBase::matrix& forcesol)
{
   
   
    //delta phi multiply by the derivative of nuclear charge density with respect to the positions
    std::complex<double> *N = new std::complex<double>[rho_basis->npw];
    std::complex<double> *vloc_at = new std::complex<double>[rho_basis->npw];
    std::complex<double> *delta_phi_g = new std::complex<double>[rho_basis->npw];
    //ModuleBase::GlobalFunc::ZEROS(delta_phi_g, rho_basis->npw);

    rho_basis->real2recip(this->delta_phi, delta_phi_g);
    // double Ael = 0;
    // double Ael1 = 0;
    // ModuleBase::GlobalFunc::ZEROS(vg, ngmc);
    int iat = 0;
    const int ig0 = rho_basis->ig_gge0; 
    for (int it = 0;it < cell.ntype;it++)
    {
        for (int ia = 0;ia < cell.atoms[it].na ; ia++)
        {
            for (int ig = 0; ig < rho_basis->npw; ig++)
            {   
                std::complex<double> phase = exp( ModuleBase::NEG_IMAG_UNIT *ModuleBase::TWO_PI * ( rho_basis->gcar[ig] * cell.atoms[it].tau[ia]));
                //vloc for each atom
                vloc_at[ig] = vloc(it, rho_basis->ig2igg[ig]) * phase;
                if(ig==ig0)
                {
                    N[ig] = cell.atoms[it].ncpp.zv / cell.omega;
                    continue; // skip G=0
                }
                const double fac = ModuleBase::e2 * ModuleBase::FOUR_PI /
                            (cell.tpiba2 * rho_basis->gg[ig]);
                N[ig] = -vloc_at[ig] / fac;
                
                //force for each atom
                forcesol(iat, 0) += rho_basis->gcar[ig][0] * imag(conj(delta_phi_g[ig]) * N[ig]);
                forcesol(iat, 1) += rho_basis->gcar[ig][1] * imag(conj(delta_phi_g[ig]) * N[ig]);
                forcesol(iat, 2) += rho_basis->gcar[ig][2] * imag(conj(delta_phi_g[ig]) * N[ig]);
            }
                
                forcesol(iat, 0) *= (cell.tpiba * cell.omega);
                forcesol(iat, 1) *= (cell.tpiba * cell.omega);
                forcesol(iat, 2) *= (cell.tpiba * cell.omega);
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

void surchem::force_cor_two(const UnitCell& cell,
                            const ModulePW::PW_Basis* rho_basis,
                            const int nspin,
                            ModuleBase::matrix& forcesol)
{
   
    std::complex<double> *n_pseudo = new std::complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(n_pseudo,rho_basis->npw);

    // this->gauss_charge(cell, pwb, n_pseudo);

    double *Vcav_sum =  new double[rho_basis->nrxx];
    ModuleBase::GlobalFunc::ZEROS(Vcav_sum, rho_basis->nrxx);
    std::complex<double> *Vcav_g = new std::complex<double>[rho_basis->npw];
    std::complex<double> *Vel_g = new std::complex<double>[rho_basis->npw];
    ModuleBase::GlobalFunc::ZEROS(Vcav_g, rho_basis->npw);
    ModuleBase::GlobalFunc::ZEROS(Vel_g, rho_basis->npw);
    for(int is=0; is<nspin; is++)
	{
		for (int ir=0; ir<rho_basis->nrxx; ir++)
		{
            Vcav_sum[ir] += this->Vcav(is, ir);
        }
    }

    rho_basis->real2recip(Vcav_sum, Vcav_g);
    rho_basis->real2recip(this->epspot, Vel_g);

    int iat = 0;
    // double Ael1 = 0;
    for (int it = 0;it < cell.ntype;it++)
    {
        double RCS = this->GetAtom.atom_RCS[cell.atoms[it].ncpp.psd];
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
                std::complex<double> phase = exp( ModuleBase::NEG_IMAG_UNIT *ModuleBase::TWO_PI * ( rho_basis->gcar[ig] * cell.atoms[it].tau[ia]));

                n_pseudo[ig].real((this->GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv)
                                  * phase.real() * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
                n_pseudo[ig].imag((this->GetAtom.atom_Z[cell.atoms[it].ncpp.psd] - cell.atoms[it].ncpp.zv)
                                  * phase.imag() * exp(-0.5 * gg * (sigma_rc_k * sigma_rc_k)));
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

                forcesol(iat, 0) *= (cell.tpiba * cell.omega);
                forcesol(iat, 1) *= (cell.tpiba * cell.omega);
                forcesol(iat, 2) *= (cell.tpiba * cell.omega);
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

void surchem::cal_force_sol(const UnitCell& cell,
                            const ModulePW::PW_Basis* rho_basis,
                            const ModuleBase::matrix& vloc,
                            const int nspin,
                            ModuleBase::matrix& forcesol)
{
    ModuleBase::TITLE("surchem", "cal_force_sol");
    ModuleBase::timer::start("surchem", "cal_force_sol");

    if (this->uses_sccs() || this->uses_pcc())
    {
        // Report SCCS/PCC failures through the standard ABACUS error path.
        try
        {
            if (forcesol.nr != cell.nat || forcesol.nc != 3)
            {
                throw std::invalid_argument("SCCS/PCC force matrix must have nat rows and three columns");
            }
            if (this->uses_sccs())
            {
                if (rho_basis == nullptr)
                {
                    throw std::invalid_argument("SCCS force requires an initialized PW basis");
                }
                this->cal_force_sccs(cell, *rho_basis, forcesol);
            }
            else
            {
                ModuleBase::GlobalFunc::ZEROS(forcesol.c, forcesol.nr * forcesol.nc);
                this->cal_force_pcc(cell, forcesol);
            }
        }
        catch (const std::exception& error)
        {
            ModuleBase::WARNING_QUIT("surchem::cal_force_sol", error.what());
        }
        ModuleBase::timer::end("surchem", "cal_force_sol");
        return;
    }

    int nat = cell.nat;
	ModuleBase::matrix force1(nat, 3);
    ModuleBase::matrix force2(nat, 3);
    
    force_cor_one(cell, rho_basis, vloc, force1);
    force_cor_two(cell, rho_basis, nspin, force2);
    
    int iat = 0;
    for (int it = 0;it < cell.ntype;it++)
	{
		for (int ia = 0;ia < cell.atoms[it].na;ia++)
		{
            for(int ipol = 0; ipol < 3; ipol++)
            {
                forcesol(iat, ipol) = 0.5*force1(iat, ipol) + force2 (iat, ipol);
            }
				
		    ++iat;
        }
    }
    
    Parallel_Reduce::reduce_pool(forcesol.c, forcesol.nr * forcesol.nc);
    ModuleBase::timer::end("surchem", "cal_force_sol");
    return;
}

// Add the smooth-ion reaction derivative and the point-ion vacuum PCC force at
// fixed converged electronic density. The reaction energy of Gaussian ions is
// independent of their width while they stay inside the epsilon=1 region, so no
// ionic-shape term is added. Electronic basis/overlap terms remain in the normal
// PW/LCAO force machinery. Reduce distributed terms before adding replicated
// point-ion terms, so neither MPI replication nor Ha-to-Ry conversion doubles them.
void surchem::cal_force_sccs(const UnitCell& cell,
                             const ModulePW::PW_Basis& rho_basis,
                             ModuleBase::matrix& forcesol) const
{
    if (forcesol.nr != cell.nat || forcesol.nc != 3)
    {
        throw std::invalid_argument("SCCS force matrix must have nat rows and three columns");
    }
    const ModuleSccs::SccsConfig& config = this->parameters_.sccs_config;
    if (!this->sccs_state_.valid)
    {
        throw std::logic_error("SCCS force requires a converged SCCS state");
    }

    // The symmetric sqrt-CG response defines the reaction energy for every
    // boundary; its derivative with respect to the ionic source is the solved
    // reaction potential (PCC included through the preconditioner). A
    // reconstructed continuous polarization source would change the force on
    // a finite grid.
    const std::vector<double>& reaction_potential
        = this->sccs_result_.electrostatic.reaction_potential;
    const ModuleBase::matrix smooth_force_hartree
        = ModuleSccs::gaussian_ionic_force(cell, rho_basis, ModuleSccs::gaussian_ion_spread,
                                          reaction_potential);
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            forcesol(atom, direction) = 2.0 * smooth_force_hartree(atom, direction);
        }
    }
    // ENVIRON 'full' mode: the core-electron Gaussians move the cavity with the
    // ions. Their force contracts the cavity potential, the derivative of the
    // electrostatic and non-electrostatic energies with respect to the cavity
    // density, with the Gaussian derivative (ENVIRON dboundary_dions).
    if (config.core_electrons)
    {
        const std::vector<double>& electrostatic_cavity = this->sccs_result_.response.cavity_potential;
        const std::vector<double>& non_electrostatic = this->sccs_result_.non_electrostatic.density_potential;
        std::vector<double> cavity_potential(electrostatic_cavity.size());
        for (std::size_t index = 0; index < cavity_potential.size(); ++index)
        {
            cavity_potential[index] = electrostatic_cavity[index] + non_electrostatic[index];
        }
        const ModuleBase::matrix core_force_hartree
            = ModuleSccs::gaussian_core_force(cell, rho_basis, config.core_spread, cavity_potential);
        for (int atom = 0; atom < cell.nat; ++atom)
        {
            for (int direction = 0; direction < 3; ++direction)
            {
                forcesol(atom, direction) += 2.0 * core_force_hartree(atom, direction);
            }
        }
    }
    Parallel_Reduce::reduce_pool(forcesol.c, forcesol.nr * forcesol.nc);
    if (config.boundary == ModuleSccs::Boundary::Periodic)
    {
        return;
    }

    this->cal_force_pcc(cell, forcesol);
}

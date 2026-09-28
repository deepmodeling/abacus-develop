#include "surchem.h"
#include "pcc/sccs_pcc_2d.h"
#include "pcc/sccs_pcc_coulomb.h"
#include "sccs/sccs_pw_charge.h"
#include "sccs/sccs_pw_force.h"
#include "sccs/experimental_gaussian.h"
#include "sccs/sccs_pw_coulomb.h"
#include "sccs/sccs_pw_reduction.h"
#include "pcc/sccs_pcc_2d_coulomb.h"
#include "source_base/parallel_reduce.h"
#include "source_base/timer.h"

#include <stdexcept>
#include <memory>

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

    if (this->uses_sccs())
    {
        try
        {
            if (rho_basis == nullptr)
            {
                throw std::invalid_argument("SCCS force requires an initialized PW basis");
            }
            this->cal_force_sccs(cell, *rho_basis, vloc, forcesol);
        }
        catch (...)
        {
            ModuleBase::timer::end("surchem", "cal_force_sol");
            throw;
        }
        ModuleBase::timer::end("surchem", "cal_force_sol");
        return;
    }
    if (this->uses_pcc())
    {
        ModuleBase::GlobalFunc::ZEROS(forcesol.c, forcesol.nr * forcesol.nc);
        this->cal_force_pcc(cell, forcesol);
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

// Add explicit smooth-ion and PCC ionic-shape derivatives at fixed converged
// electronic density. Electronic basis/overlap terms remain in the normal
// PW/LCAO force machinery. Reduce distributed terms before adding replicated
// point-ion terms, so neither MPI replication nor Ha-to-Ry conversion doubles them.
void surchem::cal_force_sccs(const UnitCell& cell,
                             const ModulePW::PW_Basis& rho_basis,
                             const ModuleBase::matrix& vloc,
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

    // Reproduce Environ's continuous dielectric polarization charge.
    // Finite-grid chain discretization does not make C[rho_pol] identical
    // to phi-phi_vac. This changes the ionic force only, not the energy.
    const std::vector<double> polarization = ModuleSccs::continuum_polarization_charge(
        this->sccs_result_.charge.solute, this->sccs_result_.response);
    const double volume_element = cell.omega / static_cast<double>(rho_basis.nxyz);
    const ModuleSccs::PoolChargeReduction charge_reduction;
    const std::vector<ModuleBase::Vector3<double>>& positions = this->fixed_source_cache_.positions;
    std::unique_ptr<ModuleSccs::CoulombOperator> coulomb;
    if (config.boundary == ModuleSccs::Boundary::Pcc0d)
    {
        coulomb.reset(new ModuleSccs::PccCoulombOperator(rho_basis, cell.tpiba,
                                                        positions, volume_element,
                                                        this->pcc_geometry_, charge_reduction));
    }
    else if (config.boundary == ModuleSccs::Boundary::Pcc2d)
    {
        coulomb.reset(new ModuleSccs::Pcc2dCoulombOperator(rho_basis, cell.tpiba,
                                                          positions, volume_element,
                                                          this->pcc_2d_geometry_, charge_reduction));
    }
    else
    {
        coulomb.reset(new ModuleSccs::PeriodicCoulombOperator(rho_basis, cell.tpiba));
    }
    std::vector<double> polarization_potential;
    coulomb->apply_potential(polarization, polarization_potential);
    const double gaussian_width = 0.5;
    const ModuleBase::matrix smooth_force_hartree
        = ModuleSccs::gaussian_ionic_force(cell, rho_basis, gaussian_width,
                                          polarization_potential);
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            forcesol(atom, direction) = 2.0 * smooth_force_hartree(atom, direction);
        }
    }
    // The 2D ionic-shape coefficient also depends explicitly on ion positions.
    // Differentiate its centered smooth and point ionic moments at fixed polarization.
    std::vector<double> point_shape_force_y(cell.nat, 0.0);
    if (config.boundary == ModuleSccs::Boundary::Pcc2d)
    {
        const ModuleSccs::Pcc2dGeometry& geometry = this->pcc_2d_geometry_;
        const double volume_element = cell.omega / static_cast<double>(rho_basis.nxyz);
        const std::vector<ModuleBase::Vector3<double>> positions
            = ModuleSccs::pw_grid_positions(rho_basis, cell.latvec, cell.lat0);
        const ModuleSccs::PoolChargeReduction reduction;
        const ModuleSccs::Pcc2dMoments smooth
            = ModuleSccs::reduced_pcc_2d_density_moments(this->sccs_result_.charge.ionic,
                                                        positions, volume_element, geometry, reduction);
        const ModuleSccs::Pcc2dMoments& point = this->pcc_ionic_moments_2d_;
        const double center = point.dipole_y / point.charge;
        const double shift = (smooth.dipole_y - point.dipole_y
                              - center * (smooth.charge - point.charge)) / point.charge;
        const double factor = ModuleBase::PI * this->sccs_result_.polarization_moments_2d.charge
                              / (geometry.parameters.periodic_area * geometry.parameters.cell_length_y);
        std::vector<double> shape_potential(rho_basis.nrxx);
        for (int ir = 0; ir < rho_basis.nrxx; ++ir)
        {
            const double relative_y = ModuleSccs::pcc_2d_relative_y(positions[ir].y, geometry);
            const double displacement = relative_y - center;
            shape_potential[ir] = factor * displacement * displacement;
        }
        const ModuleBase::matrix shape_force
            = ModuleSccs::gaussian_ionic_force(cell, rho_basis, gaussian_width, shape_potential);
        int iat = 0;
        for (int type = 0; type < cell.ntype; ++type)
        {
            for (int atom = 0; atom < cell.atoms[type].na; ++atom)
            {
                const double position_y = cell.atoms[type].tau[atom].y * cell.lat0;
                const double relative_y = ModuleSccs::pcc_2d_relative_y(position_y, geometry);
                point_shape_force_y[iat]
                    = 2.0 * factor * cell.atoms[type].ncpp.zv * (relative_y - center + shift);
                for (int direction = 0; direction < 3; ++direction)
                {
                    forcesol(iat, direction) += 2.0 * shape_force(iat, direction);
                }
                ++iat;
            }
        }
    }
    Parallel_Reduce::reduce_pool(forcesol.c, forcesol.nr * forcesol.nc);
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        forcesol(atom, 1) += 2.0 * point_shape_force_y[atom];
    }
    if (config.boundary == ModuleSccs::Boundary::Periodic)
    {
        return;
    }

    this->cal_force_pcc(cell, forcesol);
}

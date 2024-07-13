#include "esolver_lj.h"

#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_io/output_log.h"

namespace ModuleESolver
{

    void ESolver_LJ::before_all_runners(Input& inp, UnitCell& ucell)
    {
        ucell_ = &ucell;
        lj_potential = 0;
        lj_force.create(ucell.nat, 3);
        lj_virial.create(3, 3);

        // determine the maximum rcut and lj_rcut
        rcut_search_radius(inp.mdp.lj_rcut);

        // determine the LJ parameters
        set_c6_c12(inp.mdp.lj_rule, inp.mdp.lj_epsilon, inp.mdp.lj_sigma);

        // calculate the energy shift so that LJ energy is zero at rcut
        cal_en_shift(inp.mdp.lj_eshift);
    }

    void ESolver_LJ::runner(const int istep, UnitCell& ucell)
    {
        Grid_Driver grid_neigh(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
        atom_arrange::search(
            GlobalV::SEARCH_PBC,
            GlobalV::ofs_running,
            grid_neigh,
            ucell,
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);

        double distance=0.0;
        int index = 0;

        // Important! potential, force, virial must be zero per step
        lj_potential = 0;
        lj_force.zero_out();
        lj_virial.zero_out();

        ModuleBase::Vector3<double> tau1, tau2, dtau;
        for (int it = 0; it < ucell.ntype; ++it)
        {
            Atom* atom1 = &ucell.atoms[it];
            for (int ia = 0; ia < atom1->na; ++ia)
            {
                tau1 = atom1->tau[ia];
                grid_neigh.Find_atom(ucell, tau1, it, ia);
                for (int ad = 0; ad < grid_neigh.getAdjacentNum(); ++ad)
                {
                    tau2 = grid_neigh.getAdjacentTau(ad);
                    int it2 = grid_neigh.getType(ad);
                    dtau = (tau1 - tau2) * ucell.lat0;
                    distance = dtau.norm();
                    if (distance < lj_rcut(it, it2))
                    {
                        lj_potential += LJ_energy(distance, it, it2) - en_shift(it, it2);
                        ModuleBase::Vector3<double> f_ij = LJ_force(dtau, it, it2);
                        lj_force(index, 0) += f_ij.x;
                        lj_force(index, 1) += f_ij.y;
                        lj_force(index, 2) += f_ij.z;
                        LJ_virial(f_ij, dtau);
                    }
                }
                index++;
            }
        }

        lj_potential /= 2.0;
        GlobalV::ofs_running << " final etot is " << std::setprecision(11) << lj_potential * ModuleBase::Ry_to_eV
                             << " eV" << std::endl;

        // Post treatment for virial
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                lj_virial(i, j) /= (2.0 * ucell.omega);
            }
        }
#ifdef __MPI
        atom_arrange::delete_vector(
            GlobalV::ofs_running,
            GlobalV::SEARCH_PBC,
            grid_neigh,
            ucell, 
            GlobalV::SEARCH_RADIUS,
            GlobalV::test_atom_input);
#endif
    }

    double ESolver_LJ::cal_energy()
    {
        return lj_potential;
    }

    void ESolver_LJ::cal_force(ModuleBase::matrix& force)
    {
        force = lj_force;
        ModuleIO::print_force(GlobalV::ofs_running, *ucell_, "TOTAL-FORCE (eV/Angstrom)", force, false);
    }

    void ESolver_LJ::cal_stress(ModuleBase::matrix& stress)
    {
        stress = lj_virial;

        // external stress
        double unit_transform = ModuleBase::RYDBERG_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
        double external_stress[3] = {GlobalV::PRESS1, GlobalV::PRESS2, GlobalV::PRESS3};
        for (int i = 0; i < 3; i++)
        {
            stress(i, i) -= external_stress[i] / unit_transform;
        }

        ModuleIO::print_stress("TOTAL-STRESS", stress, true, false);
    }

    void ESolver_LJ::after_all_runners()
    {
        GlobalV::ofs_running << "\n\n --------------------------------------------" << std::endl;
        GlobalV::ofs_running << std::setprecision(16);
        GlobalV::ofs_running << " !FINAL_ETOT_IS " << lj_potential * ModuleBase::Ry_to_eV << " eV" << std::endl;
        GlobalV::ofs_running << " --------------------------------------------\n\n" << std::endl;
    }

    double ESolver_LJ::LJ_energy(const double d, const int i, const int j)
    {
        assert(d > 1e-6); // avoid atom overlap
        const double r2 = d * d;
        const double r4 = r2 * r2;
        const double r6 = r2 * r4;
        return lj_c12(i, j) / (r6 * r6) - lj_c6(i, j) / r6;
    }

    ModuleBase::Vector3<double> ESolver_LJ::LJ_force(const ModuleBase::Vector3<double> dr, const int i, const int j)
    {
        const double d = dr.norm();
        assert(d > 1e-6); // avoid atom overlap
        const double r2 = d * d;
        const double r4 = r2 * r2;
        const double r8 = r4 * r4;
        const double r14 = r8 * r4 * r2;
        double coff = 12.0 * lj_c12(i, j) / r14 - 6.0 * lj_c6(i, j) / r8;
        return dr * coff;
    }

    void ESolver_LJ::LJ_virial(const ModuleBase::Vector3<double>& force,
        const ModuleBase::Vector3<double>& dtau)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                lj_virial(i, j) += dtau[i] * force[j];
            }
        }
    }

    void ESolver_LJ::rcut_search_radius(const std::vector<double>& rcut)
    {
        const int ntype = this->ucell_->ntype;
        lj_rcut.create(ntype, ntype);
        double rcut_max = 0.0;

        if (rcut.size() == 1)
        {
            rcut_max = rcut[0] * ModuleBase::ANGSTROM_AU;
            for (int i = 0; i < ntype; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    lj_rcut(i, j) = rcut_max;
                    lj_rcut(j, i) = rcut_max;
                }
            }
        }
        else if (rcut.size() == ntype * (ntype + 1) / 2)
        {
            for (int i = 0; i < ntype; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    int k = i * (i + 1) / 2 + j;
                    lj_rcut(i, j) = rcut[k] * ModuleBase::ANGSTROM_AU;
                    lj_rcut(j, i) = lj_rcut(i, j);
                    rcut_max = std::max(rcut_max, lj_rcut(i, j));
                }
            }
        }

        // set the search radius
        GlobalV::SEARCH_RADIUS = rcut_max + 0.01;
    }

    void ESolver_LJ::set_c6_c12(const int rule, const std::vector<double> epsilon, const std::vector<double> sigma)
    {
        const int ntype = this->ucell_->ntype;
        lj_c6.create(ntype, ntype);
        lj_c12.create(ntype, ntype);

        std::vector<double> lj_epsilon = epsilon;
        std::vector<double> lj_sigma = sigma;

        std::transform(begin(lj_epsilon), end(lj_epsilon), begin(lj_epsilon), [](double x) {
            return x / ModuleBase::Ry_to_eV;
        });
        std::transform(begin(lj_sigma), end(lj_sigma), begin(lj_sigma), [](double x) {
            return x * ModuleBase::ANGSTROM_AU;
        });

        if (lj_epsilon.size() != lj_sigma.size())
        {
            ModuleBase::WARNING_QUIT("ESolver_LJ", " the number of lj_epsilon should be equal to lj_sigma ");
        }
        // do not need any combination rules
        else if (lj_sigma.size() == ntype * (ntype + 1) / 2)
        {
            for (int i = 0; i < ntype; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    int k = i * (i + 1) / 2 + j;
                    double temp = pow(lj_sigma[k], 6);
                    lj_c6(i, j) = 4.0 * lj_epsilon[k] * temp;
                    lj_c12(i, j) = lj_c6(i, j) * temp;
                    lj_c6(j, i) = lj_c6(i, j);
                    lj_c12(j, i) = lj_c12(i, j);
                }
            }
        }
        // combination rule 1
        else if (lj_sigma.size() == ntype && rule == 1)
        {
            for (int i = 0; i < ntype; i++)
            {
                // first determine the diagonal elements
                double temp = pow(lj_sigma[i], 6);
                lj_c6(i, i) = 4.0 * lj_epsilon[i] * temp;
                lj_c12(i, i) = lj_c6(i, i) * temp;

                // then determine the non-diagonal elements
                for (int j = 0; j < i; j++)
                {
                    lj_c6(i, j) = std::sqrt(lj_c6(i, i) * lj_c6(j, j));
                    lj_c12(i, j) = std::sqrt(lj_c12(i, i) * lj_c12(j, j));
                    lj_c6(j, i) = lj_c6(i, j);
                    lj_c12(j, i) = lj_c12(i, j);
                }
            }
        }
        // combination rule 2
        else if (lj_sigma.size() == ntype && rule == 2)
        {
            for (int i = 0; i < ntype; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sigma_ij = (lj_sigma[i] + lj_sigma[j]) / 2.0;
                    double epsilon_ij = std::sqrt(lj_epsilon[i] * lj_epsilon[j]);

                    double temp = pow(sigma_ij, 6);
                    lj_c6(i, j) = 4.0 * epsilon_ij * temp;
                    lj_c12(i, j) = lj_c6(i, j) * temp;
                    lj_c6(j, i) = lj_c6(i, j);
                    lj_c12(j, i) = lj_c12(i, j);
                }
            }
        }
    }

    void ESolver_LJ::cal_en_shift(const bool is_shift)
    {
        const int ntype = this->ucell_->ntype;
        en_shift.create(ntype, ntype);

        if (is_shift)
        {
            for (int i = 0; i < ntype; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    en_shift(i, j) = LJ_energy(lj_rcut(i, j), i, j);
                    en_shift(j, i) = en_shift(i, j);
                }
            }
        }
    }
}

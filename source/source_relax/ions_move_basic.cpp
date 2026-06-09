#include "ions_move_basic.h"

#include <algorithm>
#include "source_io/module_parameter/parameter.h"
#include "source_base/global_function.h"
#include "source_base/global_variable.h"
#include "source_cell/update_cell.h"
#include "source_cell/print_cell.h"

// Ions-specific parameters (shared variables are in Relax_Data)
int Ions_Move_Basic::update_iter = 0;
double Ions_Move_Basic::trust_radius = 0.0;
double Ions_Move_Basic::trust_radius_old = 0.0;
double Ions_Move_Basic::relax_bfgs_rmax = -1.0; // default is 0.8
double Ions_Move_Basic::relax_bfgs_rmin = -1.0; // default is 1e-5
double Ions_Move_Basic::relax_bfgs_init = -1.0; // default is 0.5
double Ions_Move_Basic::best_xxx = 1.0;
int Ions_Move_Basic::out_stru = 0;
std::vector<std::string> Ions_Move_Basic::relax_method = {"bfgs","2"};

void Ions_Move_Basic::setup_gradient(const UnitCell &ucell, const ModuleBase::matrix &force, double *pos, double *grad, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Basic", "setup_gradient");

    assert(ucell.ntype > 0);
    assert(pos != nullptr);
    assert(grad != nullptr);
    assert(dim == 3 * ucell.nat);

    std::fill_n(pos, dim, 0.0);
    std::fill_n(grad, dim, 0.0);

    // (1) init gradient
    // the unit of pos: Bohr.
    // the unit of force: Ry/Bohr.
    // the unit of gradient:
    int iat = 0;
    for (int it = 0; it < ucell.ntype; it++)
    {
        Atom *atom = &ucell.atoms[it];
        for (int ia = 0; ia < ucell.atoms[it].na; ia++)
        {
            for (int ik = 0; ik < 3; ++ik)
            {
                pos[3 * iat + ik] = atom->tau[ia][ik] * ucell.lat0;
                if (atom->mbl[ia][ik])
                {
                    grad[3 * iat + ik] = -force(iat, ik) * ucell.lat0;
                }
            }
            ++iat;
        }
    }

    return;
}

void Ions_Move_Basic::move_atoms(UnitCell &ucell, double *move, double *pos, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Basic", "move_atoms");

    assert(move != nullptr);
    assert(pos != nullptr);

    //------------------------
    // for test only
    //------------------------
    if (PARAM.inp.test_relax_method)
    {
        int iat = 0;
        ofs << "\n movement of ions (unit is Bohr) : " << std::endl;
        ofs << " " << std::setw(12) << "Atom" << std::setw(15) << "x" << std::setw(15) << "y"
                             << std::setw(15) << "z" << std::endl;
        for (int it = 0; it < ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                std::stringstream ss;
                ss << "move_" << ucell.atoms[it].label << ia + 1;
                ofs << " " << std::setw(12) << ss.str().c_str() << std::setw(15) << move[3 * iat + 0]
                                     << std::setw(15) << move[3 * iat + 1] << std::setw(15) << move[3 * iat + 2]
                                     << std::endl;
                iat++;
            }
        }
        assert(iat == ucell.nat);
    }

    const double move_threshold = 1.0e-10;
    const int total_freedom = ucell.nat * 3;

    if (ModuleSymmetry::Symmetry::symm_flag && ucell.symm.all_mbl && ucell.symm.nrotk > 0) 
    {
        ucell.symm.symmetrize_vec3_nat(move);
    }

    for (int i = 0; i < total_freedom; i++)
    {
        if (std::abs(move[i]) > move_threshold)
        {
            pos[i] += move[i];
        }
    }
    unitcell::update_pos_tau(ucell.lat,pos,ucell.ntype,ucell.nat,ucell.atoms);

    //--------------------------------------------
    // Print out the structure file.
    //--------------------------------------------
    unitcell::print_tau(ucell.atoms,ucell.Coordinate,ucell.ntype,ucell.lat0,ofs);

    return;
}

void Ions_Move_Basic::check_converged(const UnitCell &ucell, const double *grad, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Basic", "check_converged");
    assert(dim > 0);

    //------------------------------------------------
    // check the gradient value
    //------------------------------------------------
    Ions_Move_Basic::largest_grad = 0.0;
    for (int i = 0; i < dim; i++)
    {
        if (Ions_Move_Basic::largest_grad < std::abs(grad[i]))
        {
            Ions_Move_Basic::largest_grad = std::abs(grad[i]);
        }
    }
    // mohan add 2010-08-06
    Ions_Move_Basic::largest_grad /= ucell.lat0;

    if (PARAM.inp.test_relax_method)
    {
        ModuleBase::GlobalFunc::OUT(ofs, "old total energy (ry)", etot_p);
        ModuleBase::GlobalFunc::OUT(ofs, "new total energy (ry)", etot);
        ModuleBase::GlobalFunc::OUT(ofs, "energy difference (ry)", Ions_Move_Basic::ediff);
        ModuleBase::GlobalFunc::OUT(ofs, "largest gradient (ry/bohr)", Ions_Move_Basic::largest_grad);
    }

    if (PARAM.inp.out_level == "ie")
    {
        std::cout << " ETOT DIFF (eV)       : " << Ions_Move_Basic::ediff * ModuleBase::Ry_to_eV << std::endl;
        std::cout << " LARGEST GRAD (eV/Angstrom)  : " << Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
                  << std::endl;

        ofs << "\n Largest force is " << largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A
                             << " eV/Angstrom while threshold is " 
                             << PARAM.inp.force_thr_ev << " eV/Angstrom" << std::endl;
    }

    const double etot_diff = std::abs(Ions_Move_Basic::ediff);

    // need to update, mohan 2010-07-10
    const double etot_thr = 1.0e-3; // Rydeberg.

    if (Ions_Move_Basic::largest_grad == 0.0)
    {
        ofs << " largest force is 0, no movement is possible." << std::endl;
        ofs << " it may converged, otherwise no movement of atom is allowed." << std::endl;
        Ions_Move_Basic::converged = true;
    }
    // mohan update 2011-04-21
    else if (etot_diff < etot_thr && Ions_Move_Basic::largest_grad < PARAM.inp.force_thr )
    {
        ofs << "\n Ion relaxation is converged!" << std::endl;
        ofs << "\n Energy difference (Ry) = " << etot_diff << std::endl;

        Ions_Move_Basic::converged = true;
        ++Ions_Move_Basic::update_iter;
    }
    else
    {
        ofs << "\n Ion relaxation is not converged yet (threshold is "
                             << PARAM.inp.force_thr  * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A << ")" << std::endl;
        // std::cout << "\n etot_diff=" << etot_diff << " etot_thr=" << etot_thr
        //<< " largest_grad=" << largest_grad << " force_thr=" << PARAM.inp.force_thr  << std::endl;
        Ions_Move_Basic::converged = false;
    }

    return;
}

void Ions_Move_Basic::terminate(const UnitCell &ucell, std::ofstream& ofs)
{
    ModuleBase::TITLE("Ions_Move_Basic", "terminate");
    if (Ions_Move_Basic::converged)
    {
        ofs << " end of geometry optimization" << std::endl;
        ModuleBase::GlobalFunc::OUT(ofs, "istep", Ions_Move_Basic::istep);
        ModuleBase::GlobalFunc::OUT(ofs, "update iteration", Ions_Move_Basic::update_iter);
        /*
        ofs<<"Saving the approximate inverse hessian"<<std::endl;
        std::ofstream hess("hess.out");
        for(int i=0;i<dim;i++)
        {
            for(int j=0;j<dim;j++)
            {
                hess << inv_hess(i,j);
            }
        }
        hess.close();
        */
    }
    else
    {
        ofs << " the maximum number of steps has been reached." << std::endl;
        ofs << " end of geometry optimization." << std::endl;
    }

    //-----------------------------------------------------------
    // Print the structure.
    //-----------------------------------------------------------
    unitcell::print_tau(ucell.atoms,ucell.Coordinate,ucell.ntype,ucell.lat0,ofs);
    return;
}

void Ions_Move_Basic::setup_etot(const double &energy_in, const bool judgement)
{
    if (Ions_Move_Basic::istep == 1)
    {
        // p == previous
        Ions_Move_Basic::etot_p = energy_in;
        Ions_Move_Basic::etot = energy_in;
        ediff = etot - etot_p;
    }
    else
    {
        // mohan modify 2010-07-10
        // mohan modify again 2010-07-25
        if (judgement) // for sd
        {
            Ions_Move_Basic::etot = energy_in;
            if (Ions_Move_Basic::etot_p > etot)
            {
                ediff = etot - etot_p;
                Ions_Move_Basic::etot_p = etot;
            }
            else
            {
                // this step will not be accepted
                ediff = 0.0;
            }
        }
        // note: the equlibrium point is not
        // need to be the smallest energy point.
        else // for bfgs
        {
            Ions_Move_Basic::etot_p = etot;
            Ions_Move_Basic::etot = energy_in;
            ediff = etot - etot_p;
        }
    }

    return;
}

double Ions_Move_Basic::dot_func(const double *a, const double *b, const int &dim_in)
{
    double result = 0.0;
    for (int i = 0; i < dim_in; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

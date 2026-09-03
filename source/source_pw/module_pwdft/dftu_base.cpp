#include "source_pw/module_pwdft/dftu_base.h"

#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base_io.h"
#include "source_base/global_function.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

// local inline helpers for eigenvalue calculation (JacobiRotate, CalculateEigenvalues)
// have been migrated to dftu_base_io.cpp, where they are used by DFTU_BASE::write_occup_m.
// mohan refactored 2025-11-08
// All members are now non-static; default values are in the header.


Plus_U_Base::Plus_U_Base()
{}


Plus_U_Base::~Plus_U_Base()
{}


// init_base: base-only initialization shared by PW and LCAO paths.
// LCAO-specific setup (paraV, orb, ucell pointer) stays in Plus_U::init().
void Plus_U_Base::init_base(UnitCell& cell,
                             const int npol,
                             const int nspin,
                             const std::vector<int>& orbital_corr,
                             const bool yukawa_potential,
                             const double yukawa_lambda,
                             const std::string& global_readin_dir,
                             const std::string& global_out_dir,
                             const std::string& init_chg,
                             const std::string& device,
                             const int kpar,
                             const std::vector<double>& hubbard_u,
                             const double uramping,
                             const int occ_mat_ctrl,
                             const int mixing_dftu)
{
    ModuleBase::TITLE("Plus_U_Base", "init_base");

#ifndef __MPI
    std::cout << "DFT+U module is only accessible in mpi version" << std::endl;
    exit(0);
#endif

    this->nspin = nspin;
    this->orbital_corr = orbital_corr;
    this->uramping = uramping;
    this->occ_mat_ctrl = occ_mat_ctrl;
    this->mixing_dftu = mixing_dftu;
    this->u_target = hubbard_u;
    this->u_current = hubbard_u;
    if (uramping > 0.01)
    {
        std::fill(this->u_current.begin(),
                  this->u_current.end(),
                  0.0);
    }
    this->device = device;
    this->kpar = kpar;

    this->energy_u = 0.0;

    this->occmat_.init(cell, orbital_corr, nspin, npol);

    this->pot_uterm_pw_index.resize(cell.nat);
    int pot_index = 0;

    int num_locale = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            const int iat = cell.itia2iat(it, ia);

            if(!has_correlated_orbital(it))
            {
                continue;
            }

            const int tlp1_npol = (get_orbital_corr(it)*2+1)*npol;
            const int tlp1 = 2 * get_orbital_corr(it) + 1;
            const int elem_size = tlp1 * tlp1;
            if(nspin == 4)
            {
                this->pot_uterm_pw_index[iat] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                this->pot_uterm_pw_index[iat] = pot_index;
                pot_index += elem_size;
            }

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                const int N = cell.atoms[it].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (nspin == 1 || nspin == 2)
                    {
                        num_locale += (2 * l + 1) * (2 * l + 1) * 2;
                    }
                    else if (nspin == 4)
                    {
                        num_locale += (2 * l + 1) * (2 * l + 1) * npol * npol;
                    }
                }
            }
        }
    }

    if (nspin == 2) pot_index *= 2;

    this->pot_uterm_pw.resize(pot_index, 0.0);
    this->uom_array.resize(pot_index, 0.0);
    this->uom_save.resize(pot_index, 0.0);

    if (yukawa_potential)
    {
        this->yukawa_.reset(new YukawaScreening());
        this->yukawa_->init(cell, orbital_corr, yukawa_lambda);
    }
    else
    {
        // Clear any stale object from a previous init_base() call with
        // yukawa_potential == true, preserving the old explicit-flag semantics.
        this->yukawa_.reset();
    }

    if (occ_mat_ctrl != 0)
    {
        std::stringstream sst;
        sst << global_readin_dir << "dm_onsite_ini.txt";
        DFTU_BASE::read_occup_m(cell, this->occmat_, this->orbital_corr, this->occ_mat_ctrl,
                                sst.str(), init_chg, nspin, npol);
#ifdef __MPI
        DFTU_BASE::local_occup_bcast(cell, this->occmat_, this->orbital_corr, nspin, npol);
#endif

        mark_occ_mat_initialized();
        this->occmat_.copy_to_save(cell, this->orbital_corr);
        this->occmat_.write_save_to_flat(cell, this->orbital_corr,
                                         this->pot_uterm_pw_index, this->uom_save);
    }
    else
    {
        if (init_chg == "file")
        {
            std::stringstream sst;
            sst << global_readin_dir << "dm_onsite.txt";
            DFTU_BASE::read_occup_m(cell, this->occmat_, this->orbital_corr, this->occ_mat_ctrl,
                                    sst.str(), init_chg, nspin, npol);
#ifdef __MPI
            DFTU_BASE::local_occup_bcast(cell, this->occmat_, this->orbital_corr, nspin, npol);
#endif
            mark_occ_mat_initialized();
        }
        else
        {
            this->occmat_.zero(cell, this->orbital_corr);
        }
    }

    ModuleBase::Memory::record("Plus_U_Base::occ_mat", sizeof(double) * num_locale);
    return;
}


void Plus_U_Base::uramping_update()
{
    // Yukawa calculates U directly every iteration, no need for ramping
    if (use_yukawa()) {
        return;
    }
    // if uramping < 0.1, use the original U
    if (this->uramping < 0.01) {
        return;
    }
    // loop to change U
    for (int i = 0; i < static_cast<int>(this->u_target.size()); i++)
    {
        if (this->u_current[i] + this->uramping < this->u_target[i])
        {
            this->u_current[i] += this->uramping;
        }
        else
        {
            this->u_current[i] = this->u_target[i];
        }
    }
}


bool Plus_U_Base::u_converged()
{
    // Yukawa calculates U directly every iteration, always considered converged
    if (use_yukawa()) {
        return true;
    }
    for (int i = 0; i < static_cast<int>(this->u_target.size()); i++)
    {
        if (this->u_current[i] != this->u_target[i])
        {
            return false;
        }
    }
    return true;
}


// cal_occ_pw() is implemented in source_pw/module_pwdft/dftu_base_occ.cpp
// as a Plus_U_Base method. Pure per-atom kernels live in dftu_base_tools.{h,cpp}
// as free functions in namespace DFTU_BASE.

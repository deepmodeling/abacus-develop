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
// LCAO callers validate the paraV square-matrix invariant before calling.
void Plus_U_Base::init_base(UnitCell& cell,
                             const int npol,
                             const int nspin,
                             const std::vector<int>& l_channel,
                             const bool yukawa_potential,
                             const double yukawa_lambda,
                             const std::string& global_readin_dir,
                             const std::string& global_out_dir,
                             const std::string& init_chg,
                             const std::string& device,
                             const std::vector<double>& hubbard_u,
                             const double uramping,
                             const int init_occ_mat,
                             const int mixing_dftu,
                             const DFTU_BASE::OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("Plus_U_Base", "init_base");

#ifndef __MPI
    ModuleBase::WARNING_QUIT("Plus_U_Base::init_base", "DFT+U module is only accessible in MPI version");
#endif

    this->l_channel = l_channel;
    this->uramping = uramping;
    this->init_occ_mat = init_occ_mat;
    this->u_target = hubbard_u;
    this->u_current = hubbard_u;
    if (uramping > 0.01)
    {
        std::fill(this->u_current.begin(),
                  this->u_current.end(),
                  0.0);
    }
    this->device = device;

    this->energy_u = 0.0;

    // The occupation-matrix file is read only once for the whole run.
    // occmat_.init() below reallocates and clears the storage on every
    // ionic step, so preserve the current in-memory matrix first and
    // restore it after init() when the file has already been loaded.
    OccupationMatrix occmat_kept;
    const bool restore_occmat = (init_occ_mat != 0 && this->occmat_file_loaded_);
    if (restore_occmat)
    {
        occmat_kept = this->occmat_;
    }

    this->occmat_.init(cell, l_channel, nspin, npol);

    this->uterm_mat_index.resize(cell.nat);
    int pot_index = 0;

    int num_locale = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            const int iat = cell.itia2iat(it, ia);

            const int target_l = this->l_channel[it];
            if (target_l == -1)
            {
                continue;
            }

            const int tlp1_npol = (target_l * 2 + 1) * npol;
            const int tlp1 = 2 * target_l + 1;
            const int elem_size = tlp1 * tlp1;
            if(nspin == 4)
            {
                this->uterm_mat_index[iat] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                this->uterm_mat_index[iat] = pot_index;
                pot_index += elem_size;
            }

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                // only the first radial channel (n=0) of each l is stored
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

    if (nspin == 2) pot_index *= 2;

    this->uterm_mat.resize(pot_index, 0.0);

    // construct the occupation-matrix mixer only when mixing is enabled
    if (mixing_dftu != 0)
    {
        this->occ_mixer_.reset(new OccMatMixer());
        this->occ_mixer_->init(&cell, &this->l_channel,
                               &this->uterm_mat_index, nspin, pot_index);
    }

    if (yukawa_potential)
    {
        this->yukawa_.reset(new YukawaScreening());
        this->yukawa_->init(cell, l_channel, yukawa_lambda);
    }
    else
    {
        // Clear any stale object from a previous init_base() call with
        // yukawa_potential == true, preserving the old explicit-flag semantics.
        this->yukawa_.reset();
    }

    if (init_occ_mat != 0)
    {
        if (this->occmat_file_loaded_)
        {
            // The file was read at the first ionic step; keep using the
            // in-memory occupation matrix instead of reading the file
            // again. For init_occ_mat=2 this matrix never changes; for
            // init_occ_mat=1 it carries the result of the previous
            // ionic step.
            this->occmat_ = occmat_kept;
            this->set_occmat_ready();
            this->occmat_.copy_to_save(cell, this->l_channel);
            if (this->has_occ_mixer())
            {
                // seed the freshly rebuilt mixer with the current matrix
                this->occ_mixer().seed_save(this->occmat_);
            }
        }
        else
        {
            // Try the user-prepared dm_onsite_ini.txt first, then fall back to
            // the occupation-matrix snapshot written by a previous run
            // (occ_mat.txt), and finally to the legacy dm_onsite.txt name.
            const std::vector<std::string> candidates = {"dm_onsite_ini.txt",
                                                          "occ_mat.txt",
                                                          "dm_onsite.txt"};
            const std::string readin_fn = DFTU_BASE::find_first_existing_file(global_readin_dir, candidates);
            if (readin_fn.empty())
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::init_base",
                                         "init_occ_mat is set but no occupation-matrix file found in "
                                         + global_readin_dir
                                         + ". Tried: dm_onsite_ini.txt, occ_mat.txt, dm_onsite.txt");
            }
            DFTU_BASE::read_occup_m(cell, this->occmat_, this->l_channel, this->init_occ_mat,
                                    readin_fn, init_chg, nspin, npol, soc_layout);
#ifdef __MPI
            DFTU_BASE::local_occup_bcast(cell, this->occmat_, this->l_channel, nspin, npol);
#endif

            this->set_occmat_ready();
            this->occmat_.copy_to_save(cell, this->l_channel);
            if (this->has_occ_mixer())
            {
                // seed the mixing history with the file-loaded occupation matrix
                this->occ_mixer().seed_save(this->occmat_);
            }
            this->occmat_file_loaded_ = true;
        }
    }
    else
    {
        if (init_chg == "file")
        {
            // occ_mat.txt is the current output name; fall back to the
            // legacy dm_onsite.txt so that output directories written by
            // older versions can still be used for restarts.
            const std::vector<std::string> candidates = {"occ_mat.txt", "dm_onsite.txt"};
            const std::string readin_fn = DFTU_BASE::find_first_existing_file(global_readin_dir, candidates);
            if (readin_fn.empty())
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::init_base",
                                         "init_chg is set to file but no occupation-matrix file "
                                         "found in "
                                         + global_readin_dir
                                         + ". Tried: occ_mat.txt, dm_onsite.txt. Please do an scf "
                                           "calculation first.");
            }
            DFTU_BASE::read_occup_m(cell, this->occmat_, this->l_channel, this->init_occ_mat,
                                    readin_fn, init_chg, nspin, npol, soc_layout);
#ifdef __MPI
            DFTU_BASE::local_occup_bcast(cell, this->occmat_, this->l_channel, nspin, npol);
#endif
            this->set_occmat_ready();
        }
        else
        {
            this->occmat_.zero(cell, this->l_channel);
        }
    }

    ModuleBase::Memory::record("Plus_U_Base::occ_mat", sizeof(double) * num_locale);
}


void Plus_U_Base::uramping_update()
{
    // Yukawa calculates U directly every iteration, no need for ramping
    if (this->yukawa_ != nullptr)
    {
        return;
    }
    // if uramping < 0.1, use the original U
    if (this->uramping < 0.01)
    {
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
    if (this->yukawa_ != nullptr)
    {
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


// cal_occ_pw() is implemented as free function DFTU_BASE::cal_occ_pw
// in source_pw/module_pwdft/dftu_pw.cpp.
// All pure per-atom kernels also live in dftu_pw.{h,cpp}
// as free functions in namespace DFTU_BASE.

#include "source_pw/module_pwdft/dftu_base.h"

#include "source_base/global_function.h"
#include "source_base/memory_recorder.h"
#include "source_base/parallel_global.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cmath>
#include <complex>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

// local inline helpers for eigenvalue calculation (JacobiRotate, CalculateEigenvalues)
// have been migrated to dftu_output.cpp, where they are used by dftu_io::write_occup_m.
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
    this->use_yukawa = yukawa_potential;
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

    this->occ_mat.resize(cell.nat);
    this->occ_mat_save.resize(cell.nat);
    this->eff_pot_pw_index.resize(cell.nat);
    int pot_index = 0;

    this->iatlnmipol2iwt.resize(cell.nat);

    int num_locale = 0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        for (int ia = 0; ia < cell.atoms[it].na; ia++)
        {
            const int iat = cell.itia2iat(it, ia);

            occ_mat[iat].resize(cell.atoms[it].nwl + 1);
            occ_mat_save[iat].resize(cell.atoms[it].nwl + 1);

            this->iatlnmipol2iwt[iat].resize(cell.atoms[it].nwl + 1);

            if(!has_correlated_orbital(it))
            {
                continue;
            }

            const int tlp1_npol = (get_orbital_corr(it)*2+1)*npol;
            const int tlp1 = 2 * get_orbital_corr(it) + 1;
            const int elem_size = tlp1 * tlp1;
            if(nspin == 4)
            {
                this->eff_pot_pw_index[iat] = pot_index;
                pot_index += tlp1_npol * tlp1_npol;
            }
            else
            {
                this->eff_pot_pw_index[iat] = pot_index;
                pot_index += elem_size;
            }

            for (int l = 0; l <= cell.atoms[it].nwl; l++)
            {
                const int N = cell.atoms[it].l_nchi[l];

                occ_mat[iat][l].resize(N);
                occ_mat_save[iat][l].resize(N);

                for (int n = 0; n < N; n++)
                {
                    if (nspin == 1 || nspin == 2)
                    {
                        occ_mat[iat][l][n].resize(2);
                        occ_mat_save[iat][l][n].resize(2);

                        occ_mat[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_mat[iat][l][n][1].create(2 * l + 1, 2 * l + 1);

                        occ_mat_save[iat][l][n][0].create(2 * l + 1, 2 * l + 1);
                        occ_mat_save[iat][l][n][1].create(2 * l + 1, 2 * l + 1);
                        num_locale += (2 * l + 1) * (2 * l + 1) * 2;
                    }
                    else if (nspin == 4)
                    {
                        occ_mat[iat][l][n].resize(1);
                        occ_mat_save[iat][l][n].resize(1);

                        occ_mat[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        occ_mat_save[iat][l][n][0].create((2 * l + 1) * npol, (2 * l + 1) * npol);
                        num_locale += (2 * l + 1) * (2 * l + 1) * npol * npol;
                    }
                }
            }

            this->iatlnmipol2iwt[iat].resize(cell.atoms[it].nwl + 1);
            for (int L = 0; L <= cell.atoms[it].nwl; L++)
            {
                this->iatlnmipol2iwt[iat][L].resize(cell.atoms[it].l_nchi[L]);

                for (int n = 0; n < cell.atoms[it].l_nchi[L]; n++)
                {
                    this->iatlnmipol2iwt[iat][L][n].resize(2 * L + 1);

                    for (int m = 0; m < 2 * L + 1; m++)
                    {
                        this->iatlnmipol2iwt[iat][L][n][m].resize(npol);
                    }
                }
            }

            for (int iw = 0; iw < cell.atoms[it].nw * npol; iw++)
            {
                int iw0 = iw / npol;
                int ipol = iw % npol;
                int iwt = cell.itiaiw2iwt(it, ia, iw);
                int l = cell.atoms[it].iw2l[iw0];
                int n = cell.atoms[it].iw2n[iw0];
                int m = cell.atoms[it].iw2m[iw0];

                this->iatlnmipol2iwt[iat][l][n][m][ipol] = iwt;
            }
        }
    }

    if (nspin == 2) pot_index *= 2;

    this->eff_pot_pw.resize(pot_index, 0.0);
    this->uom_array.resize(pot_index, 0.0);
    this->uom_save.resize(pot_index, 0.0);

    if (use_yukawa)
    {
        this->Fk.resize(cell.ntype);

        this->U_Yukawa.resize(cell.ntype);
        this->J_Yukawa.resize(cell.ntype);

        for (int it = 0; it < cell.ntype; it++)
        {
            const int NL = cell.atoms[it].nwl + 1;

            this->Fk[it].resize(NL);
            this->U_Yukawa[it].resize(NL);
            this->J_Yukawa[it].resize(NL);

            for (int l = 0; l < NL; l++)
            {
                int N = cell.atoms[it].l_nchi[l];

                this->Fk[it][l].resize(N);
                for (int n = 0; n < N; n++)
                {
                    this->Fk[it][l][n].resize(l + 1, 0.0);
                }

                this->U_Yukawa[it][l].resize(N, 0.0);
                this->J_Yukawa[it][l].resize(N, 0.0);
            }
        }
    }

    if (occ_mat_ctrl != 0)
    {
        std::stringstream sst;
        sst << global_readin_dir << "dm_onsite_ini.txt";
        this->read_occup_m(cell, sst.str(), init_chg, nspin, npol);
#ifdef __MPI
        this->local_occup_bcast(cell, nspin, npol);
#endif

        mark_occ_mat_initialized();
        this->copy_occ_mat(cell);
    }
    else
    {
        if (init_chg == "file")
        {
            std::stringstream sst;
            sst << global_readin_dir << "dm_onsite.txt";
            this->read_occup_m(cell, sst.str(), init_chg, nspin, npol);
#ifdef __MPI
            this->local_occup_bcast(cell, nspin, npol);
#endif
            mark_occ_mat_initialized();
        }
        else
        {
            this->zero_occ_mat(cell);
        }
    }

    ModuleBase::Memory::record("Plus_U_Base::occ_mat", sizeof(double) * num_locale);
    return;
}


void Plus_U_Base::uramping_update()
{
    // Yukawa calculates U directly every iteration, no need for ramping
    if (use_yukawa) {
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
    if (use_yukawa) {
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


// copy_locale — save current occ_mat to occ_mat_save and uom_save
void Plus_U_Base::copy_occ_mat(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "copy_occ_mat");
    ModuleBase::timer::start("Plus_U_Base", "copy_occ_mat");

    for (int T = 0; T < ucell.ntype; T++)
    {
        int target_l = get_orbital_corr(T);
        if (target_l == -1)
            continue;

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            if (this->nspin == 4)
            {
                occ_mat_save[iat][target_l][0][0] = occ_mat[iat][target_l][0][0];
                if(this->uom_save.size() != 0)
                {
                    const int size = occ_mat[iat][target_l][0][0].nr * occ_mat[iat][target_l][0][0].nc;
                    for(int mm=0; mm<size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat]+mm] = occ_mat[iat][target_l][0][0].c[mm];
                    }
                }
            }
            else if (this->nspin == 1 || this->nspin == 2)
            {
                occ_mat_save[iat][target_l][0][0] = occ_mat[iat][target_l][0][0];
                occ_mat_save[iat][target_l][0][1] = occ_mat[iat][target_l][0][1];
                if(this->uom_save.size() != 0)
                {
                    const int size = occ_mat[iat][target_l][0][0].nr * occ_mat[iat][target_l][0][0].nc;
                    const int half_size = this->uom_save.size() / 2;
                    for(int mm=0; mm<size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat]+mm] = occ_mat[iat][target_l][0][0].c[mm];
                        this->uom_save[half_size + eff_pot_pw_index[iat]+mm] = occ_mat[iat][target_l][0][1].c[mm];
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "copy_occ_mat");
}


void Plus_U_Base::zero_occ_mat(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "zero_occ_mat");
    ModuleBase::timer::start("Plus_U_Base", "zero_occ_mat");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            for (int l = 0; l < ucell.atoms[T].nwl + 1; l++)
            {
                const int N = ucell.atoms[T].l_nchi[l];

                for (int n = 0; n < N; n++)
                {
                    if (this->nspin == 4)
                    {
                        occ_mat[iat][l][n][0].zero_out();
                    }
                    else if (this->nspin == 1 || this->nspin == 2)
                    {
                        occ_mat[iat][l][n][0].zero_out();
                        occ_mat[iat][l][n][1].zero_out();
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "zero_occ_mat");
}


void Plus_U_Base::mix_occ_mat(const UnitCell& ucell,
                              const double& mixing_beta)
{
    ModuleBase::TITLE("Plus_U_Base", "mix_occ_mat");
    ModuleBase::timer::start("Plus_U_Base", "mix_occ_mat");

    double beta = mixing_beta;

    for (int T = 0; T < ucell.ntype; T++)
    {
        int target_l = get_orbital_corr(T);
        if (target_l == -1)
            continue;

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);

            if (this->nspin == 4)
            {
                const int size = occ_mat[iat][target_l][0][0].nr * occ_mat[iat][target_l][0][0].nc;
                for (int mm = 0; mm < size; mm++)
                {
                    occ_mat[iat][target_l][0][0].c[mm] = occ_mat[iat][target_l][0][0].c[mm] * beta + occ_mat_save[iat][target_l][0][0].c[mm] * (1.0 - beta);
                }
                if (this->uom_save.size() != 0)
                {
                    for (int mm = 0; mm < size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat] + mm] = occ_mat[iat][target_l][0][0].c[mm];
                    }
                }
            }
            else if (this->nspin == 1 || this->nspin == 2)
            {
                const int size = occ_mat[iat][target_l][0][0].nr * occ_mat[iat][target_l][0][0].nc;
                const int half_size = this->uom_save.size() / 2;
                for (int mm = 0; mm < size; mm++)
                {
                    occ_mat[iat][target_l][0][0].c[mm] = occ_mat[iat][target_l][0][0].c[mm] * beta + occ_mat_save[iat][target_l][0][0].c[mm] * (1.0 - beta);
                    occ_mat[iat][target_l][0][1].c[mm] = occ_mat[iat][target_l][0][1].c[mm] * beta + occ_mat_save[iat][target_l][0][1].c[mm] * (1.0 - beta);
                }
                if (this->uom_save.size() != 0)
                {
                    for (int mm = 0; mm < size; mm++)
                    {
                        this->uom_save[eff_pot_pw_index[iat] + mm] = occ_mat[iat][target_l][0][0].c[mm];
                        this->uom_save[half_size + eff_pot_pw_index[iat] + mm] = occ_mat[iat][target_l][0][1].c[mm];
                    }
                }
            }
        }
    }
    ModuleBase::timer::end("Plus_U_Base", "mix_occ_mat");
}


void Plus_U_Base::set_occ_mat(const UnitCell& ucell)
{
    ModuleBase::TITLE("Plus_U_Base", "set_occ_mat");
    ModuleBase::timer::start("Plus_U_Base", "set_occ_mat");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T)) continue;
        const int l = get_orbital_corr(T);
        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            if (this->nspin == 4)
            {
                for(int mm = 0; mm < occ_mat[iat][l][0][0].nr * occ_mat[iat][l][0][0].nc; mm++)
                    occ_mat[iat][l][0][0].c[mm] = this->uom_array[eff_pot_pw_index[iat] + mm];
            }
            else if (this->nspin == 1 || this->nspin == 2)
            {
                const int half_size = this->uom_array.size() / 2;
                for(int mm = 0; mm < occ_mat[iat][l][0][0].nr * occ_mat[iat][l][0][0].nc; mm++)
                {
                    occ_mat[iat][l][0][0].c[mm] = this->uom_array[eff_pot_pw_index[iat] + mm];
                    if (this->nspin == 2)
                    {
                        occ_mat[iat][l][0][1].c[mm] = this->uom_array[half_size + eff_pot_pw_index[iat] + mm];
                    }
                }
            }
        }
    }

    ModuleBase::timer::end("Plus_U_Base", "set_occ_mat");
}


void Plus_U_Base::get_occ_mat_flat(const int iat, const int l, std::vector<double>& occ) const
{
    const int tlp1 = 2 * l + 1;
    const int size = tlp1 * tlp1;
    if (nspin == 2)
    {
        for (int is = 0; is < 2; is++)
        {
            for (int i = 0; i < size; i++)
            {
                occ[is * size + i] = occ_mat[iat][l][0][is].c[i];
            }
        }
    }
    else
    {
        for (int i = 0; i < static_cast<int>(occ.size()); i++)
        {
            occ[i] = occ_mat[iat][l][0][0].c[i];
        }
    }
}


void Plus_U_Base::set_occ_mat_flat(const int iat, const int l, const int spin,
                                   const std::vector<double>& occ)
{
    for (int i = 0; i < static_cast<int>(occ.size()); i++)
    {
        occ_mat[iat][l][0][spin].c[i] = occ[i];
    }
}


void Plus_U_Base::read_occup_m(const UnitCell& ucell,
                               const std::string& fn,
                               const std::string& init_chg,
                               int nspin,
                               int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "read_occup_m");

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    std::ifstream ifdftu(fn.c_str(), std::ios::in);

    if (!ifdftu)
    {
        if (occ_mat_ctrl > 0)
        {
            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not find the file dm_onsite_ini.txt. Please check your dm_onsite_ini.txt");
        }
        else
        {
            if (init_chg == "file")
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not find the file dm_onsite.txt. Please do scf calculation first");
            }
        }
        ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "Can not open dm_onsite.txt file");
    }

    ifdftu.clear();
    ifdftu.seekg(0);

    char word[20];

    int T = 0;
    int iat = 0;
    int spin = 0;
    int L = 0;
    int zeta = 0;

    ifdftu.rdstate();

    while (ifdftu.good())
    {
        ifdftu >> word;
        if (ifdftu.eof())
        {
            break;
        }

        if (strcmp("Atom=", word) == 0)
        {
            ifdftu >> iat;
            iat -= 1;
            ifdftu >> word;

            if (strcmp("L=", word) != 0)
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> L;
            ifdftu >> word;

            if (strcmp("ORBITAL=", word) != 0)
            {
                ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
            }
            ifdftu >> zeta;
            ifdftu.ignore(150, '\n');

            T = ucell.iat2it[iat];
            const int NL = ucell.atoms[T].nwl + 1;
            const int LC = get_orbital_corr(T);

            for (int l = 0; l < NL; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                if (nspin == 1 || nspin == 2)
                {
                    for (int is = 0; is < 2; is++)
                    {
                        ifdftu >> word;
                        if (strcmp("spin=", word) == 0)
                        {
                            ifdftu >> spin;
                            spin -= 1;
                            ifdftu.ignore(150, '\n');

                            double value = 0.0;
                            for (int m0 = 0; m0 < 2 * L + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                                {
                                    ifdftu >> value;
                                    occ_mat[iat][L][zeta][spin](m0, m1) = value;
                                }
                                ifdftu.ignore(150, '\n');
                            }
                        }
                        else
                        {
                            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
                        }
                    }
                }
                else if (nspin == 4) // SOC
                {
                    double value = 0.0;
                    for (int m0 = 0; m0 < 2 * L + 1; m0++)
                    {
                        for (int ipol0 = 0; ipol0 < npol; ipol0++)
                        {
                            const int m0_all = m0 + (2 * L + 1) * ipol0;

                            for (int m1 = 0; m1 < 2 * L + 1; m1++)
                            {
                                for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                {
                                    int m1_all = m1 + (2 * L + 1) * ipol1;
                                    ifdftu >> value;
                                    occ_mat[iat][L][zeta][0](m0_all, m1_all) = value;
                                }
                            }
                            ifdftu.ignore(150, '\n');
                        }
                    }
                }
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("Plus_U_Base::read_occup_m", "WRONG IN READING LOCAL OCCUPATION NUMBER MATRIX FROM Plus_U FILE");
        }

        ifdftu.rdstate();

        if (ifdftu.eof() != 0)
        {
            break;
        }
    }

    return;
}


void Plus_U_Base::local_occup_bcast(const UnitCell& ucell,
                                    int nspin,
                                    int npol)
{
    ModuleBase::TITLE("Plus_U_Base", "local_occup_bcast");

    for (int T = 0; T < ucell.ntype; T++)
    {
        if (!has_correlated_orbital(T))
        {
            continue;
        }

        for (int I = 0; I < ucell.atoms[T].na; I++)
        {
            const int iat = ucell.itia2iat(T, I);
            const int L = get_orbital_corr(T);

            for (int l = 0; l <= ucell.atoms[T].nwl; l++)
            {
                if (l != get_orbital_corr(T))
                {
                    continue;
                }

                for (int n = 0; n < ucell.atoms[T].l_nchi[l]; n++)
                {
                    if (n != 0)
                    {
                        continue;
                    }

                    if (nspin == 1 || nspin == 2)
                    {
                        for (int spin = 0; spin < 2; spin++)
                        {
                            for (int m0 = 0; m0 < 2 * l + 1; m0++)
                            {
                                for (int m1 = 0; m1 < 2 * l + 1; m1++)
                                {
#ifdef __MPI
                                    MPI_Bcast(&occ_mat[iat][l][n][spin](m0, m1), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
                                }
                            }
                        }
                    }
                    else if (nspin == 4) // SOC
                    {
                        for (int m0 = 0; m0 < 2 * L + 1; m0++)
                        {
                            for (int ipol0 = 0; ipol0 < npol; ipol0++)
                            {
                                const int m0_all = m0 + (2 * L + 1) * ipol0;

                                for (int m1 = 0; m1 < 2 * L + 1; m1++)
                                {
                                    for (int ipol1 = 0; ipol1 < npol; ipol1++)
                                    {
                                        int m1_all = m1 + (2 * L + 1) * ipol1;
#ifdef __MPI
                                        MPI_Bcast(&occ_mat[iat][l][n][0](m0_all, m1_all),
                                                  1,
                                                  MPI_DOUBLE,
                                                  0,
                                                  MPI_COMM_WORLD);
#endif
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}


// cal_occ_pw() is implemented in source_lcao/module_dftu/dftu_pw.cpp
// as a Plus_U_Base method. It will be relocated to this directory
// in Phase 5 of the class-split refactor.

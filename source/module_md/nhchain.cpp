#include "nhchain.h"

#include "md_func.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include "module_base/timer.h"

Nose_Hoover::Nose_Hoover(MD_para& MD_para_in, UnitCell& unit_in) : MD_base(MD_para_in, unit_in)
{
    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    mdp.md_pfirst /= unit_transform;
    mdp.md_plast /= unit_transform;
    mdp.md_pfreq *= ModuleBase::AU_to_FS;

    if (mdp.md_tfirst == 0)
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", " md_tfirst must be larger than 0 in NHC");
    }

    /// init NPT related variables
    for (int i = 0; i < 6; ++i)
    {
        pstart[i] = pstop[i] = pfreq[i] = p_target[i] = pflag[i] = 0;
    }

    if (mdp.md_type == "npt")
    {
        /// determine the NPT methods
        if (mdp.md_pmode == "iso")
        {
            mdp.md_pcouple = "xyz";
            pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
            pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
            pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
            pflag[0] = pflag[1] = pflag[2] = 1;
        }
        else if (mdp.md_pmode == "aniso")
        {
            if (mdp.md_pcouple == "xyz")
            {
                ModuleBase::WARNING_QUIT("Nose_Hoover", "md_pcouple==xyz will convert aniso to iso!");
            }
            pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
            pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
            pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
            pflag[0] = pflag[1] = pflag[2] = 1;
        }
        /**
         * The lattice must be lower-triangular under tri mode.
         * e11  0    0
         * e21  e22  0
         * e31  e32  e33
         * Under Voigt notation, xx, yy, zz, yz, xz, xy.
         */
        else if (mdp.md_pmode == "tri")
        {
            if (ucell.latvec.e12 || ucell.latvec.e13 || ucell.latvec.e23)
            {
                ModuleBase::WARNING_QUIT("Nose_Hoover", "the lattice must be lower-triangular when md_pmode == tri!");
            }
            pstart[0] = pstart[1] = pstart[2] = mdp.md_pfirst;
            pstop[0] = pstop[1] = pstop[2] = mdp.md_plast;
            pfreq[0] = pfreq[1] = pfreq[2] = mdp.md_pfreq;
            pflag[0] = pflag[1] = pflag[2] = 1;

            pstart[3] = pstart[4] = pstart[5] = 0;
            pstop[3] = pstop[4] = pstop[5] = 0;
            pfreq[3] = pfreq[4] = pfreq[5] = mdp.md_pfreq;
            pflag[3] = pflag[4] = pflag[5] = 1;
        }
        else
        {
            ModuleBase::WARNING_QUIT("Nose_Hoover", "No such md_pmode yet!");
        }
    }

    /// determine whether NPT ensemble
    npt_flag = 0;
    for (int i = 0; i < 6; ++i)
    {
        npt_flag += pflag[i];
    }
    pdim = pflag[0] + pflag[1] + pflag[2];

    tdof = 3 * ucell.nat - frozen_freedom_;

    /// allocate thermostats coupled with particles
    mass_eta = new double[mdp.md_tchain];
    eta = new double[mdp.md_tchain];
    v_eta = new double[mdp.md_tchain + 1];
    g_eta = new double[mdp.md_tchain];

    v_eta[mdp.md_tchain] = 0;
    for (int i = 0; i < mdp.md_tchain; ++i)
    {
        eta[i] = v_eta[i] = g_eta[i] = 0;
    }

    /// allocate barostat and thermostats coupled with barostat
    if (npt_flag)
    {
        for (int i = 0; i < 6; ++i)
        {
            v_omega[i] = mass_omega[i] = 0;
        }

        if (mdp.md_pchain)
        {
            mass_peta = new double[mdp.md_pchain];
            peta = new double[mdp.md_pchain];
            v_peta = new double[mdp.md_pchain + 1];
            g_peta = new double[mdp.md_pchain];

            v_peta[mdp.md_pchain] = 0;
            for (int i = 0; i < mdp.md_pchain; ++i)
            {
                peta[i] = v_peta[i] = g_peta[i] = 0;
            }
        }
    }

    /// w[0] = 1;

    w[0] = 0.784513610477560;
    w[6] = 0.784513610477560;
    w[1] = 0.235573213359357;
    w[5] = 0.235573213359357;
    w[2] = -1.17767998417887;
    w[4] = -1.17767998417887;
    w[3] = 1 - w[0] - w[1] - w[2] - w[4] - w[5] - w[6];
}

Nose_Hoover::~Nose_Hoover()
{
    delete[] mass_eta;
    delete[] eta;
    delete[] v_eta;
    delete[] g_eta;

    if (npt_flag && mdp.md_pchain)
    {
        delete[] mass_peta;
        delete[] peta;
        delete[] v_peta;
        delete[] g_peta;
    }
}

void Nose_Hoover::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Nose_Hoover", "setup");
    ModuleBase::timer::tick("Nose_Hoover", "setup");

    MD_base::setup(p_esolver, global_readin_dir);
    if (mdp.md_type == "npt")
    {
        ucell.cell_parameter_updated = true;
    }

    /// determine target temperature
    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, mdp.md_tfirst, mdp.md_tlast);

    /// init thermostats coupled with particles
    mass_eta[0] = tdof * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for (int m = 1; m < mdp.md_tchain; ++m)
    {
        mass_eta[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
        g_eta[m] = (mass_eta[m - 1] * v_eta[m - 1] * v_eta[m - 1] - t_target) / mass_eta[m];
    }

    /// NPT ensemble
    if (npt_flag)
    {
        /// determine target stress
        target_stress();

        /// couple stress component due to md_pcouple
        couple_stress();

        /// init barostat
        double nkt = (ucell.nat + 1) * t_target;

        for (int i = 0; i < 6; ++i)
        {
            if (pflag[i])
            {
                mass_omega[i] = nkt / pfreq[i] / pfreq[i];
            }
        }

        /// init thermostats coupled with barostat
        if (mdp.md_pchain)
        {
            mass_peta[0] = t_target / mdp.md_pfreq / mdp.md_pfreq;
            for (int m = 1; m < mdp.md_pchain; ++m)
            {
                mass_peta[m] = t_target / mdp.md_pfreq / mdp.md_pfreq;
                g_peta[m] = (mass_peta[m - 1] * v_peta[m - 1] * v_peta[m - 1] - t_target) / mass_peta[m];
            }
        }
    }

    ModuleBase::timer::tick("Nose_Hoover", "setup");
}

void Nose_Hoover::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("Nose_Hoover", "first_half");
    ModuleBase::timer::tick("Nose_Hoover", "first_half");

    /// update thermostats coupled with barostat if NPT ensemble
    if (npt_flag && mdp.md_pchain)
    {
        baro_thermo();
    }

    /// update target T
    t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, mdp.md_tfirst, mdp.md_tlast);

    /// update thermostats coupled with particles
    particle_thermo();

    if (npt_flag)
    {
        /// update temperature and stress due to velocity rescaling
        t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);
        MD_func::compute_stress(ucell, vel, allmass, mdp.cal_stress,virial, stress);

        /// couple stress component due to md_pcouple
        couple_stress();

        /// determine target stress
        target_stress();

        /// update v_omega
        update_baro();

        /// update vel due to barostat
        vel_baro();
    }

    /// perform half-step update of vel due to atomic force
    MD_base::update_vel(force);

    if (npt_flag)
    {
        /// perform half-step update of volume
        update_volume(ofs);
    }

    /// perform one step update of pos due to atomic velocity
    MD_base::update_pos();

    if (npt_flag)
    {
        /// perform half-step update of volume
        update_volume(ofs);
    }

    ModuleBase::timer::tick("Nose_Hoover", "first_half");
}

void Nose_Hoover::second_half()
{
    ModuleBase::TITLE("Nose_Hoover", "second_half");
    ModuleBase::timer::tick("Nose_Hoover", "second_half");

    /// perform half-step update of vel due to atomic force
    MD_base::update_vel(force);

    if (npt_flag)
    {
        /// update vel due to barostat
        vel_baro();
    }

    /// update temperature and kinetic energy due to velocity rescaling
    t_current = MD_func::current_temp(kinetic, ucell.nat, frozen_freedom_, allmass, vel);

    if (npt_flag)
    {
        /// update stress due to velocity rescaling
        MD_func::compute_stress(ucell, vel, allmass, mdp.cal_stress,virial, stress);

        /// couple stress component due to md_pcouple
        couple_stress();

        /// update v_omega
        update_baro();
    }

    /// update thermostats coupled with particles
    particle_thermo();

    /// update thermostats coupled with barostat if NPT ensemble
    if (npt_flag && mdp.md_pchain)
    {
        baro_thermo();
    }

    ModuleBase::timer::tick("Nose_Hoover", "second_half");
}

void Nose_Hoover::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
}

void Nose_Hoover::write_restart(const std::string& global_out_dir)
{
    if (!mdp.my_rank)
    {
        std::stringstream ssc;
        ssc << global_out_dir << "Restart_md.dat";
        std::ofstream file(ssc.str().c_str());

        file << step_ + step_rst_ << std::endl;
        file << mdp.md_tchain << std::endl;
        for (int i = 0; i < mdp.md_tchain; ++i)
        {
            file << eta[i] << "   ";
        }
        file << std::endl;
        for (int i = 0; i < mdp.md_tchain; ++i)
        {
            file << v_eta[i] << "   ";
        }
        file << std::endl;

        /// npt
        if (npt_flag)
        {
            for (int i = 0; i < 6; ++i)
            {
                file << v_omega[i] << "   ";
            }
            file << std::endl;

            file << mdp.md_pchain << std::endl;
            for (int i = 0; i < mdp.md_pchain; ++i)
            {
                file << peta[i] << "   ";
            }
            file << std::endl;
            for (int i = 0; i < mdp.md_pchain; ++i)
            {
                file << v_peta[i] << "   ";
            }
            file << std::endl;
        }
        file.close();
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Nose_Hoover::restart(const std::string& global_readin_dir)
{
    bool ok = true;
    bool ok2 = true;
    bool ok3 = true;

    if (!mdp.my_rank)
    {
        std::stringstream ssc;
        ssc << global_readin_dir << "Restart_md.dat";
        std::ifstream file(ssc.str().c_str());

        if (!file)
        {
            ok = false;
        }

        if (ok)
        {
            double Mnum;
            file >> step_rst_ >> Mnum;

            if (Mnum != mdp.md_tchain)
            {
                ok2 = false;
            }

            if (ok2)
            {
                for (int i = 0; i < mdp.md_tchain; ++i)
                {
                    file >> eta[i];
                }
                for (int i = 0; i < mdp.md_tchain; ++i)
                {
                    file >> v_eta[i];
                }
            }

            /// npt
            if (npt_flag)
            {
                for (int i = 0; i < 6; ++i)
                {
                    file >> v_omega[i];
                }

                file >> Mnum;
                if (Mnum != mdp.md_pchain)
                {
                    ok3 = false;
                }

                if (ok3)
                {
                    for (int i = 0; i < mdp.md_pchain; ++i)
                    {
                        file >> peta[i];
                    }
                    for (int i = 0; i < mdp.md_pchain; ++i)
                    {
                        file >> v_peta[i];
                    }
                }
            }

            file.close();
        }
    }

#ifdef __MPI
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ok2, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ok3, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if (!ok)
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "no Restart_md.dat !");
    }
    if (!ok2)
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Num of thermostats coupled with particles is not the same !");
    }
    if (!ok3)
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Num of thermostats coupled with barostat is not the same !");
    }

#ifdef __MPI
    MPI_Bcast(&step_rst_, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(eta, mdp.md_tchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v_eta, mdp.md_tchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (npt_flag)
    {
        MPI_Bcast(v_omega, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(peta, mdp.md_pchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(v_peta, mdp.md_pchain, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
#endif
}

void Nose_Hoover::particle_thermo()
{
    /// update mass_eta
    mass_eta[0] = tdof * t_target / mdp.md_tfreq / mdp.md_tfreq;
    for (int m = 1; m < mdp.md_tchain; ++m)
    {
        mass_eta[m] = t_target / mdp.md_tfreq / mdp.md_tfreq;
    }

    /// propogate g_eta
    if (mass_eta[0] > 0)
    {
        g_eta[0] = (2 * kinetic - tdof * t_target) / mass_eta[0];
    }
    else
    {
        g_eta[0] = 0;
    }

    /// integrate loop
    double factor;
    double scale = 1.0;
    double KE = kinetic;
    for (int i = 0; i < nc_tchain; ++i)
    {
        for (int j = 0; j < nys; ++j)
        {
            double delta = w[j] * mdp.md_dt / nc_tchain;

            /// propogate v_eta
            for (int m = mdp.md_tchain - 1; m >= 0; --m)
            {
                factor = exp(-v_eta[m + 1] * delta / 8.0);
                v_eta[m] *= factor;
                v_eta[m] += g_eta[m] * delta / 4.0;
                v_eta[m] *= factor;
            }

            /// propogate eta
            for (int m = 0; m < mdp.md_tchain; ++m)
            {
                eta[m] += v_eta[m] * delta / 2.0;
            }

            /// update rescale factor of particle velocity
            scale *= exp(-v_eta[0] * delta / 2.0);
            if (!std::isfinite(scale))
            {
                ModuleBase::WARNING_QUIT("Nose_Hoover", "Please set a proper md_tfreq!");
            }
            KE = kinetic * scale * scale;

            /// propogate g_eta
            if (mass_eta[0] > 0)
            {
                g_eta[0] = (2 * KE - tdof * t_target) / mass_eta[0];
            }
            else
            {
                g_eta[0] = 0;
            }

            /// propogate v_eta
            v_eta[0] *= factor;
            v_eta[0] += g_eta[0] * delta / 4.0;
            v_eta[0] *= factor;

            for (int m = 1; m < mdp.md_tchain; ++m)
            {
                factor = exp(-v_eta[m + 1] * delta / 8.0);
                v_eta[m] *= factor;
                g_eta[m] = (mass_eta[m - 1] * v_eta[m - 1] * v_eta[m - 1] - t_target) / mass_eta[m];
                v_eta[m] += g_eta[m] * delta / 4.0;
                v_eta[m] *= factor;
            }
        }
    }

    /// rescale velocity due to thermostats
    for (int i = 0; i < ucell.nat; ++i)
    {
        vel[i] *= scale;
    }
}

void Nose_Hoover::baro_thermo()
{
    /// the freedom of lattice
    int pdof = npt_flag;

    /// update kenetic energy of lattice
    double ke_omega = 0;
    for (int i = 0; i < 6; ++i)
    {
        if (pflag[i])
        {
            ke_omega += mass_omega[i] * v_omega[i] * v_omega[i];
        }
    }

    /// update force
    double lkt_press = t_target;
    if (mdp.md_pmode != "iso")
    {
        lkt_press *= pdof;
    }
    g_peta[0] = (ke_omega - lkt_press) / mass_peta[0];

    /// integrate loop
    double factor;
    double scale = 1.0;
    double kecurrent = ke_omega;
    for (int i = 0; i < nc_pchain; ++i)
    {
        for (int j = 0; j < nys; ++j)
        {
            double delta = w[j] * mdp.md_dt / nc_pchain;

            /// propogate v_peta
            for (int m = mdp.md_pchain - 1; m >= 0; --m)
            {
                factor = exp(-v_peta[m + 1] * delta / 8.0);
                v_peta[m] *= factor;
                v_peta[m] += g_peta[m] * delta / 4.0;
                v_peta[m] *= factor;
            }

            /// propogate peta
            for (int m = 0; m < mdp.md_pchain; ++m)
            {
                peta[m] += v_peta[m] * delta / 2.0;
            }

            /// update rescale factor of lattice velocity
            scale *= exp(-v_peta[0] * delta / 2.0);
            kecurrent = ke_omega * scale * scale;

            /// propogate g_peta
            g_peta[0] = (kecurrent - lkt_press) / mass_peta[0];

            /// propogate v_peta
            v_peta[0] *= factor;
            v_peta[0] += g_peta[0] * delta / 4.0;
            v_peta[0] *= factor;

            for (int m = 1; m < mdp.md_pchain; ++m)
            {
                factor = exp(-v_peta[m + 1] * delta / 8.0);
                v_peta[m] *= factor;
                g_peta[m] = (mass_peta[m - 1] * v_peta[m - 1] * v_peta[m - 1] - t_target) / mass_peta[m];
                v_peta[m] += g_eta[m] * delta / 4.0;
                v_peta[m] *= factor;
            }
        }
    }

    /// rescale lattice due to thermostats
    for (int i = 0; i < 6; ++i)
    {
        if (pflag[i])
        {
            v_omega[i] *= scale;
        }
    }
}

void Nose_Hoover::update_baro()
{
    double term_one = 0;
    if (mdp.md_pmode == "iso")
    {
        term_one = tdof * t_current;
    }
    else
    {
        ModuleBase::matrix t_vector;
        MD_func::temp_vector(ucell.nat, vel, allmass, t_vector);

        for (int i = 0; i < 3; ++i)
        {
            if (pflag[i])
            {
                term_one += t_vector(i, i);
            }
        }
    }
    term_one /= pdim * ucell.nat;

    double g_omega;
    double term_two = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (pflag[i])
        {
            g_omega = (p_current[i] - p_hydro) * ucell.omega / mass_omega[i] + term_one / mass_omega[i];
            v_omega[i] += g_omega * mdp.md_dt / 2.0;
            term_two += v_omega[i];
        }
    }
    term_two /= pdim * ucell.nat;

    for (int i = 3; i < 6; ++i)
    {
        if (pflag[i])
        {
            g_omega = p_current[i] * ucell.omega / mass_omega[i];
            v_omega[i] += g_omega * mdp.md_dt / 2.0;
        }
    }

    mtk_term = term_two;
}

void Nose_Hoover::vel_baro()
{
    double factor[3];
    for (int i = 0; i < 3; ++i)
    {
        factor[i] = exp(-(v_omega[i] + mtk_term) * mdp.md_dt / 4);
    }

    for (int i = 0; i < ucell.nat; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            vel[i][j] *= factor[j];
        }

        /// Note: I am not sure whether fixed atoms should update here
        if (ionmbl[i][0])
        {
            vel[i][0] -= (vel[i][1] * v_omega[5] + vel[i][2] * v_omega[4]) * mdp.md_dt / 2;
        }
        if (ionmbl[i][1])
        {
            vel[i][1] -= vel[i][2] * v_omega[3] * mdp.md_dt / 2;
        }

        for (int j = 0; j < 3; ++j)
        {
            vel[i][j] *= factor[j];
        }
    }
}

void Nose_Hoover::update_volume(std::ofstream& ofs)
{
    double factor;

    /// tri mode, off-diagonal components, first half
    if (pflag[4])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 16);
        ucell.latvec.e31 *= factor;
        ucell.latvec.e31 += (v_omega[5] * ucell.latvec.e32 + v_omega[4] * ucell.latvec.e33);
        ucell.latvec.e31 *= factor;
    }

    if (pflag[3])
    {
        factor = exp(v_omega[1] * mdp.md_dt / 8);
        ucell.latvec.e32 *= factor;
        ucell.latvec.e32 += (v_omega[3] * ucell.latvec.e33);
        ucell.latvec.e32 *= factor;
    }

    if (pflag[5])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 8);
        ucell.latvec.e21 *= factor;
        ucell.latvec.e21 += (v_omega[5] * ucell.latvec.e22);
        ucell.latvec.e21 *= factor;
    }

    if (pflag[4])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 16);
        ucell.latvec.e31 *= factor;
        ucell.latvec.e31 += (v_omega[5] * ucell.latvec.e32 + v_omega[4] * ucell.latvec.e33);
        ucell.latvec.e31 *= factor;
    }

    /// Diagonal components
    if (pflag[0])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 2);
        ucell.latvec.e11 *= factor;
    }

    if (pflag[1])
    {
        factor = exp(v_omega[1] * mdp.md_dt / 2);
        ucell.latvec.e22 *= factor;
    }

    if (pflag[2])
    {
        factor = exp(v_omega[2] * mdp.md_dt / 2);
        ucell.latvec.e33 *= factor;
    }

    /// tri mode, off-diagonal components, second half
    if (pflag[4])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 16);
        ucell.latvec.e31 *= factor;
        ucell.latvec.e31 += (v_omega[5] * ucell.latvec.e32 + v_omega[4] * ucell.latvec.e33);
        ucell.latvec.e31 *= factor;
    }

    if (pflag[3])
    {
        factor = exp(v_omega[1] * mdp.md_dt / 8);
        ucell.latvec.e32 *= factor;
        ucell.latvec.e32 += (v_omega[3] * ucell.latvec.e33);
        ucell.latvec.e32 *= factor;
    }

    if (pflag[5])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 8);
        ucell.latvec.e21 *= factor;
        ucell.latvec.e21 += (v_omega[5] * ucell.latvec.e22);
        ucell.latvec.e21 *= factor;
    }

    if (pflag[4])
    {
        factor = exp(v_omega[0] * mdp.md_dt / 16);
        ucell.latvec.e31 *= factor;
        ucell.latvec.e31 += (v_omega[5] * ucell.latvec.e32 + v_omega[4] * ucell.latvec.e33);
        ucell.latvec.e31 *= factor;
    }

    /// reset ucell and pos due to change of lattice
    ucell.setup_cell_after_vc(ofs);
}

void Nose_Hoover::target_stress()
{
    double delta = static_cast<double>(step_ + step_rst_) / mdp.md_nstep;

    p_hydro = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (pflag[i])
        {
            p_target[i] = pstart[i] + delta * (pstop[i] - pstart[i]);
            p_hydro += p_target[i];
        }
    }
    if (pdim)
    {
        p_hydro /= pdim;
    }

    for (int i = 3; i < 6; ++i)
    {
        if (pflag[i])
        {
            p_target[i] = pstart[i] + delta * (pstop[i] - pstart[i]);
        }
    }
}

void Nose_Hoover::couple_stress()
{
    if (mdp.md_pcouple == "xyz")
    {
        double ave = (stress(0, 0) + stress(1, 1) + stress(2, 2)) / 3.0;
        p_current[0] = p_current[1] = p_current[2] = ave;
    }
    else if (mdp.md_pcouple == "xy")
    {
        double ave = (stress(0, 0) + stress(1, 1)) / 2.0;
        p_current[0] = p_current[1] = ave;
        p_current[2] = stress(2, 2);
    }
    else if (mdp.md_pcouple == "yz")
    {
        double ave = (stress(1, 1) + stress(2, 2)) / 2.0;
        p_current[1] = p_current[2] = ave;
        p_current[0] = stress(0, 0);
    }
    else if (mdp.md_pcouple == "xz")
    {
        double ave = (stress(0, 0) + stress(2, 2)) / 2.0;
        p_current[0] = p_current[2] = ave;
        p_current[1] = stress(1, 1);
    }
    else
    {
        p_current[0] = stress(0, 0);
        p_current[1] = stress(1, 1);
        p_current[2] = stress(2, 2);
    }

    if (!std::isfinite(p_current[0]) || !std::isfinite(p_current[1]) || !std::isfinite(p_current[2]))
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Non-numeric stress component!");
    }

    p_current[3] = stress(1, 2);
    p_current[4] = stress(0, 2);
    p_current[5] = stress(0, 1);

    if (!std::isfinite(p_current[3]) || !std::isfinite(p_current[4]) || !std::isfinite(p_current[5]))
    {
        ModuleBase::WARNING_QUIT("Nose_Hoover", "Non-numeric stress component!");
    }
}
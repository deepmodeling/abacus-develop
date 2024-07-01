#include "module_io/io_dmk.h"

#include "module_base/parallel_common.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"

void ModuleIO::read_dmk(
#ifdef __MPI
    const int nnrg,
    const int* trace_lo,
#endif
    const bool gamma_only_local,
    const int nlocal,
    const int nspin,
    const int& is,
    const std::string& fn,
    double*** DM,
    double** DM_R,
    double& ef,
    const UnitCell* ucell)
{
    ModuleBase::TITLE("ModuleIO", "read_dmk");
    ModuleBase::timer::tick("ModuleIO", "read_dmk");

    GlobalV::ofs_running << "\n processor 0 is reading density matrix from file < " << fn << " > " << std::endl;
    // xiaohui modify 2015-03-25
    // bool quit_mesia = false;
    bool quit_abacus = false;

    std::ifstream ifs;
    if (GlobalV::MY_RANK == 0)
    {
        ifs.open(fn.c_str());
        if (!ifs)
        {
            // xiaohui modify 2015-03-25
            // quit_mesia = true;
            quit_abacus = true;
        }
        else
        {
            // if the number is not match,
            // quit the program or not.
            bool quit = false;

            std::string name;
            ifs >> name;

            // check lattice constant, unit is Angstrom
            ModuleBase::CHECK_DOUBLE(ifs, ucell->lat0 * ModuleBase::BOHR_TO_A, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e11, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e12, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e13, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e21, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e22, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e23, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e31, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e32, quit);
            ModuleBase::CHECK_DOUBLE(ifs, ucell->latvec.e33, quit);

            for (int it = 0; it < ucell->ntype; it++)
            {
                ModuleBase::CHECK_STRING(ifs, ucell->atoms[it].label, quit);
            }

            for (int it = 0; it < ucell->ntype; it++)
            {
                ModuleBase::CHECK_DOUBLE(ifs, ucell->atoms[it].na, quit);
            }

            std::string coordinate;
            ifs >> coordinate;

            for (int it = 0; it < ucell->ntype; it++)
            {
                for (int ia = 0; ia < ucell->atoms[it].na; ia++)
                {
                    ModuleBase::CHECK_DOUBLE(ifs, ucell->atoms[it].taud[ia].x, quit);
                    ModuleBase::CHECK_DOUBLE(ifs, ucell->atoms[it].taud[ia].y, quit);
                    ModuleBase::CHECK_DOUBLE(ifs, ucell->atoms[it].taud[ia].z, quit);
                }
            }

            ModuleBase::CHECK_INT(ifs, nspin);
            ModuleBase::GlobalFunc::READ_VALUE(ifs, ef);
            ModuleBase::CHECK_INT(ifs, nlocal);
            ModuleBase::CHECK_INT(ifs, nlocal);
        } // If file exist, read in data.
    }     // Finish reading the first part of density matrix.

#ifndef __MPI
    GlobalV::ofs_running << " Read SPIN = " << is + 1 << " density matrix now." << std::endl;

    if (gamma_only_local)
    {
        for (int i = 0; i < nlocal; ++i)
        {
            for (int j = 0; j < nlocal; ++j)
            {
                ifs >> DM[is][i][j];
            }
        }
    }
    else
    {
#ifdef __MPI
        ModuleBase::WARNING_QUIT("ModuleIO::read_dmk", "The nnrg should not be update");
        ModuleBase::CHECK_INT(ifs, nnrg);

        for (int i = 0; i < nnrg; ++i)
        {
            ifs >> DM_R[is][i];
        }
#endif
    }
#else

    // distribution of necessary data
    // xiaohui modify 2015-03-25
    // Parallel_Common::bcast_bool(quit_mesia);
    Parallel_Common::bcast_bool(quit_abacus);
    // xiaohui modify 2015-03-25
    // if(quit_mesia)
    if (quit_abacus)
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_dmk", "Can not find the density matrix file.");
    }

    Parallel_Common::bcast_double(ef);

    if (gamma_only_local)
    {
        std::vector<double> tmp(nlocal);
        for (int i = 0; i < nlocal; ++i)
        {
            // GlobalV::ofs_running << " i=" << i << std::endl;
            ModuleBase::GlobalFunc::ZEROS(tmp, nlocal);
            if (GlobalV::MY_RANK == 0)
            {
                for (int j = 0; j < nlocal; ++j)
                {
                    ifs >> tmp[j];
                }
            }
            Parallel_Common::bcast_double(tmp.data(), nlocal);

            const int mu = trace_lo[i];
            if (mu >= 0)
            {
                for (int j = 0; j < nlocal; ++j)
                {
                    const int nu = trace_lo[j];
                    if (nu >= 0)
                    {
                        DM[is][mu][nu] = tmp[j];
                    }
                }
            }
        } // i
    }
    else
    {
        ModuleBase::WARNING_QUIT("ModuleIO::read_dmk", "not ready to readin DM_R");
    }
#endif
    if (GlobalV::MY_RANK == 0)
        ifs.close();

    GlobalV::ofs_running << " Finish reading density matrix." << std::endl;

    ModuleBase::timer::tick("ModuleIO", "read_dmk");
    return;
}

std::string ModuleIO::dmk_gen_fname(const bool gamma_only, const int ispin, const int ik)
{
    if (gamma_only)
    {
        return "SPIN" + std::to_string(ispin + 1) + "_DM";
    }
    else
    {
        // this case is not implemented now.
        ModuleBase::WARNING_QUIT("dmk_gen_fname", "Not implemented for non-gamma_only case.");
    }
}

template <typename T>
void ModuleIO::write_dmk(const std::vector<std::vector<T>> dmk,
                         const int precision,
                         const std::vector<double> efs,
                         const UnitCell* ucell,
                         const Parallel_2D& pv)
{
    ModuleBase::TITLE("ModuleIO", "write_dmk");
    ModuleBase::timer::tick("ModuleIO", "write_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm_2D, &my_rank);
#endif

    bool gamma_only = std::is_same<double, T>::value;
    int nlocal = pv.get_global_row_size();
    int nspin = efs.size();
    int nk = dmk.size() / nspin;
    if (nk * nspin != dmk.size())
    {
        ModuleBase::WARNING_QUIT("write_dmk", "The size of dmk is not consistent with nspin and nk.");
    }
    Parallel_2D pv_glb;

    // when nspin == 2, assume the order of K in dmk is K1_up, K2_up, ..., K1_down, K2_down, ...
    for (int ispin = 0; ispin < nspin; ispin++)
    {
        for (int ik = 0; ik < nk; ik++)
        {
            // gather dmk[ik] to dmk_global
            std::vector<T> dmk_global(my_rank == 0 ? nlocal * nlocal : 0);
#ifdef __MPI
            pv_glb.set(nlocal, nlocal, nlocal, pv.comm_2D, pv.blacs_ctxt);
            Cpxgemr2d(nlocal,
                      nlocal,
                      const_cast<T*>(dmk[ik + nk * ispin].data()),
                      1,
                      1,
                      const_cast<int*>(pv.desc),
                      dmk_global.data(),
                      1,
                      1,
                      pv_glb.desc,
                      pv_glb.blacs_ctxt);
#else
            dmk_global = dmk[ik + nk * ispin];
#endif

            if (my_rank == 0)
            {
                std::string fn = GlobalV::global_out_dir + dmk_gen_fname(gamma_only, ispin, ik);
                std::ofstream ofs(fn.c_str());

                if (!ofs)
                {
                    ModuleBase::WARNING("ModuleIO::write_dmk", "Can't create DENSITY MATRIX File!");
                }

                // write the UnitCell information
                ofs << ucell->latName << std::endl;
                ofs << " " << ucell->lat0 * ModuleBase::BOHR_TO_A << std::endl;
                ofs << " " << ucell->latvec.e11 << " " << ucell->latvec.e12 << " " << ucell->latvec.e13 << std::endl;
                ofs << " " << ucell->latvec.e21 << " " << ucell->latvec.e22 << " " << ucell->latvec.e23 << std::endl;
                ofs << " " << ucell->latvec.e31 << " " << ucell->latvec.e32 << " " << ucell->latvec.e33 << std::endl;
                for (int it = 0; it < ucell->ntype; it++)
                {
                    ofs << " " << ucell->atoms[it].label;
                }
                ofs << std::endl;
                for (int it = 0; it < ucell->ntype; it++)
                {
                    ofs << " " << ucell->atoms[it].na;
                }
                ofs << std::endl;
                ofs << "Direct" << std::endl;
                for (int it = 0; it < ucell->ntype; it++)
                {
                    Atom* atom = &ucell->atoms[it];
                    ofs << std::setprecision(15);
                    for (int ia = 0; ia < ucell->atoms[it].na; ia++)
                    {
                        ofs << " " << atom->taud[ia].x << " " << atom->taud[ia].y << " " << atom->taud[ia].z
                            << std::endl;
                    }
                }
                ofs << "\n " << dmk.size(); // nspin
                ofs << "\n " << efs[ispin] << " (fermi energy)";
                ofs << "\n  " << nlocal << " " << nlocal << std::endl;

                ofs << std::setprecision(precision);
                ofs << std::scientific;
                for (int i = 0; i < nlocal; ++i)
                {
                    for (int j = 0; j < nlocal; ++j)
                    {
                        if (j % 8 == 0)
                        {
                            ofs << "\n";
                        }
                        if (std::is_same<double, T>::value)
                        {
                            ofs << " " << dmk_global[i * nlocal + j];
                        }
                        else if (std::is_same<std::complex<double>, T>::value)
                        {
                            ofs << " (" << std::real(dmk_global[i * nlocal + j]) << ","
                                << std::imag(dmk_global[i * nlocal + j]) << ")";
                        }
                    }
                }
                ofs.close();
            } // rank0
        }     // ik
    }         // ispin

    ModuleBase::timer::tick("ModuleIO", "write_dmk");
}

template void ModuleIO::write_dmk<double>(const std::vector<std::vector<double>> dmk,
                                          const int precision,
                                          const std::vector<double> efs,
                                          const UnitCell* ucell,
                                          const Parallel_2D& pv);

template void ModuleIO::write_dmk<std::complex<double>>(const std::vector<std::vector<std::complex<double>>> dmk,
                                                        const int precision,
                                                        const std::vector<double> efs,
                                                        const UnitCell* ucell,
                                                        const Parallel_2D& pv);
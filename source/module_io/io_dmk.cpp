#include "module_io/io_dmk.h"

#include "module_parameter/parameter.h"
#include "module_base/parallel_common.h"
#include "module_base/scalapack_connector.h"
#include "module_base/timer.h"

/*
The format of the DMK file is as follows:
'''
<latName>
<lat0>
<latvec a>
<latvec b>
<latvec c>
<label1> <label2> ...
<na1> <na2> ...
Direct
<label1-atom1-x> <label1-atom1-y> <label1-atom1-z>
<label1-atom2-x> <label1-atom2-y> <label1-atom2-z>
...
<label2-atom1-x> <label2-atom1-y> <label2-atom1-z>
<label2-atom2-x> <label2-atom2-y> <label2-atom2-z>
...

<ispin>
<Efermi> (fermi energy)
<nlocal> <nlocal>

<dmk data>
...
'''


Example:
'''
sc
 5.29177
 1 0 0
 0 1 0
 0 0 1
 H
 2
Direct
 0 0 0.933859999999186
 0 0 0.0661400000008143

 1
 -0.0883978533958687 (fermi energy)
  10 10

 5.773e-01 3.902e-02 1.661e-02 4.797e-17 -2.255e-17 5.773e-01 3.902e-02
-1.661e-02 -1.461e-17 -4.414e-17
 ...
 '''
 */

std::string ModuleIO::dmk_gen_fname(const bool gamma_only,
                                    const int ispin,
                                    const int ik) {
    if (gamma_only) {
        return "SPIN" + std::to_string(ispin + 1) + "_DM";
    } else {
        // this case is not implemented now.
        ModuleBase::WARNING_QUIT("dmk_gen_fname",
                                 "Not implemented for non-gamma_only case.");
    }
}

void ModuleIO::dmk_write_ucell(std::ofstream& ofs, const UnitCell* ucell) {
    // write the UnitCell information
    ofs << ucell->latName << std::endl;
    ofs << " " << ucell->lat0 * ModuleBase::BOHR_TO_A << std::endl;
    ofs << " " << ucell->latvec.e11 << " " << ucell->latvec.e12 << " "
        << ucell->latvec.e13 << std::endl;
    ofs << " " << ucell->latvec.e21 << " " << ucell->latvec.e22 << " "
        << ucell->latvec.e23 << std::endl;
    ofs << " " << ucell->latvec.e31 << " " << ucell->latvec.e32 << " "
        << ucell->latvec.e33 << std::endl;
    for (int it = 0; it < ucell->ntype; it++) {
        ofs << " " << ucell->atoms[it].label;
    }
    ofs << std::endl;
    for (int it = 0; it < ucell->ntype; it++) {
        ofs << " " << ucell->atoms[it].na;
    }
    ofs << std::endl;
    ofs << "Direct" << std::endl;
    for (int it = 0; it < ucell->ntype; it++) {
        Atom* atom = &ucell->atoms[it];
        ofs << std::setprecision(15);
        for (int ia = 0; ia < ucell->atoms[it].na; ia++) {
            ofs << " " << atom->taud[ia].x << " " << atom->taud[ia].y << " "
                << atom->taud[ia].z << std::endl;
        }
    }
}

void ModuleIO::dmk_read_ucell(std::ifstream& ifs) {
    std::string tmp;
    for (int i = 0; i < 6; i++) {
        std::getline(ifs, tmp); // latName + lat0 + latvec + atom label
    }
    std::getline(ifs, tmp); // atom number of each type

    std::istringstream iss(tmp);
    int natom = 0;
    int total_natom = 0;
    while (iss >> natom) {
        total_natom += natom;
    }
    for (int i = 0; i < total_natom + 1; i++) {
        std::getline(ifs, tmp); // Direct + atom coordinates
    }
}

void ModuleIO::dmk_readData(std::ifstream& ifs, double& data) { ifs >> data; }

void ModuleIO::dmk_readData(std::ifstream& ifs, std::complex<double>& data) {
    double real, imag;
    ifs >> real;
    ifs >> imag;
    data = std::complex<double>(real, imag);
}

template <typename T>
bool ModuleIO::read_dmk(const int nspin,
                        const int nk,
                        const Parallel_2D& pv,
                        const std::string& dmk_dir,
                        std::vector<std::vector<T>>& dmk) {
    ModuleBase::TITLE("ModuleIO", "read_dmk");
    ModuleBase::timer::tick("ModuleIO", "read_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    int nlocal = pv.get_global_row_size();
    bool gamma_only = std::is_same<double, T>::value;
    std::vector<std::vector<T>> dmk_global(nspin * nk, std::vector<T>(nlocal * nlocal,0));

    // write a lambda function to check the consistency of the data
    auto check_consistency = [&](const std::string& fn,
                                 const std::string& name,
                                 const std::string& value,
                                 const int& target) {
        if (std::stoi(value) != target) {
            ModuleBase::WARNING("ModuleIO::read_dmk",
                                name + " is not consistent in file < " + fn
                                    + " >.");
            std::cout << name << " = " << target << ", " << name
                      << " in file = " << value << std::endl;
            return false;
        }
        return true;
    };

    bool read_success = true;
    std::string tmp;
    if (my_rank == 0) {
        for (int ispin = 0; ispin < nspin; ispin++) {
            for (int ik = 0; ik < nk; ik++) {
                std::string fn = dmk_dir + dmk_gen_fname(gamma_only, ispin, ik);
                std::ifstream ifs(fn.c_str());

                if (!ifs) {
                    ModuleBase::WARNING("ModuleIO::read_dmk",
                                        "Can't open DENSITY MATRIX File < " + fn
                                            + " >.");
                    read_success = false;
                    break;
                }

                // read the UnitCell
                dmk_read_ucell(ifs);

                ifs >> tmp; // nspin
                if (!check_consistency(fn, "ispin", tmp, ispin + 1)) {
                    read_success = false;
                    ifs.close();
                    break;
                }
                ifs >> tmp;
                ifs >> tmp;
                ifs >> tmp; // fermi energy
                ifs >> tmp; // nlocal
                if (!check_consistency(fn, "nlocal", tmp, nlocal)) {
                    read_success = false;
                    ifs.close();
                    break;
                }
                ifs >> tmp; // nlocal
                if (!check_consistency(fn, "nlocal", tmp, nlocal)) {
                    read_success = false;
                    ifs.close();
                    break;
                }

                // read the DMK data
                for (int i = 0; i < nlocal; ++i) {
                    for (int j = 0; j < nlocal; ++j) {
                        dmk_readData(
                            ifs,
                            dmk_global[ik + nk * ispin][i * nlocal + j]);
                    }
                }
                ifs.close();
            } // ik
            if (!read_success) {
                break;
            }
        } // ispin
    }     // rank0

#ifdef __MPI
    MPI_Bcast(&read_success, 1, MPI_C_BOOL, 0, pv.comm());
#endif

    if (read_success) {
#ifdef __MPI
        // seperate dmk data to each processor with 2D block distribution
        dmk.resize(nspin * nk,
                   std::vector<T>(pv.get_row_size() * pv.get_col_size()));
        Parallel_2D pv_glb;
        pv_glb.set(nlocal, nlocal, nlocal, pv.blacs_ctxt);
        for (int ik = 0; ik < nspin * nk; ik++) {
            Cpxgemr2d(nlocal,
                      nlocal,
                      dmk_global[ik].data(),
                      1,
                      1,
                      pv_glb.desc,
                      dmk[ik].data(),
                      1,
                      1,
                      const_cast<int*>(pv.desc),
                      pv_glb.blacs_ctxt);
        }
#else
        dmk = dmk_global;
#endif
    }
    ModuleBase::timer::tick("ModuleIO", "read_dmk");
    return read_success;
}

template <typename T>
void ModuleIO::write_dmk(const std::vector<std::vector<T>>& dmk,
                         const int precision,
                         const std::vector<double>& efs,
                         const UnitCell* ucell,
                         const Parallel_2D& pv) {
    ModuleBase::TITLE("ModuleIO", "write_dmk");
    ModuleBase::timer::tick("ModuleIO", "write_dmk");

    int my_rank = 0;
#ifdef __MPI
    MPI_Comm_rank(pv.comm(), &my_rank);
#endif

    bool gamma_only = std::is_same<double, T>::value;
    int nlocal = pv.get_global_row_size();
    int nspin = efs.size();
    int nk = dmk.size() / nspin;
    if (nk * nspin != dmk.size()) {
        ModuleBase::WARNING_QUIT(
            "write_dmk",
            "The size of dmk is not consistent with nspin and nk.");
    }
    Parallel_2D pv_glb;

    // when nspin == 2, assume the order of K in dmk is K1_up, K2_up, ...,
    // K1_down, K2_down, ...
    for (int ispin = 0; ispin < nspin; ispin++) {
        for (int ik = 0; ik < nk; ik++) {
            // gather dmk[ik] to dmk_global
            std::vector<T> dmk_global(my_rank == 0 ? nlocal * nlocal : 0);
#ifdef __MPI
            pv_glb.set(nlocal, nlocal, nlocal, pv.blacs_ctxt);
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

            if (my_rank == 0) {
                std::string fn = PARAM.globalv.global_out_dir
                                 + dmk_gen_fname(gamma_only, ispin, ik);
                std::ofstream ofs(fn.c_str());

                if (!ofs) {
                    ModuleBase::WARNING("ModuleIO::write_dmk",
                                        "Can't create DENSITY MATRIX File < "
                                            + fn + " >.");
                    continue;
                }

                // write the UnitCell information
                dmk_write_ucell(ofs, ucell);

                ofs << "\n " << dmk.size(); // nspin
                ofs << "\n " << std::fixed << std::setprecision(5) << efs[ispin]
                    << " (fermi energy)";
                ofs << "\n  " << nlocal << " " << nlocal << std::endl;

                ofs << std::setprecision(precision);
                ofs << std::scientific;
                for (int i = 0; i < nlocal; ++i) {
                    for (int j = 0; j < nlocal; ++j) {
                        if (j % 8 == 0) {
                            ofs << "\n";
                        }
                        if (std::is_same<double, T>::value) {
                            ofs << " " << dmk_global[i * nlocal + j];
                        } else if (std::is_same<std::complex<double>,
                                                T>::value) {
                            ofs << " (" << std::real(dmk_global[i * nlocal + j])
                                << "," << std::imag(dmk_global[i * nlocal + j])
                                << ")";
                        }
                    }
                }
                ofs.close();
            } // rank0
        }     // ik
    }         // ispin

    ModuleBase::timer::tick("ModuleIO", "write_dmk");
}

template bool ModuleIO::read_dmk<double>(const int nspin,
                                         const int nk,
                                         const Parallel_2D& pv,
                                         const std::string& dmk_dir,
                                         std::vector<std::vector<double>>& dmk);

template bool ModuleIO::read_dmk<std::complex<double>>(
    const int nspin,
    const int nk,
    const Parallel_2D& pv,
    const std::string& dmk_dir,
    std::vector<std::vector<std::complex<double>>>& dmk);

template void
    ModuleIO::write_dmk<double>(const std::vector<std::vector<double>>& dmk,
                                const int precision,
                                const std::vector<double>& efs,
                                const UnitCell* ucell,
                                const Parallel_2D& pv);

template void ModuleIO::write_dmk<std::complex<double>>(
    const std::vector<std::vector<std::complex<double>>>& dmk,
    const int precision,
    const std::vector<double>& efs,
    const UnitCell* ucell,
    const Parallel_2D& pv);

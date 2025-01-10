#include "elecstate.h"
#include "elecstate_getters.h"
#include "module_parameter/parameter.h"
#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "occupy.h"
namespace elecstate
{
/**
 * Notes on refactor of ESolver's functions
 *
 * the print of SCF iteration on-the-fly information.
 * 1. Previously it is expected for nspin 1, 2, and 4, also with xc_type 3/5 or not, the information will organized in
 * different ways. This brings inconsistencies between patterns of print and make it hard to vectorize information.
 * 2. the function print_etot actually do two kinds of things, 1) print information into running_*.log, 2) print
 * information onto screen. These two tasks are, in no way should be placed/implemented in one function directly
 * 3. there are information redundance: the istep of SCF can provide information determing whether print out the SCF
 * iteration info. table header or not, rather than dividing into two functions and hard code the format.
 *
 * For nspin 1, print: ITER, ETOT, EDIFF, DRHO, TIME
 * 	   nspin 2, print: ITER, TMAG, AMAG, ETOT, EDIFF, DRHO, TIME
 * 	   nspin 4 with nlcc, print: ITER, TMAGX, TMAGY, TMAGZ, AMAG, ETOT, EDIFF, DRHO, TIME
 * xc type_id 3/5: DKIN
 *
 * Based on summary above, there are several groups of info:
 * 1. counting: ITER
 * 2. (optional) magnetization: TMAG or TMAGX-TMAGY-TMAGZ, AMAG
 * 3. energies: ETOT, EDIFF
 * 4. densities: DRHO, DKIN(optional)
 * 5. time: TIME
 */
void print_scf_iterinfo(const std::string& ks_solver,
                        const int& istep,
                        const int& witer,
                        const std::vector<double>& mag,
                        const int& wmag,
                        const double& etot,
                        const double& ediff,
                        const int& wener,
                        const std::vector<double>& drho,
                        const int& wrho,
                        const double& time,
                        const int& wtime)
{
    std::map<std::string, std::string> iter_header_dict
        = {{"cg", "CG"},
           {"cg_in_lcao", "CG"},
           {"lapack", "LA"},
           {"genelpa", "GE"},
           {"elpa", "EL"},
           {"dav", "DA"},
           {"dav_subspace", "DS"},
           {"scalapack_gvx", "GV"},
           {"cusolver", "CU"},
           {"bpcg", "BP"},
           {"pexsi", "PE"}}; // I change the key of "cg_in_lcao" to "CG" because all the other are only two letters
    // ITER column
    std::vector<std::string> th_fmt = {" %-" + std::to_string(witer) + "s"}; // table header: th: ITER
    std::vector<std::string> td_fmt
        = {" " + iter_header_dict[ks_solver] + "%-" + std::to_string(witer - 2) + ".0f"}; // table data: td: GE10086
    // magnetization column, might be non-exist, but size of mag can only be 0, 2 or 4
    for (int i = 0; i < mag.size(); i++)
    {
        th_fmt.emplace_back(" %" + std::to_string(wmag) + "s");
    }
    for (int i = 0; i < mag.size(); i++)
    {
        td_fmt.emplace_back(" %" + std::to_string(wmag) + ".2e");
    } // hard-code precision here
    // energies
    for (int i = 0; i < 2; i++)
    {
        th_fmt.emplace_back(" %" + std::to_string(wener) + "s");
    }
    for (int i = 0; i < 2; i++)
    {
        td_fmt.emplace_back(" %" + std::to_string(wener) + ".8e");
    }
    // densities column, size can be 1 or 2, DRHO or DRHO, DKIN
    for (int i = 0; i < drho.size(); i++)
    {
        th_fmt.emplace_back(" %" + std::to_string(wrho) + "s");
    }
    for (int i = 0; i < drho.size(); i++)
    {
        td_fmt.emplace_back(" %" + std::to_string(wrho) + ".4e");
    }
    // time column, trivial
    th_fmt.emplace_back(" %" + std::to_string(wtime) + "s\n");
    td_fmt.emplace_back(" %" + std::to_string(wtime) + ".2f\n");
    // contents
    std::vector<std::string> titles;
    std::vector<double> values;
    switch (mag.size())
    {
    case 2:
        titles = {"ITER",
                  FmtCore::center("TMAG", wmag),
                  FmtCore::center("AMAG", wmag),
                  FmtCore::center("ETOT/eV", wener),
                  FmtCore::center("EDIFF/eV", wener),
                  FmtCore::center("DRHO", wrho)};
        values = {double(istep), mag[0], mag[1], etot, ediff, drho[0]};
        break;
    case 4:
        titles = {"ITER",
                  FmtCore::center("TMAGX", wmag),
                  FmtCore::center("TMAGY", wmag),
                  FmtCore::center("TMAGZ", wmag),
                  FmtCore::center("AMAG", wmag),
                  FmtCore::center("ETOT/eV", wener),
                  FmtCore::center("EDIFF/eV", wener),
                  FmtCore::center("DRHO", wrho)};
        values = {double(istep), mag[0], mag[1], mag[2], mag[3], etot, ediff, drho[0]};
        break;
    default:
        titles = {"ITER",
                  FmtCore::center("ETOT/eV", wener),
                  FmtCore::center("EDIFF/eV", wener),
                  FmtCore::center("DRHO", wrho)};
        values = {double(istep), etot, ediff, drho[0]};
        break;
    }
    if (drho.size() > 1)
    {
        titles.push_back(FmtCore::center("DKIN", wrho));
        values.push_back(drho[1]);
    }
    titles.push_back(FmtCore::center("TIME/s", wtime));
    values.push_back(time);
    std::string buf;
    if (istep == 1)
    {
        for (int i = 0; i < titles.size(); i++)
        {
            buf += FmtCore::format(th_fmt[i].c_str(), titles[i]);
        }
    }
    for (int i = 0; i < values.size(); i++)
    {
        buf += FmtCore::format(td_fmt[i].c_str(), values[i]);
    }
    std::cout << buf;
}
/// @brief print and check for band energy and occupations
/// @param ofs
void ElecState::print_eigenvalue(std::ofstream& ofs)
{
    bool wrong = false;
    for (int ik = 0; ik < this->klist->get_nks(); ++ik)
    {
        for (int ib = 0; ib < this->ekb.nc; ++ib)
        {
            if (std::abs(this->ekb(ik, ib)) > 1.0e10)
            {
                GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << this->ekb(ik, ib) << " Ry"
                                     << std::endl;
                wrong = true;
            }
        }
    }
    if (wrong)
    {
        ModuleBase::WARNING_QUIT("print_eigenvalue", "Eigenvalues are too large!");
    }

    if (GlobalV::MY_RANK != 0)
    {
        return;
    }

    ModuleBase::TITLE("ESolver_KS_PW", "print_eigenvalue");

    ofs << "\n STATE ENERGY(eV) AND OCCUPATIONS ";
    for (int ik = 0; ik < this->klist->get_nks(); ik++)
    {
        ofs << std::setprecision(5);
        ofs << std::setiosflags(std::ios::showpoint);
        if (ik == 0)
        {
            ofs << "   NSPIN == " << PARAM.inp.nspin << std::endl;
            if (PARAM.inp.nspin == 2)
            {
                ofs << "SPIN UP : " << std::endl;
            }
        }
        else if (ik == this->klist->get_nks() / 2)
        {
            if (PARAM.inp.nspin == 2)
            {
                ofs << "SPIN DOWN : " << std::endl;
            }
        }

        if (PARAM.inp.nspin == 2)
        {
            if (this->klist->isk[ik] == 0)
            {
                ofs << " " << ik + 1 << "/" << this->klist->get_nks() / 2
                    << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x << " " << this->klist->kvec_c[ik].y << " "
                    << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);
            }
            if (this->klist->isk[ik] == 1)
            {
                ofs << " " << ik + 1 - this->klist->get_nks() / 2 << "/" << this->klist->get_nks() / 2
                    << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x << " " << this->klist->kvec_c[ik].y << " "
                    << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);
            }
        } // Pengfei Li  added  14-9-9
        else
        {
            ofs << " " << ik + 1 << "/" << this->klist->get_nks()
                << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x << " " << this->klist->kvec_c[ik].y << " "
                << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik] << " pws)" << std::endl;

            ofs << std::setprecision(6);
        }

        ofs << std::setprecision(6);
        ofs << std::setiosflags(std::ios::showpoint);
        for (int ib = 0; ib < this->ekb.nc; ib++)
        {
            ofs << std::setw(8) << ib + 1 << std::setw(15) << this->ekb(ik, ib) * ModuleBase::Ry_to_eV << std::setw(15)
                << this->wg(ik, ib) << std::endl;
        }
        ofs << std::endl;
    } // end ik
    return;
}

/// @brief function for printing eigenvalues : ekb
/// @param ik: index of kpoints
/// @param printe: print energy every 'printe' electron iteration.
/// @param iter: index of iterations
void ElecState::print_band(const int& ik, const int& printe, const int& iter)
{
    // check the band energy.
    bool wrong = false;
    for (int ib = 0; ib < PARAM.inp.nbands; ++ib)
    {
        if (std::abs(this->ekb(ik, ib)) > 1.0e10)
        {
            GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << this->ekb(ik, ib) << " Ry"
                                 << std::endl;
            wrong = true;
        }
    }
    if (wrong)
    {
        ModuleBase::WARNING_QUIT("print_eigenvalue", "Eigenvalues are too large!");
    }

    if (GlobalV::MY_RANK == 0)
    {
        if (printe > 0 && ((iter + 1) % printe == 0))
        {
            GlobalV::ofs_running << std::setprecision(6);
            GlobalV::ofs_running << " Energy (eV) & Occupations  for spin=" << this->klist->isk[ik] + 1
                                 << " K-point=" << ik + 1 << std::endl;
            GlobalV::ofs_running << std::setiosflags(std::ios::showpoint);
            for (int ib = 0; ib < PARAM.inp.nbands; ib++)
            {
                GlobalV::ofs_running << " " << std::setw(6) << ib + 1 << std::setw(15)
                                     << this->ekb(ik, ib) * ModuleBase::Ry_to_eV;
                // for the first electron iteration, we don't have the energy
                // spectrum, so we can't get the occupations.
                GlobalV::ofs_running << std::setw(15) << this->wg(ik, ib);
                GlobalV::ofs_running << std::endl;
            }
        }
    }
    return;
}

/// @brief print total free energy and other energies
/// @param ucell: unit cell
/// @param converged: if converged
/// @param iter_in: iter
/// @param scf_thr: threshold for scf
/// @param duration: time of each iteration
/// @param pw_diag_thr: threshold for diagonalization
/// @param avg_iter: averaged diagonalization iteration of each scf iteration
/// @param print: if print to screen
void ElecState::print_etot(const Magnetism& magnet,
                           const bool converged,
                           const int& iter_in,
                           const double& scf_thr,
                           const double& scf_thr_kin,
                           const double& duration,
                           const int printe,
                           const double& pw_diag_thr,
                           const double& avg_iter,
                           const bool print)
{
    ModuleBase::TITLE("energy", "print_etot");
    const int iter = iter_in;
    const int nrxx = this->charge->nrxx;
    const int nxyz = this->charge->nxyz;

    GlobalV::ofs_running << std::setprecision(12);
    GlobalV::ofs_running << std::setiosflags(std::ios::right);

    GlobalV::ofs_running << "\n Density error is " << scf_thr << std::endl;

    if (PARAM.inp.basis_type == "pw") {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Error Threshold", pw_diag_thr); // xiaohui add 2013-09-02
}

    std::vector<std::string> titles;
    std::vector<double> energies_Ry;
    std::vector<double> energies_eV;
    if (printe > 0 && ((iter + 1) % printe == 0 || converged || iter == PARAM.inp.scf_nmax))
    {
        int n_order = std::max(0, Occupy::gaussian_type);
        titles.push_back("E_KohnSham");
        energies_Ry.push_back(this->f_en.etot);
        titles.push_back("E_KS(sigma->0)");
        energies_Ry.push_back(this->f_en.etot - this->f_en.demet / (2 + n_order));
        titles.push_back("E_Harris");
        energies_Ry.push_back(this->f_en.etot_harris);
        titles.push_back("E_band");
        energies_Ry.push_back(this->f_en.eband);
        titles.push_back("E_one_elec");
        energies_Ry.push_back(this->f_en.eband + this->f_en.deband);
        titles.push_back("E_Hartree");
        energies_Ry.push_back(this->f_en.hartree_energy);
        titles.push_back("E_xc");
        energies_Ry.push_back(this->f_en.etxc - this->f_en.etxcc);
        titles.push_back("E_Ewald");
        energies_Ry.push_back(this->f_en.ewald_energy);
        titles.push_back("E_entropy(-TS)");
        energies_Ry.push_back(this->f_en.demet);
        titles.push_back("E_descf");
        energies_Ry.push_back(this->f_en.descf);
        titles.push_back("E_LocalPP");
        energies_Ry.push_back(this->f_en.e_local_pp);
        std::string vdw_method = PARAM.inp.vdw_method;
        if (vdw_method == "d2") // Peize Lin add 2014-04, update 2021-03-09
        {
            titles.push_back("E_vdwD2");
            energies_Ry.push_back(this->f_en.evdw);
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj") // jiyy add 2019-05, update 2021-05-02
        {
            titles.push_back("E_vdwD3");
            energies_Ry.push_back(this->f_en.evdw);
        }
        titles.push_back("E_exx");
        energies_Ry.push_back(this->f_en.exx);
        if (PARAM.inp.imp_sol)
        {
            titles.push_back("E_sol_el");
            energies_Ry.push_back(this->f_en.esol_el);
            titles.push_back("E_sol_cav");
            energies_Ry.push_back(this->f_en.esol_cav);
        }
        if (PARAM.inp.efield_flag)
        {
            titles.push_back("E_efield");
            energies_Ry.push_back(elecstate::Efield::etotefield);
        }
        if (PARAM.inp.gate_flag)
        {
            titles.push_back("E_gatefield");
            energies_Ry.push_back(elecstate::Gatefield::etotgatefield);
        }

#ifdef __DEEPKS
        if (PARAM.inp.deepks_scf) // caoyu add 2021-08-10
        {
            titles.push_back("E_DeePKS");
            energies_Ry.push_back(GlobalC::ld.E_delta);
        }
#endif
    }
    else
    {
        titles.push_back("E_KohnSham");
        energies_Ry.push_back(this->f_en.etot);
        titles.push_back("E_Harris");
        energies_Ry.push_back(this->f_en.etot_harris);
    }

    if (PARAM.globalv.two_fermi)
    {
        titles.push_back("E_Fermi_up");
        energies_Ry.push_back(this->eferm.ef_up);
        titles.push_back("E_Fermi_dw");
        energies_Ry.push_back(this->eferm.ef_dw);
    }
    else
    {
        titles.push_back("E_Fermi");
        energies_Ry.push_back(this->eferm.ef);
    }
    if (PARAM.inp.out_bandgap)
    {
        if (!PARAM.globalv.two_fermi)
        {
            titles.push_back("E_bandgap");
            energies_Ry.push_back(this->bandgap);
        }
        else
        {
            titles.push_back("E_bandgap_up");
            energies_Ry.push_back(this->bandgap_up);
            titles.push_back("E_bandgap_dw");
            energies_Ry.push_back(this->bandgap_dw);
        }
    }
    energies_eV.resize(energies_Ry.size());
    std::transform(energies_Ry.begin(), energies_Ry.end(), energies_eV.begin(), [](double ener) {
        return ener * ModuleBase::Ry_to_eV;
    });
    FmtTable table({"Energy", "Rydberg", "eV"},
                   titles.size(),
                   {"%-14s", "%20.10f", "%20.10f"},
                   {FmtTable::Align::LEFT, FmtTable::Align::CENTER});
    table << titles << energies_Ry << energies_eV;
    GlobalV::ofs_running << table.str() << std::endl;
    if (PARAM.inp.out_level == "ie" || PARAM.inp.out_level == "m") // xiaohui add 'm' option, 2015-09-16
    {
        std::vector<double> mag;
        switch (PARAM.inp.nspin)
        {
        case 2:
            mag = {magnet.tot_magnetization, magnet.abs_magnetization};
            break;
        case 4:
            mag = {magnet.tot_magnetization_nc[0],
                   magnet.tot_magnetization_nc[1],
                   magnet.tot_magnetization_nc[2],
                   magnet.abs_magnetization};
            break;
        default:
            mag = {};
            break;
        }
        std::vector<double> drho = {scf_thr};
        if (elecstate::get_xc_func_type() == 3 || elecstate::get_xc_func_type() == 5)
        {
            drho.push_back(scf_thr_kin);
        }
        elecstate::print_scf_iterinfo(PARAM.inp.ks_solver,
                                      iter,
                                      6,
                                      mag,
                                      10,
                                      this->f_en.etot * ModuleBase::Ry_to_eV,
                                      this->f_en.etot_delta * ModuleBase::Ry_to_eV,
                                      16,
                                      drho,
                                      12,
                                      duration,
                                      6);
    }
    return;
}

/// @brief function to print name, value and value*Ry_to_eV
/// @param name: name
/// @param value: value
void ElecState::print_format(const std::string& name, const double& value)
{
    GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
    std::stringstream name2;
    name2 << name;
    GlobalV::ofs_running << " " << std::setw(16) << name2.str() << std::setw(30) << value << std::setw(30)
                         << value * ModuleBase::Ry_to_eV << std::endl;
    GlobalV::ofs_running << std::resetiosflags(std::ios::showpos);
    return;
}
} // namespace elecstate
#include "elecstate.h"
#include "module_base/formatter.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_common.h"
#include "module_elecstate/module_pot/H_Hartree_pw.h"
#include "module_elecstate/module_pot/efield.h"
#include "module_elecstate/module_pot/gatefield.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
#include "module_parameter/parameter.h"
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
           {"pexsi", "PE"},
           {"cusolvermp", "CM"}}; // I change the key of "cg_in_lcao" to "CG" because all the other are only two letters
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

/// @brief function for printing eigenvalues : ekb
/// @param ik: index of kpoints
/// @param printe: print energy every 'printe' electron iteration.
/// @param iter: index of iterations
void print_band(const ModuleBase::matrix& ekb,
                const ModuleBase::matrix& wg,
                const K_Vectors* klist,
                const int& ik,
                const int& printe,
                const int& iter,
                std::ofstream &ofs)
{
    const double largest_eig = 1.0e10;

    // check the band energy.
    bool wrong = false;
    for (int ib = 0; ib < PARAM.globalv.nbands_l; ++ib)
    {
        if (std::abs(ekb(ik, ib)) > largest_eig)
        {
            GlobalV::ofs_warning << " ik=" << ik + 1 << " ib=" << ib + 1 << " " << ekb(ik, ib) << " Ry" << std::endl;
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
            ofs << std::setprecision(6);
            ofs << " Energy (eV) & Occupations for spin=" << klist->isk[ik] + 1
                                 << " k-point=" << ik + 1 << std::endl;
            ofs << std::setiosflags(std::ios::showpoint);
            for (int ib = 0; ib < PARAM.globalv.nbands_l; ib++)
            {
                ofs << " " << std::setw(6) << ib + 1 << std::setw(15)
                                     << ekb(ik, ib) * ModuleBase::Ry_to_eV;
                // for the first electron iteration, we don't have the energy
                // spectrum, so we can't get the occupations.
                ofs << std::setw(15) << wg(ik, ib);
                ofs << std::endl;
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
void print_etot(const Magnetism& magnet,
                const ElecState& elec,
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
    const int nrxx = elec.charge->nrxx;
    const int nxyz = elec.charge->nxyz;

    GlobalV::ofs_running << std::setprecision(12);
    GlobalV::ofs_running << std::setiosflags(std::ios::right);
    GlobalV::ofs_running << " Electron density deviation is " << scf_thr << std::endl;

    if (PARAM.inp.basis_type == "pw")
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Diago Threshold", pw_diag_thr);
    }

    std::vector<std::string> titles;
    std::vector<double> energies_Ry;
    std::vector<double> energies_eV;

    if (printe > 0 && ((iter + 1) % printe == 0 || converged || iter == PARAM.inp.scf_nmax))
    {
        int n_order = std::max(0, Occupy::gaussian_type);
        titles.push_back("E_KohnSham");
        energies_Ry.push_back(elec.f_en.etot);
        titles.push_back("E_KS(sigma->0)");
        energies_Ry.push_back(elec.f_en.etot - elec.f_en.demet / (2 + n_order));
        titles.push_back("E_Harris");
        energies_Ry.push_back(elec.f_en.etot_harris);
        titles.push_back("E_band");
        energies_Ry.push_back(elec.f_en.eband);
        titles.push_back("E_one_elec");
        energies_Ry.push_back(elec.f_en.eband + elec.f_en.deband);
        titles.push_back("E_Hartree");
        energies_Ry.push_back(elec.f_en.hartree_energy);
        titles.push_back("E_xc");
        energies_Ry.push_back(elec.f_en.etxc - elec.f_en.etxcc);
        titles.push_back("E_Ewald");
        energies_Ry.push_back(elec.f_en.ewald_energy);
        titles.push_back("E_entropy(-TS)");
        energies_Ry.push_back(elec.f_en.demet);
        titles.push_back("E_descf");
        energies_Ry.push_back(elec.f_en.descf);
        titles.push_back("E_LocalPP");
        energies_Ry.push_back(elec.f_en.e_local_pp);
        std::string vdw_method = PARAM.inp.vdw_method;
        if (vdw_method == "d2") // Peize Lin add 2014-04, update 2021-03-09
        {
            titles.push_back("E_vdwD2");
            energies_Ry.push_back(elec.f_en.evdw);
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj") // jiyy add 2019-05, update 2021-05-02
        {
            titles.push_back("E_vdwD3");
            energies_Ry.push_back(elec.f_en.evdw);
        }
        titles.push_back("E_exx");
        energies_Ry.push_back(elec.f_en.exx);
        if (PARAM.inp.imp_sol)
        {
            titles.push_back("E_sol_el");
            energies_Ry.push_back(elec.f_en.esol_el);
            titles.push_back("E_sol_cav");
            energies_Ry.push_back(elec.f_en.esol_cav);
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
        if (PARAM.inp.deepks_scf)
        {
            titles.push_back("E_DeePKS");
            energies_Ry.push_back(elec.f_en.edeepks_delta);
        }
#endif
    }
    else
    {
        titles.push_back("E_KohnSham");
        energies_Ry.push_back(elec.f_en.etot);
        titles.push_back("E_Harris");
        energies_Ry.push_back(elec.f_en.etot_harris);
    }

    // print out the Fermi energy if needed
    if (PARAM.globalv.two_fermi)
    {
        titles.push_back("E_Fermi_up");
        energies_Ry.push_back(elec.eferm.ef_up);
        titles.push_back("E_Fermi_dw");
        energies_Ry.push_back(elec.eferm.ef_dw);
    }
    else
    {
        titles.push_back("E_Fermi");
        energies_Ry.push_back(elec.eferm.ef);
    }

    // print out the band gap if needed
    if (PARAM.inp.out_bandgap)
    {
        if (!PARAM.globalv.two_fermi)
        {
            titles.push_back("E_bandgap");
            energies_Ry.push_back(elec.bandgap);
        }
        else
        {
            titles.push_back("E_bandgap_up");
            energies_Ry.push_back(elec.bandgap_up);
            titles.push_back("E_bandgap_dw");
            energies_Ry.push_back(elec.bandgap_dw);
        }
    }
    energies_eV.resize(energies_Ry.size());
    std::transform(energies_Ry.begin(), energies_Ry.end(), energies_eV.begin(), [](double ener) {
        return ener * ModuleBase::Ry_to_eV;
    });

    // for each SCF step, we print out energy
    FmtTable table(/*titles=*/{"Energy", "Rydberg", "eV"},
                   /*nrows=*/titles.size(),
                   /*formats=*/{"%-14s", "%20.10f", "%20.10f"}, 
                   /*indents=*/0,
                   /*align=*/{/*value*/FmtTable::Align::LEFT, /*title*/FmtTable::Align::CENTER});
    // print out the titles
    table << titles << energies_Ry << energies_eV;

    GlobalV::ofs_running << table.str() << std::endl;


    
    if (PARAM.inp.out_level == "ie" || PARAM.inp.out_level == "m")
    {
        std::vector<double> mag;
        switch (PARAM.inp.nspin)
        {
        case 2:
            mag = {magnet.tot_mag, magnet.abs_mag};
            break;
        case 4:
            mag = {magnet.tot_mag_nc[0],
                   magnet.tot_mag_nc[1],
                   magnet.tot_mag_nc[2],
                   magnet.abs_mag};
            break;
        default:
            mag = {};
            break;
        }
        std::vector<double> drho = {scf_thr};
        if (XC_Functional::get_ked_flag())
        {
            drho.push_back(scf_thr_kin);
        }
        elecstate::print_scf_iterinfo(PARAM.inp.ks_solver,
                                      iter,
                                      6,
                                      mag,
                                      10,
                                      elec.f_en.etot * ModuleBase::Ry_to_eV,
                                      elec.f_en.etot_delta * ModuleBase::Ry_to_eV,
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
void print_format(const std::string& name, const double& value)
{
    GlobalV::ofs_running << std::setiosflags(std::ios::showpos);
    GlobalV::ofs_running << " " << std::setw(16) << name << std::setw(30) << value << std::setw(30)
                         << value * ModuleBase::Ry_to_eV << std::endl;
    GlobalV::ofs_running << std::resetiosflags(std::ios::showpos);
    return;
}
} // namespace elecstate

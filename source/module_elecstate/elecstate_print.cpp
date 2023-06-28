#include "elecstate.h"
#include "occupy.h"
#include "elecstate_getters.h"
#include "module_base/global_variable.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_elecstate/potentials/efield.h"
#include "module_elecstate/potentials/gatefield.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include "module_hamilt_lcao/module_deepks/LCAO_deepks.h"
namespace elecstate
{
/// @brief print and check for band energy and occupations
/// @param ofs
void ElecState::print_eigenvalue(std::ofstream& ofs)
{
    bool wrong = false;
    for (int ik = 0; ik < this->klist->nks; ++ik)
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
    for (int ik = 0; ik < this->klist->nks; ik++)
    {
        ofs << std::setprecision(5);
        ofs << std::setiosflags(ios::showpoint);
        if (ik == 0)
        {
            ofs << "   NSPIN == " << GlobalV::NSPIN << std::endl;
            if (GlobalV::NSPIN == 2)
            {
                ofs << "SPIN UP : " << std::endl;
            }
        }
        else if (ik == this->klist->nks / 2)
        {
            if (GlobalV::NSPIN == 2)
            {
                ofs << "SPIN DOWN : " << std::endl;
            }
        }

        if (GlobalV::NSPIN == 2)
        {
            if (this->klist->isk[ik] == 0)
            {
                ofs << " " << ik + 1 << "/" << this->klist->nks / 2 
                << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x << " " << this->klist->kvec_c[ik].y << " " << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);
            }
            if (this->klist->isk[ik] == 1)
            {
                ofs << " " << ik + 1 - this->klist->nks / 2 << "/" << this->klist->nks / 2
                    << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x << " " << this->klist->kvec_c[ik].y << " "
                    << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik] << " pws)" << std::endl;

                ofs << std::setprecision(6);
            }
        } // Pengfei Li  added  14-9-9
        else
        {
            ofs << " " << ik + 1 << "/" << this->klist->nks << " kpoint (Cartesian) = " << this->klist->kvec_c[ik].x
                << " " << this->klist->kvec_c[ik].y << " " << this->klist->kvec_c[ik].z << " (" << this->klist->ngk[ik]
                << " pws)" << std::endl;

            ofs << std::setprecision(6);
        }

        ofs << std::setprecision(6);
        ofs << std::setiosflags(ios::showpoint);
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
    for (int ib = 0; ib < GlobalV::NBANDS; ++ib)
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
            GlobalV::ofs_running << " Energy (eV) & Occupations  for spin=" << GlobalV::CURRENT_SPIN + 1
                                 << " K-point=" << ik + 1 << std::endl;
            GlobalV::ofs_running << std::setiosflags(ios::showpoint);
            for (int ib = 0; ib < GlobalV::NBANDS; ib++)
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
/// @param converged: if converged
/// @param iter_in: iter
/// @param scf_thr: threshold for scf
/// @param duration: time of each iteration
/// @param pw_diag_thr: threshold for diagonalization
/// @param avg_iter: averaged diagonalization iteration of each scf iteration
/// @param print: if print to screen
void ElecState::print_etot(const bool converged,
                           const int& iter_in,
                           const double& scf_thr,
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
    GlobalV::ofs_running << std::setiosflags(ios::right);

    GlobalV::ofs_running << "\n Density error is " << scf_thr << std::endl;

    if (GlobalV::BASIS_TYPE == "pw")
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "Error Threshold", pw_diag_thr); // xiaohui add 2013-09-02

    if (printe > 0 && ((iter + 1) % printe == 0 || converged || iter == GlobalV::SCF_NMAX))
    {
        int n_order = std::max(0, Occupy::gaussian_type);
        GlobalV::ofs_running << "\n " << std::setw(16) << "Energy" << std::setw(30) << "Rydberg" << std::setw(30)
                             << "eV" << std::endl;
        this->print_format("E_KohnSham", this->f_en.etot);
        this->print_format("E_KS(sigma->0)", this->f_en.etot - this->f_en.demet/(2+n_order));
        this->print_format("E_Harris", this->f_en.etot_harris);
        this->print_format("E_band", this->f_en.eband);
        this->print_format("E_one_elec", this->f_en.eband + this->f_en.deband);
        this->print_format("E_Hartree", this->f_en.hartree_energy);
        this->print_format("E_xc", this->f_en.etxc - this->f_en.etxcc);
        this->print_format("E_Ewald", this->f_en.ewald_energy);
        this->print_format("E_entropy(-TS)", this->f_en.demet); // mohan add 2011-12-02
        this->print_format("E_descf", this->f_en.descf);
        std::string vdw_method = get_input_vdw_method();
        if (vdw_method == "d2") // Peize Lin add 2014-04, update 2021-03-09
        {
            this->print_format("E_vdwD2", this->f_en.evdw);
        }
        else if (vdw_method == "d3_0" || vdw_method == "d3_bj") // jiyy add 2019-05, update 2021-05-02
        {
            this->print_format("E_vdwD3", this->f_en.evdw);
        }
        this->print_format("E_exx", this->f_en.exx);
        if (GlobalV::imp_sol)
        {
            this->print_format("E_sol_el", this->f_en.esol_el);
            this->print_format("E_sol_cav", this->f_en.esol_cav);
        }
        if (GlobalV::EFIELD_FLAG)
        {
            this->print_format("E_efield", elecstate::Efield::etotefield);
        }
        if (GlobalV::GATE_FLAG)
        {
            this->print_format("E_gatefield", elecstate::Gatefield::etotgatefield);
        }

#ifdef __DEEPKS
        if (GlobalV::deepks_scf) // caoyu add 2021-08-10
        {
            this->print_format("E_DeePKS", GlobalC::ld.E_delta);
        }
#endif
    }
    else
    {
        GlobalV::ofs_running << "\n " << std::setw(12) << "Energy" << std::setw(30) << "Rydberg" << std::setw(30)
                             << "eV" << std::endl;
        this->print_format("E_KohnSham", this->f_en.etot);
        this->print_format("E_Harris", this->f_en.etot_harris);
    }

    if (GlobalV::TWO_EFERMI)
    {
        this->print_format("E_Fermi_up", this->eferm.ef_up);
        this->print_format("E_Fermi_dw", this->eferm.ef_dw);
    }
    else
    {
        this->print_format("E_Fermi", this->eferm.ef);
    }
    if (GlobalV::out_bandgap)
    {
        if (!GlobalV::TWO_EFERMI)
        {
            this->print_format("E_bandgap", this->bandgap);
        }
        else
        {
            this->print_format("E_bandgap_up", this->bandgap_up);
            this->print_format("E_bandgap_dw", this->bandgap_dw);
        }
    }

    if (iter_in == 1) // pengfei Li added 2015-1-31
    {
        this->f_en.etot_old = this->f_en.etot;
    }

    // mohan update 2011-02-26
    std::stringstream ss;

    // xiaohui add 2013-09-02, Peize Lin update 2020.11.14
    std::string label;
    std::string ks_solver_type = get_ks_solver_type();
    if (ks_solver_type == "cg")
    {
        label = "CG";
    }
    else if (ks_solver_type == "lapack")
    {
        label = "LA";
    }
    else if (ks_solver_type == "genelpa")
    {
        label = "GE";
    }
    else if (ks_solver_type == "dav")
    {
        label = "DA";
    }
    else if (ks_solver_type == "scalapack_gvx")
    {
        label = "GV";
    }
    else if (ks_solver_type == "cusolver")
    {
        label = "CU";
    }
    else
    {
        ModuleBase::WARNING_QUIT("Energy", "print_etot found unknown ks_solver_type");
    }
    ss << label << iter;
    // xiaohui add 2013-09-02

    bool scientific = true;
    int prec = 6;

    if (!print)
        return;

    if (GlobalV::OUT_LEVEL == "ie" || GlobalV::OUT_LEVEL == "m") // xiaohui add 'm' option, 2015-09-16
    {
        std::cout << " " << std::setw(7) << ss.str();
        // std::cout << std::setiosflags(ios::fixed);
        // std::cout << std::setiosflags(ios::showpos);
        if (scientific)
        {
            std::cout << std::setiosflags(ios::scientific);
        }

        if (GlobalV::COLOUR)
        {
            if (GlobalV::MY_RANK == 0)
            {
                printf("\e[36m%-15f\e[0m", this->f_en.etot);
                if (GlobalV::NSPIN == 2)
                {
                    std::cout << std::setprecision(2);
                    std::cout << std::setw(10) << get_ucell_tot_magnetization();
                    std::cout << std::setw(10) << get_ucell_abs_magnetization();
                }
                else if (GlobalV::NSPIN == 4 && GlobalV::NONCOLIN)
                {
                    std::cout << std::setprecision(2);
                    std::cout << std::setw(10) << get_ucell_tot_magnetization_nc_x() << std::setw(10)
                              << get_ucell_tot_magnetization_nc_y() << std::setw(10)
                              << get_ucell_tot_magnetization_nc_z();
                    std::cout << std::setw(10) << get_ucell_abs_magnetization();
                }
                if (scf_thr > 1.0)
                {
                    // 31 is red
                    printf("\e[31m%-14e\e[0m", scf_thr);
                    // printf( "[31m%-14e[0m", scf_thr);
                }
                else
                {
                    // 32 is green
                    printf("\e[32m%-14e\e[0m", scf_thr);
                    // printf( "[32m%-14e[0m", scf_thr);
                }
                // 34 is blue
                printf("\e[36m%-15f\e[0m", this->f_en.etot * ModuleBase::Ry_to_eV);
                std::cout << std::setprecision(3);
                std::cout << std::resetiosflags(ios::scientific);

                std::cout << std::setw(11) << duration;
                std::cout << std::endl;
            }
        }
        else
        {
            std::cout << std::setprecision(prec);
            if (GlobalV::NSPIN == 2)
            {
                std::cout << std::setprecision(2);
                std::cout << std::setw(10) << get_ucell_tot_magnetization();
                std::cout << std::setw(10) << get_ucell_abs_magnetization();
            }
            std::cout << std::setprecision(6);
            std::cout << std::setw(15) << this->f_en.etot * ModuleBase::Ry_to_eV;
            std::cout << std::setw(15) << (this->f_en.etot - this->f_en.etot_old) * ModuleBase::Ry_to_eV;
            std::cout << std::setprecision(3);
            std::cout << std::setw(11) << scf_thr;
            std::cout << std::setprecision(3);
            std::cout << std::setw(11) << duration;
            std::cout << std::endl;
        }
    }
    else
    {
    }

    this->f_en.etot_old = this->f_en.etot;
    return;
}

/// @brief function to print name, value and value*Ry_to_eV
/// @param name: name
/// @param value: value
void ElecState::print_format(const std::string& name, const double& value)
{
    GlobalV::ofs_running << std::setiosflags(ios::showpos);
    std::stringstream name2;
    name2 << name;
    GlobalV::ofs_running << " " << std::setw(16) << name2.str() << std::setw(30) << value << std::setw(30)
                         << value * ModuleBase::Ry_to_eV << std::endl;
    GlobalV::ofs_running << std::resetiosflags(ios::showpos);
    return;
}

} // namespace elecstate
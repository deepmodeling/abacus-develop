#include "read_pseudo.h"
#include "source_base/global_file.h"
#include "cal_atoms_info.h"
#include "read_pp.h"
#include "bcast_cell.h"
#include "source_base/element_elec_config.h"
#include "source_base/parallel_common.h"

#include <cstring>

namespace unitcell {
AtomsInfoResult read_pseudo(std::ofstream& ofs, UnitCell& ucell,
                   const std::string& pseudo_dir,
                   const std::string& global_out_dir,
                   const bool out_element_info,
                   const std::string& dft_functional,
                   const bool lspinorb,
                   const double pseudo_rcut,
                   const double soc_lambda,
                   const int nspin,
                   const int npol,
                   const std::string& basis_type,
                   const std::string& esolver_type,
                   const std::string& init_wfc,
                   const int nbands,
                   const bool two_fermi,
                   const double nelec_delta,
                   const std::string& smearing_method,
                   const std::string& ks_solver,
                   const int bndpar,
                   const double nelec,
                   const double nupdown) {
    // read in non-local pseudopotential and ouput the projectors.
    ofs << "\n\n";
    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " |                 #Read Pseudopotentials Files#                      |" << std::endl;
    ofs << " | ABACUS supports norm-conserving (NC) pseudopotentials for both     |" << std::endl;
    ofs << " | plane wave basis set and numerical atomic orbital basis set.       |" << std::endl;
    ofs << " | In addition, ABACUS supports ultrasoft pseudopotentials (USPP)     |" << std::endl;
    ofs << " | for plane wave basis set.                                          |" << std::endl;
    ofs << " |                                                                    |" << std::endl;
    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
    ofs << "\n";

    const std::string pseudo_dir_ = pseudo_dir;
    const std::string global_out_dir_ = global_out_dir;
    const bool out_element_info_ = out_element_info;
    const std::string dft_functional_ = dft_functional;
    read_cell_pseudopots(pseudo_dir_, ofs, ucell, global_out_dir_, dft_functional_, lspinorb, pseudo_rcut, soc_lambda);

	if (GlobalV::MY_RANK == 0) 
	{
		for (int it = 0; it < ucell.ntype; it++) 
		{
			Atom* atom = &ucell.atoms[it];
			if (!(atom->label_orb.empty())) 
			{
                ucell.compare_atom_labels(atom->label_orb, atom->ncpp.psd);
            }
        }

		if (out_element_info_) 
		{
			for (int i = 0; i < ucell.ntype; i++) 
			{
				ModuleBase::Global_File::make_dir_atom(ucell.atoms[i].label, global_out_dir_);
            }
			for (int it = 0; it < ucell.ntype; it++) 
			{
				Atom* atom = &ucell.atoms[it];
                std::stringstream ss;
                ss << global_out_dir_ << atom->label << "/"
                   << atom->label << ".NONLOCAL";
                std::ofstream ofs(ss.str().c_str());

                ofs << "<HEADER>" << std::endl;
                ofs << std::setw(10) << atom->label << "\t"
                    << "label" << std::endl;
                ofs << std::setw(10) << atom->ncpp.pp_type << "\t"
                    << "Pseudopotential type" << std::endl;
                ofs << std::setw(10) << atom->ncpp.lmax << "\t"
                    << "lmax" << std::endl;
                ofs << "</HEADER>" << std::endl;

                ofs << "\n<DIJ>" << std::endl;
                ofs << std::setw(10) << atom->ncpp.nbeta << "\t"
                    << "nummber of projectors." << std::endl;
                for (int ib = 0; ib < atom->ncpp.nbeta; ib++) {
                    for (int ib2 = 0; ib2 < atom->ncpp.nbeta; ib2++) {
                        ofs << std::setw(10) << atom->ncpp.lll[ib] << " "
                            << atom->ncpp.lll[ib2] << " "
                            << atom->ncpp.dion(ib, ib2) << std::endl;
                    }
                }
                ofs << "</DIJ>" << std::endl;

                for (int i = 0; i < atom->ncpp.nbeta; i++) {
                    ofs << "<PP_BETA>" << std::endl;
                    ofs << std::setw(10) << i << "\t"
                        << "the index of projectors." << std::endl;
                    ofs << std::setw(10) << atom->ncpp.lll[i] << "\t"
                        << "the angular momentum." << std::endl;

                    // mohan add
                    // only keep the nonzero part.
                    int cut_mesh = atom->ncpp.mesh;
                    for (int j = atom->ncpp.mesh - 1; j >= 0; --j) {
                        if (std::abs(atom->ncpp.betar(i, j)) > 1.0e-10) {
                            cut_mesh = j;
                            break;
                        }
                    }
                    if (cut_mesh % 2 == 0) {
                        ++cut_mesh;
                    }

                    ofs << std::setw(10) << cut_mesh << "\t"
                        << "the number of mesh points." << std::endl;

                    for (int j = 0; j < cut_mesh; ++j) {
                        ofs << std::setw(15) << atom->ncpp.r[j] << std::setw(15)
                            << atom->ncpp.betar(i, j) << std::setw(15)
                            << atom->ncpp.rab[j] << std::endl;
                    }
                    ofs << "</PP_BETA>" << std::endl;
                }

                ofs.close();
            }
        }
    }

#ifdef __MPI
    bcast_atoms_pseudo(ucell.atoms,ucell.ntype);
#endif

    for (int it = 0; it < ucell.ntype; it++) {
        if (ucell.atoms[0].ncpp.xc_func != ucell.atoms[it].ncpp.xc_func) {
            GlobalV::ofs_warning << "\n type " << ucell.atoms[0].label
                                 << " functional is " << ucell.atoms[0].ncpp.xc_func;

            GlobalV::ofs_warning << "\n type " << ucell.atoms[it].label
                                 << " functional is " << ucell.atoms[it].ncpp.xc_func
                                 << std::endl;

            ModuleBase::WARNING_QUIT("setup_cell",
                                     "All DFT functional must consistent.");
        }
    }

    // setup the total number of PAOs
    cal_natomwfc(ofs,ucell.natomwfc,ucell.ntype,ucell.atoms,nspin);

    // Calculate the information of atoms from the pseudopotential
    // CRITICAL: Must pass the user-specified nbands and nelec parameters to cal_atoms_info().
    // Previously, nbands and nelec were not passed, causing cal_atoms_info() to use default 0,
    // which triggered cal_nbands() and cal_nelec() to auto-calculate regardless of user input.
    // This led to incorrect energy calculations (deviation ~139 eV in test 006_PW_UPF201_Eu,
    // and ~7-8 eV in tests 076_PW_elec_add, 078_PW_S2_elec_add, 082_PW_gatefield).
    CalAtomsInfo ca;
    AtomsInfoResult atoms_info = ca.cal_atoms_info(ucell.atoms, ucell.ntype,
                                                    nspin, two_fermi, nelec_delta,
                                                    esolver_type, lspinorb,
                                                    basis_type, smearing_method,
                                                    ks_solver, bndpar,
                                                    nbands,
                                                    nelec,
                                                    nupdown);

    // setup nlocal
    // nlocal is calculated by CalAtomsInfo::cal_atoms_info() above
    // Use the input nbands parameter (from user specification) instead of atoms_info.nbands
    cal_nwfc(ofs, ucell, ucell.atoms, nspin, atoms_info.nlocal, npol,
              basis_type, esolver_type, init_wfc, nbands);

    // Check whether the number of valence is minimum
	if (GlobalV::MY_RANK == 0) 
	{
		int abtype = 0;
		for (int it = 0; it < ucell.ntype; it++) 
		{
			if (ModuleBase::MinZval.find(ucell.atoms[it].ncpp.psd)
					!= ModuleBase::MinZval.end()) 
			{
				if (ucell.atoms[it].ncpp.zv
						> ModuleBase::MinZval.at(ucell.atoms[it].ncpp.psd)) 
				{
                    abtype += 1;
					if (abtype == 1) 
					{
                        std::cout << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                     "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                     "%%%%%%%%%%%%%%%%%%%%%%%%%%"
                                  << std::endl;
                        ofs << "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                               "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                               "%%%%%%%%%%%%%%%%%%%%%"
                            << std::endl;
                    }
                    std::cout << " Warning: the number of valence electrons in "
                                 "pseudopotential > "
                              << ModuleBase::MinZval.at(ucell.atoms[it].ncpp.psd);
                    std::cout << " for " << ucell.atoms[it].ncpp.psd << ": "
                              << ModuleBase::EleConfig.at(ucell.atoms[it].ncpp.psd)
                              << std::endl;
                    ofs << " Warning: the number of valence electrons in "
                           "pseudopotential > "
                        << ModuleBase::MinZval.at(ucell.atoms[it].ncpp.psd);
                    ofs << " for " << ucell.atoms[it].ncpp.psd << ": "
                        << ModuleBase::EleConfig.at(ucell.atoms[it].ncpp.psd)
                        << std::endl;
                }
            }
        }
		if (abtype > 0) 
		{
			std::cout << " Pseudopotentials with additional electrons can "
                         "yield (more) accurate outcomes, but may be "
                         "less efficient."
                      << std::endl;
            std::cout
                << " If you're confident that your chosen pseudopotential is "
                   "appropriate, you can safely ignore "
                   "this warning."
                << std::endl;
            std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                         "%%%%%%%%%%%%\n"
                      << std::endl;
            ofs << " Pseudopotentials with additional electrons can yield "
                   "(more) accurate outcomes, but may be less "
                   "efficient."
                << std::endl;
            ofs << " If you're confident that your chosen pseudopotential is "
                   "appropriate, you can safely ignore this "
                   "warning."
                << std::endl;
            ofs << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                   "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
                   "%%%%%%%";
            ModuleBase::GlobalFunc::OUT(ofs, "");
        }
    }

    cal_meshx(ucell.meshx,ucell.atoms,ucell.ntype);

#ifdef __MPI
    Parallel_Common::bcast_int(ucell.meshx);
    Parallel_Common::bcast_int(ucell.natomwfc);
    Parallel_Common::bcast_int(ucell.lmax);
    Parallel_Common::bcast_int(ucell.lmax_ppwf);
#endif

    return atoms_info;
}

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void read_cell_pseudopots(const std::string& pp_dir, std::ofstream& log, UnitCell& ucell,
                          const std::string& global_out_dir,
                          const std::string& dft_functional,
                          const bool lspinorb,
                          const double pseudo_rcut,
                          const double soc_lambda)
{
    ModuleBase::TITLE("UnitCell", "read_cell_pseudopots");
    // setup reading log for pseudopot_upf
    const std::string global_out_dir_ = global_out_dir;
    const std::string dft_functional_ = dft_functional;
    const bool lspinorb_ = lspinorb;
    const double pseudo_rcut_ = pseudo_rcut;
    const double soc_lambda_ = soc_lambda;
    std::stringstream ss;
    ss << global_out_dir_ << "atom_pseudo.log";

    // Read in the atomic pseudo potentials
    std::string pp_address;
    for (int i = 0; i < ucell.ntype; i++)
    {
        Pseudopot_upf upf;
        upf.coulomb_potential = ucell.atoms[i].coulomb_potential;

        // mohan update 2010-09-12
        int error = 0;
        int error_ap = 0;

        if (GlobalV::MY_RANK == 0)
        {
            pp_address = pp_dir + ucell.pseudo_fn[i];
            error = upf.init_pseudo_reader(pp_address, ucell.pseudo_type[i], ucell.atoms[i].ncpp); // xiaohui add 2013-06-23

            if (error == 0) // mohan add 2021-04-16
            {
                if (ucell.atoms[i].flag_empty_element) // Peize Lin add for bsse 2021.04.07
                {
                    upf.set_empty_element(ucell.atoms[i].ncpp);
                }
                upf.set_upf_q(ucell.atoms[i].ncpp); // liuyu add 2023-09-21
                // average pseudopotential if needed
                error_ap = upf.average_p(soc_lambda_, ucell.atoms[i].ncpp, lspinorb_);
            }
            ucell.atoms[i].coulomb_potential = upf.coulomb_potential;
        }

#ifdef __MPI
        Parallel_Common::bcast_int(error);
        Parallel_Common::bcast_int(error_ap);
        Parallel_Common::bcast_bool(ucell.atoms[i].coulomb_potential);
#endif

        if (error_ap)
        {
            ModuleBase::WARNING_QUIT("read_cell_pseudopots", "error when average the pseudopotential.");
        }

        if (error == 1)
        {
            std::cout << " Pseudopotential directory now is : " << pp_address << std::endl;
            GlobalV::ofs_warning << " Pseudopotential directory now is : " << pp_address << std::endl;
            ModuleBase::WARNING_QUIT("read_cell_pseudopots", "Couldn't find pseudopotential file.");
        }
        else if (error == 2)
        {
            ModuleBase::WARNING_QUIT("read_cell_pseudopots", "Pseudopotential data do not match.");
        }
        else if (error == 3)
        {
            ModuleBase::WARNING_QUIT(
                "read_cell_pseudopots",
                "Check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave "
                "functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n "
                "unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave "
                "functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
        }
        else if (error == 4)
        {
            ModuleBase::WARNING_QUIT("read_cell_pseudopots", "Unknown pseudopotential type.");
        }

        if (GlobalV::MY_RANK == 0)
        {
		    upf.complete_default(ucell.atoms[i].ncpp, pseudo_rcut_);

            log << std::endl;
            ModuleBase::GlobalFunc::OUT(log, "Pseudopotential file", ucell.pseudo_fn[i]);
            ModuleBase::GlobalFunc::OUT(log, "Pseudopotential type", ucell.atoms[i].ncpp.pp_type);
            ModuleBase::GlobalFunc::OUT(log, "Exchange-correlation functional", ucell.atoms[i].ncpp.xc_func);
            ModuleBase::GlobalFunc::OUT(log, "Nonlocal core correction", ucell.atoms[i].ncpp.nlcc);
            // ModuleBase::GlobalFunc::OUT(log, "spin orbital", ucell.atoms[i].has_so);
            ModuleBase::GlobalFunc::OUT(log, "Valence electrons", ucell.atoms[i].ncpp.zv);
            ModuleBase::GlobalFunc::OUT(log, "Lmax", ucell.atoms[i].ncpp.lmax);
            ModuleBase::GlobalFunc::OUT(log, "Number of zeta", ucell.atoms[i].ncpp.nchi);
            ModuleBase::GlobalFunc::OUT(log, "Number of projectors", ucell.atoms[i].ncpp.nbeta);
            for (int ib = 0; ib < ucell.atoms[i].ncpp.nbeta; ib++)
            {
                ModuleBase::GlobalFunc::OUT(log, "L of projector", ucell.atoms[i].ncpp.lll[ib]);
            }
            //			ModuleBase::GlobalFunc::OUT(log,"Grid Mesh Number", atoms[i].mesh);
            if (dft_functional_ != "default")
            {
                std::string xc_func1 = dft_functional_;
                transform(xc_func1.begin(), xc_func1.end(), xc_func1.begin(), (::toupper));
                if (xc_func1 != ucell.atoms[i].ncpp.xc_func)
                {
                    std::cout << " NAME OF ELEMENT      : " << ucell.atoms[i].label << std::endl;
                    std::cout << " DFT FUNC. (PSEUDO)   : " << ucell.atoms[i].ncpp.xc_func << std::endl;
                    std::cout << " DFT FUNC. (SET TO)   : " << xc_func1 << std::endl; 
                    std::cout << " MAKE SURE THIS DFT FUNCTIONAL IS WHAT YOU NEED" << std::endl;


                    GlobalV::ofs_warning << " NAME OF ELEMENT      : " << ucell.atoms[i].label << std::endl;
                    GlobalV::ofs_warning << " DFT FUNC. (PSEUDO)   : " << ucell.atoms[i].ncpp.xc_func << std::endl;
                    GlobalV::ofs_warning << " DFT FUNC. (SET TO)   : " << xc_func1 << std::endl; 
                    GlobalV::ofs_warning << " MAKE SURE THIS DFT FUNCTIONAL IS WHAT YOU NEED" << std::endl;

                    ucell.atoms[i].ncpp.xc_func = xc_func1;
                    ModuleBase::GlobalFunc::OUT(log, "DFT functional set to", xc_func1);
                }
            }
        }
    }
    return;
}

void print_unitcell_pseudo(const std::string& fn, UnitCell& ucell)
{
    ModuleBase::TITLE("unitcell", "print_unitcell_pseudo");
    std::ofstream ofs(fn.c_str());

    ucell.print_cell(ofs);
    for (int i = 0; i < ucell.ntype; i++)
    {
        ucell.atoms[i].print_Atom(ofs);
    }

    ofs.close();
    return;
}

}

#include "read_input.h"

#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"

#include <fstream>
#include <iostream>
#include <string.h>
namespace ModuleIO
{

void strtolower(char* sa, char* sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

void read_information(std::ifstream& ifs, std::vector<std::string>& output, const std::string& delimiters)
{
    std::string value1;
    ifs >> value1;
    output.push_back(value1);

    std::string line;
    getline(ifs, line);

    std::istringstream iss(line);
    std::string word;
    while (iss >> word)
    {
        if (delimiters.find(word[0]) != std::string::npos)
            break;
        output.push_back(word);
    }
}

bool ReadInput::check_mode = false;

ReadInput::ReadInput(const int& rank)
{
    this->rank = rank;
    this->item_general();
    this->item_pw();
    this->item_sdft();
    this->item_relax();
    this->item_lcao();
    this->item_postprocess();
    this->item_md();
    this->item_others();
}

void ReadInput::readin_parameters(Parameter& param, const std::string& filename_in, const std::string& filename_out)
{
    // 1. only rank 0 read the input file
    if (this->rank == 0)
    {
        this->read_txt_input(param, filename_in);
    }
    if (this->check_mode)
    {
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "  INPUT parameters have been successfully checked!" << std::endl;
        std::cout << "----------------------------------------------------------" << std::endl;
        exit(0);
        return;
    }
    if (filename_out != "")
    {
        this->write_txt_input(param, filename_out);
    }
#ifdef __MPI
    // 2. broadcast the parameters
    for (auto& bcastfunc: this->bcastfuncs)
    {
        bcastfunc(param);
    }
#endif
}

void ReadInput::read_txt_input(Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("ReadInput", "read_txt_input");

    std::ifstream ifs(filename.c_str(), std::ios::in);

    if (!ifs)
    {
        std::cout << " Can't find the INPUT file." << std::endl;
        ModuleBase::WARNING_QUIT("Input::Init", "Error during readin parameters.", 1);
    }

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    // ifs >> std::setiosflags(ios::uppercase);
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n');
        if (strcmp(word, "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
        ifs.rdstate();
    }

    if (ierr == 0)
    {
        std::cout << " Error parameter list. "
                  << " The parameter list always starts with key word 'INPUT_PARAMETERS'. " << std::endl;
        ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
    }

    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word1;
        if (ifs.eof())
            break;
        strtolower(word1, word);

        auto it = this->input_lists.find(word);
        if (it != this->input_lists.end())
        {
            Input_Item* p_item = &(it->second);
            this->readvalue_items.push_back(p_item);
            if (it->second.checkvalue)
                this->checkvalue_items.push_back(p_item);
            if (it->second.resetvalue)
                this->resetvalue_items.push_back(p_item);
            read_information(ifs, it->second.str_values, "#/!");
        }
        else
        {
            if (word[0] != '#' && word[0] != '/' && word[0] != '!')
            {
                std::cout << " THE PARAMETER NAME '" << word << "' IS NOT USED!" << std::endl;
                ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
            }
            ifs.ignore(150, '\n');
        }
    }

    // 1) read the value of the parameters
    for (auto& readvalue_item: this->readvalue_items)
    {
        readvalue_item->readvalue(*readvalue_item, param);
    }

    // 2) reset the value of some parameters based on readin values
    //    e.g. if (calulation_type == "nscf") then set "init_chg" to "file".
    for (auto& resetvalue_item: this->resetvalue_items)
    {
        resetvalue_item->resetvalue(*resetvalue_item, param);
    }

    // 3) check the value of the parameters
    for (auto& checkvalue_item: this->checkvalue_items)
    {
        checkvalue_item->checkvalue(*checkvalue_item, param);
    }
}

void ReadInput::write_txt_input(const Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("ReadInput", "write_txt_input");
    for (auto& item: this->input_lists)
    {
        item.second.getfinalvalue(item.second, param);
    }
    std::ofstream ofs(filename.c_str(), std::ios::out);
    ofs << "INPUT_PARAMETERS" << std::endl;
    ofs << std::setiosflags(std::ios::left);

    ofs << "#Parameters (1.General)" << std::endl;
    for (auto& item: this->input_lists)
    {
        const Input_Item* p_item = &(item.second);
        if(p_item->label == "ecutwfc")
        {
            ofs << "\n#Parameters (2.PW)" << std::endl;
        }
        else if(p_item->label == "method_sto")
        {
            ofs << "\n#Parameters (3.Stochastic DFT)" << std::endl;
        }
        else if(p_item->label == "ks_solver")
        {
            ofs << "\n#Parameters (4.Relaxation)" << std::endl;
        }
        else if(p_item->label == "basis_type")
        {
            ofs << "\n#Parameters (5.LCAO)" << std::endl;
        }
        else if(p_item->label == "smearing_method")
        {
            ofs << "\n#Parameters (6.Smearing)" << std::endl;
        }
        else if(p_item->label == "mixing_type")
        {
            ofs << "\n#Parameters (7.Charge Mixing)" << std::endl;
        }
        else if(p_item->label == "dos_emin_ev")
        {
            ofs << "\n#Parameters (8.DOS)" << std::endl;
        }
        else if(p_item->label == "md_type")
        {
            ofs << "\n#Parameters (9.Molecular dynamics)" << std::endl;
        }
        else if(p_item->label == "efield_flag")
        {
            ofs << "\n#Parameters (10.Electric field and dipole correction)" << std::endl;
        }
        else if(p_item->label == "gate_flag")
        {
            ofs << "\n#Parameters (11.Gate field)" << std::endl;
        }
        else if(p_item->label == "out_alllog")
        {
            ofs << "\n#Parameters (12.Test)" << std::endl;
        }
        else if(p_item->label == "vdw_method")
        {
            ofs << "\n#Parameters (13.vdW Correction)" << std::endl;
        }
        else if(p_item->label == "exx_hybrid_alpha")
        {
            ofs << "\n#Parameters (14.exx)" << std::endl;
        }
        else if(p_item->label == "td_force_dt")
        {
            ofs << "\n#Parameters (16.tddft)" << std::endl;
        }
        else if(p_item->label == "berry_phase")
        {
            ofs << "\n#Parameters (17.berry_wannier)" << std::endl;
        }
        else if(p_item->label == "imp_sol")
        {
            ofs << "\n#Parameters (18.implicit_solvation)" << std::endl;
        }
        else if(p_item->label == "of_kinetic")
        {
            ofs << "\n#Parameters (19.orbital free density functional theory)" << std::endl;
        }
        else if(p_item->label == "dft_plus_u")
        {
            ofs << "\n#Parameters (20.dft+u)" << std::endl;
        }
        else if(p_item->label == "bessel_nao_ecut")
        {
            ofs << "\n#Parameters (21.spherical bessel)" << std::endl;
        }
        else if(p_item->label == "sc_mag_switch")
        {
            ofs << "\n#Parameters (22.non-collinear spin-constrained DFT)" << std::endl;
        }
        else if(p_item->label == "qo_switch")
        {
            ofs << "\n#Parameters (23.Quasiatomic Orbital analysis)" << std::endl;
        }
        else if(p_item->label == "pexsi_npole")
        {
            ofs << "\n#Parameters (24.PEXSI)" << std::endl;
        }
        ModuleBase::GlobalFunc::OUTP(ofs, p_item->label, p_item->final_value.str(), p_item->annotation);
    }
}

void ReadInput::add_item(const Input_Item& item)
{
    // only rank 0 read the input file
    // So only rank 0 add the item to the input list
    if (this->rank == 0)
    {
        this->input_lists.insert(make_pair(item.label, item));
    }
}


} // namespace ModuleIO
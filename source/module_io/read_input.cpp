#include "read_input.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>

#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"
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

std::string longstring(const std::vector<std::string>& str_values, const int length)
{
    std::string output;
    output = "";
    for (int i = 0; i < length; ++i)
    {
        output += str_values[i];
        if (i != length - 1)
        {
            output += " ";
        }
    }
    return output;
}

bool convert_bool(std::string str)
{
    for (auto& i: str)
    {
        i = tolower(i);
    }
    if (str == "true")
    {
        return true;
    }
    else if (str == "false")
    {
        return false;
    }
    else if (str == "1")
    {
        return true;
    }
    else if (str == "0")
    {
        return false;
    }
    else if (str == "t")
    {
        return true;
    }
    else if (str == "f")
    {
        return false;
    }
    else
    {
        std::string warningstr = "Bad boolean parameter ";
        warningstr.append(str);
        warningstr.append(", please check the input parameters in file INPUT");
        ModuleBase::WARNING_QUIT("Input", warningstr);
    }
}

void read_information(std::ifstream& ifs, std::vector<std::string>& output, const std::string& delimiters)
{
    std::string line;
    getline(ifs, line);

    std::istringstream iss(line);
    std::string word;
    while (iss >> word)
    {
        if (delimiters.find(word[0]) != std::string::npos) {
            break;
}
        output.push_back(word);
    }
    if (output.size() == 0)
    {
        output.push_back("");
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

void ReadInput::read_parameters(Parameter& param, const std::string& filename_in)
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

    // 2. broadcast the parameters
    for (auto& bcastfunc: this->bcastfuncs)
    {
        bcastfunc(param);
    }
}

void ReadInput::create_directory(const Parameter& param)
{
    ModuleBase::TITLE("ReadInput", "create_directory");

    // mohan move forward 2011-02-26
    //----------------------------------------------------------
    // OTHRE CLASS MEMBER FUNCTION :
    // NAME : Run::make_dir( dir name : OUT.suffix)
    //----------------------------------------------------------
    bool out_dir = false;
    if (!param.input.out_app_flag
        && (param.input.out_mat_hs2 || param.input.out_mat_r || param.input.out_mat_t || param.input.out_mat_dh))
    {
        out_dir = true;
    }
    // NOTE: "make_dir_out" must be called by all processes!!!
    //       Maybe it is not good, because only rank 0 can create the directory.
    ModuleBase::Global_File::make_dir_out(param.input.suffix,
                                          param.input.calculation,
                                          out_dir,
                                          this->rank,
                                          param.input.mdp.md_restart,
                                          param.input.out_alllog); // xiaohui add 2013-09-01

    const std::string ss = "test -d " + PARAM.inp.read_file_dir;
    if (system(ss.c_str()))
    {
        ModuleBase::WARNING_QUIT("ReadInput", "please set right files directory for reading in.");
    }

    return;
}

void ReadInput::write_parameters(const Parameter& param, const std::string& filename_out)
{
    if (this->rank == 0)
    {
        this->write_txt_input(param, filename_out);
    }
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
                  << " The parameter list always starts with key word "
                     "'INPUT_PARAMETERS'. "
                  << std::endl;
        ModuleBase::WARNING_QUIT("Input", "Bad parameter, please check the input parameters in file INPUT", 1);
    }

    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word1;
        if (ifs.eof()) {
            break;
}
        strtolower(word1, word);
        auto it = std::find_if(input_lists.begin(),
                               input_lists.end(),
                               [&word](const std::pair<std::string, Input_Item>& item) { return item.first == word; });
        if (it != this->input_lists.end())
        {
            Input_Item* p_item = &(it->second);
            this->readvalue_items.push_back(p_item);
            if(p_item->is_read())
            {
                std::string warningstr = "The parameter " + p_item->label + " has been read twice.";
                ModuleBase::WARNING_QUIT("ReadInput", warningstr);
            }
            // qianrui delete '/' 2024-07-10, because path has '/' head.
            read_information(ifs, p_item->str_values, "#!");
        }
        else
        {
            if (word[0] != '#' && word[0] != '/' && word[0] != '!')
            {
                std::cout << " THE PARAMETER NAME '" << word << "' IS NOT USED!" << std::endl;
                ModuleBase::WARNING_QUIT("Input",
                                         "Bad parameter, please check the "
                                         "input parameters in file INPUT",
                                         1);
            }
            ifs.ignore(150, '\n');
        }

        ifs.rdstate();
        if (ifs.eof())
        {
            break;
        }
        else if (ifs.bad())
        {
            std::cout << " Bad input parameters. " << std::endl;
            exit(1);
        }
        else if (ifs.fail())
        {
            std::cout << " word = " << word << std::endl;
            std::cout << " Fail to read parameters. " << std::endl;
            ifs.clear();
            exit(1);
        }
    }

    // 1) read the value of the parameters
    for (auto& readvalue_item: this->readvalue_items)
    {
        readvalue_item->read_value(*readvalue_item, param);
    }

    // 2) count the number of atom types from STRU file
    if (this->check_ntype_flag) {
        check_ntype(param.input.stru_file, param.input.ntype);
}

    // 3) reset this value when some conditions are met
    //    e.g. if (calulation_type == "nscf") then set "init_chg" to "file".
    for (auto& input_item: this->input_lists)
    {
        Input_Item* resetvalue_item = &(input_item.second);
        if (resetvalue_item->reset_value != nullptr) {
            resetvalue_item->reset_value(*resetvalue_item, param);
}
    }

    // 4) check the value of the parameters
    for (auto& input_item: this->input_lists)
    {
        Input_Item* checkvalue_item = &(input_item.second);
        if (checkvalue_item->check_value != nullptr) {
            checkvalue_item->check_value(*checkvalue_item, param);
}
    }
}

void ReadInput::write_txt_input(const Parameter& param, const std::string& filename)
{
    ModuleBase::TITLE("ReadInput", "write_txt_input");
    std::ofstream ofs(filename.c_str(), std::ios::out);
    ofs << "INPUT_PARAMETERS" << std::endl;
    ofs << std::setiosflags(std::ios::left);

    ofs << "#Parameters (1.General)" << std::endl;
    for (auto& item: this->input_lists)
    {
        Input_Item* p_item = &(item.second);
        if (p_item->get_final_value == nullptr) {
            continue;
}
        p_item->get_final_value(*p_item, param);
        if (p_item->label == "ecutwfc")
        {
            ofs << "\n#Parameters (2.PW)" << std::endl;
        }
        else if (p_item->label == "method_sto")
        {
            ofs << "\n#Parameters (3.Stochastic DFT)" << std::endl;
        }
        else if (p_item->label == "ks_solver")
        {
            ofs << "\n#Parameters (4.Relaxation)" << std::endl;
        }
        else if (p_item->label == "basis_type")
        {
            ofs << "\n#Parameters (5.LCAO)" << std::endl;
        }
        else if (p_item->label == "smearing_method")
        {
            ofs << "\n#Parameters (6.Smearing)" << std::endl;
        }
        else if (p_item->label == "mixing_type")
        {
            ofs << "\n#Parameters (7.Charge Mixing)" << std::endl;
        }
        else if (p_item->label == "dos_emin_ev")
        {
            ofs << "\n#Parameters (8.DOS)" << std::endl;
        }
        else if (p_item->label == "md_type")
        {
            ofs << "\n#Parameters (9.Molecular dynamics)" << std::endl;
        }
        else if (p_item->label == "efield_flag")
        {
            ofs << "\n#Parameters (10.Electric field and dipole correction)" << std::endl;
        }
        else if (p_item->label == "gate_flag")
        {
            ofs << "\n#Parameters (11.Gate field)" << std::endl;
        }
        else if (p_item->label == "out_alllog")
        {
            ofs << "\n#Parameters (12.Test)" << std::endl;
        }
        else if (p_item->label == "vdw_method")
        {
            ofs << "\n#Parameters (13.vdW Correction)" << std::endl;
        }
        else if (p_item->label == "exx_hybrid_alpha")
        {
            ofs << "\n#Parameters (14.exx)" << std::endl;
        }
        else if (p_item->label == "td_force_dt")
        {
            ofs << "\n#Parameters (16.tddft)" << std::endl;
        }
        else if (p_item->label == "berry_phase")
        {
            ofs << "\n#Parameters (17.berry_wannier)" << std::endl;
        }
        else if (p_item->label == "imp_sol")
        {
            ofs << "\n#Parameters (18.implicit_solvation)" << std::endl;
        }
        else if (p_item->label == "of_kinetic")
        {
            ofs << "\n#Parameters (19.orbital free density functional theory)" << std::endl;
        }
        else if (p_item->label == "dft_plus_u")
        {
            ofs << "\n#Parameters (20.dft+u)" << std::endl;
        }
        else if (p_item->label == "bessel_nao_ecut")
        {
            ofs << "\n#Parameters (21.spherical bessel)" << std::endl;
        }
        else if (p_item->label == "sc_mag_switch")
        {
            ofs << "\n#Parameters (22.non-collinear spin-constrained DFT)" << std::endl;
        }
        else if (p_item->label == "qo_switch")
        {
            ofs << "\n#Parameters (23.Quasiatomic Orbital analysis)" << std::endl;
        }
        else if (p_item->label == "pexsi_npole")
        {
            ofs << "\n#Parameters (24.PEXSI)" << std::endl;
        }
        ModuleBase::GlobalFunc::OUTP(ofs, p_item->label, p_item->final_value.str(), p_item->annotation);
    }
}

void ReadInput::check_ntype(const std::string& fn, int& param_ntype)
{
    std::ifstream ifa(fn.c_str(), std::ios::in);
    if (!ifa)
    {
        GlobalV::ofs_warning << fn;
        ModuleBase::WARNING_QUIT("Input::count_ntype", "Can not find the file containing atom positions.!");
    }

    int ntype_stru = 0;
    std::string temp;
    if (ModuleBase::GlobalFunc::SCAN_BEGIN(ifa, "ATOMIC_SPECIES"))
    {
        while (true)
        {
            ModuleBase::GlobalFunc::READ_VALUE(ifa, temp);
            if (temp == "LATTICE_CONSTANT" || temp == "NUMERICAL_ORBITAL" || temp == "NUMERICAL_DESCRIPTOR"
                || temp == "PAW_FILES" || ifa.eof())
            {
                break;
            }
            else if (isalpha(temp[0]))
            {
                ntype_stru += 1;
            }
        }
    }

    if (param_ntype == 0)
    {
        param_ntype = ntype_stru;
        GlobalV::ofs_running << " 'ntype' is no longer required in INPUT, and "
                                "it will be ignored. "
                             << std::endl;
    }
    else if (param_ntype != ntype_stru)
    {
        ModuleBase::WARNING_QUIT("ReadInput",
                                 "The ntype in INPUT is not equal to the ntype "
                                 "counted in STRU, check it.");
    }
    if (param_ntype <= 0)
    {
        ModuleBase::WARNING_QUIT("ReadInput", "ntype should be greater than 0.");
    }
    else {
        GlobalV::ofs_running << " 'ntype' is automatically set to " << param_ntype << std::endl;
}
}

void ReadInput::add_item(const Input_Item& item)
{
    // only rank 0 read the input file
    // So only rank 0 add the item to the input list
    if (this->rank == 0)
    {
        this->input_lists.push_back(make_pair(item.label, item));
    }
}
bool find_str(const std::vector<std::string>& strings, const std::string& strToFind)
{
    auto it = std::find(strings.begin(), strings.end(), strToFind);
    return it != strings.end();
};
} // namespace ModuleIO
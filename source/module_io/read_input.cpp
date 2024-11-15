#include "read_input.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <cassert>
#include "module_base/formatter.h"
#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/tool_quit.h"
#include "module_base/tool_title.h"
namespace ModuleIO
{

std::string longstring(const std::vector<std::string>& words)
{
    return FmtCore::join(" ", words);
}

bool assume_as_boolean(const std::string& val)
{
    const std::string val_ = FmtCore::lower(val);

    const std::array<std::string, 7> t_ = {"true", "1", "t", "yes", "y", "on", ".true."};
    const std::array<std::string, 7> f_ = {"false", "0", "f", "no", "n", "off", ".false."};
    // This will work because std::array<T, N>::size() is a constexpr function
    // Ouch it is of C++17 standard...
    // static_assert(t_.size() == f_.size(), "t_ and f_ must have the same lengths");
#ifdef __DEBUG // C++11 can do this
    assert(t_.size() == f_.size());
#endif

    if (std::find(t_.begin(), t_.end(), val_) != t_.end())
    {
        return true;
    }
    else if (std::find(f_.begin(), f_.end(), val_) != f_.end())
    {
        return false;
    }
    else
    {
        std::string warnmsg = "Bad boolean parameter ";
        warnmsg.append(val);
        warnmsg.append(", please check the input parameters in file INPUT");
        ModuleBase::WARNING_QUIT("Input", warnmsg);
    }
}

std::string to_dir(const std::string& str)
{
    std::string str_dir = str;
    if (str_dir.empty())
    {
        return "./";
    }
    else if (str_dir.back() != '/')
    {
        str_dir += "/";
    }
    
    return str_dir;
}

void read_information(std::ifstream& ifs, std::vector<std::string>& output, const std::string& delimiters)
{
    std::string line;
    getline(ifs, line);

    std::istringstream iss(line);
    std::string word;
    while (iss >> word)
    {
        if (delimiters.find(word[0]) != std::string::npos)
        {
            break;
        }
        output.push_back(word);
    }
}

bool ReadInput::check_mode = false;

ReadInput::ReadInput(const int& rank)
{
    this->rank = rank;
    
    // add items
    this->item_system();
    this->item_elec_stru();
    this->item_relax();
    this->item_md();
    this->item_ofdft();
    this->item_sdft();
    this->item_deepks();
    this->item_rt_tddft();
    this->item_lr_tddft();
    this->item_output();
    this->item_postprocess();
    this->item_model();
    this->item_exx();
    this->item_dftu();
    this->item_others();

    // set globalv functions
    this->set_globalv_bcast();
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

    std::string word, word1;
    int ierr = 0;

    // ifs >> std::setiosflags(ios::uppercase);
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n');
        if (word == "INPUT_PARAMETERS")
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
        if (ifs.eof()) { break; }
        word = FmtCore::lower(word1);
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

    // 2) reset this value when some conditions are met
    //    e.g. if (calulation_type == "nscf") then set "init_chg" to "file".
    for (auto& input_item: this->input_lists)
    {
        Input_Item* resetvalue_item = &(input_item.second);
        if (resetvalue_item->reset_value != nullptr) {
            resetvalue_item->reset_value(*resetvalue_item, param);
        }
    }
    this->set_globalv(param);

    // 3) count the number of atom types from STRU file
    if (this->check_ntype_flag)
    {
        check_ntype(param.globalv.global_in_stru, param.input.ntype);
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

    ofs << "#Parameters (1.System)" << std::endl;
    for (auto& item: this->input_lists)
    {
        Input_Item* p_item = &(item.second);
        if (p_item->get_final_value == nullptr) {
            continue;
}
        p_item->get_final_value(*p_item, param);
        if (p_item->label == "ks_solver")
        {
            ofs << "\n#Parameters (2.Electronic structure)" << std::endl;
        }
        else if (p_item->label == "nb2d")
        {
            ofs << "\n#Parameters (3.LCAO)" << std::endl;
        }
        else if (p_item->label == "relax_method")
        {
            ofs << "\n#Parameters (4.Relaxation)" << std::endl;
        }
        else if (p_item->label == "md_type")
        {
            ofs << "\n#Parameters (5.Molecular dynamics)" << std::endl;
        }
        else if (p_item->label == "of_kinetic")
        {
            ofs << "\n#Parameters (6.orbital free density functional theory)" << std::endl;
        }
        else if (p_item->label == "method_sto")
        {
            ofs << "\n#Parameters (7.Stochastic DFT)" << std::endl;
        }
        else if (p_item->label == "deepks_out_labels")
        {
            ofs << "\n#Parameters (8.DeepKS)" << std::endl;
        }
        else if (p_item->label == "td_force_dt")
        {
            ofs << "\n#Parameters (9.rt-tddft)" << std::endl;
        }
        else if (p_item->label == "lr_nstates")
        {
            ofs << "\n#Parameters (10.lr-tddft)" << std::endl;
        }
        else if (p_item->label == "out_stru")
        {
            ofs << "\n#Parameters (11.Output)" << std::endl;
        }
        else if (p_item->label == "dos_emin_ev")
        {
            ofs << "\n#Parameters (12.Postprocess)" << std::endl;
        }
        else if (p_item->label == "efield_flag")
        {
            ofs << "\n#Parameters (13.Model)" << std::endl;
        }
        else if (p_item->label == "vdw_method")
        {
            ofs << "\n#Parameters (14.vdW Correction)" << std::endl;
        }
        else if (p_item->label == "exx_hybrid_alpha")
        {
            ofs << "\n#Parameters (15.exx)" << std::endl;
        }
        else if (p_item->label == "dft_plus_u")
        {
            ofs << "\n#Parameters (16.dft+u)" << std::endl;
        }
        else if (p_item->label == "sc_mag_switch")
        {
            ofs << "\n#Parameters (17.non-collinear spin-constrained DFT)" << std::endl;
        }
        else if (p_item->label == "qo_switch")
        {
            ofs << "\n#Parameters (18.Quasiatomic Orbital analysis)" << std::endl;
        }
        else if (p_item->label == "pexsi_npole")
        {
            ofs << "\n#Parameters (19.PEXSI)" << std::endl;
        }
        else if (p_item->label == "out_alllog")
        {
            ofs << "\n#Parameters (20.Test)" << std::endl;
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
        ModuleBase::WARNING_QUIT("ReadInput::check_ntype", "Can not find the file: " + fn);
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

int ReadInput::current_md_step(const std::string& file_dir)
{
    std::stringstream ssc;
    ssc << file_dir << "Restart_md.dat";
    std::ifstream file(ssc.str().c_str());

    if (!file)
    {
        ModuleBase::WARNING_QUIT("current_md_step", "no Restart_md.dat");
    }

    int md_step;
    file >> md_step;
    file.close();

    return md_step;
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

std::string nofound_str(std::vector<std::string> init_chgs, const std::string& str)
{
    std::string warningstr = "The parameter ";
    warningstr.append(str);
    warningstr.append(" must be ");
    for(int i = 0; i < init_chgs.size(); i++)
    {
        warningstr.append("'");
        warningstr.append(init_chgs[i]);
        warningstr.append("'");
        if(i < init_chgs.size() - 2)
        {
            warningstr.append(", ");
        }
        else if(i == init_chgs.size() - 2)
        {
            warningstr.append(" or ");
        }
    }
    warningstr.append("!");

    return warningstr;
}

} // namespace ModuleIO
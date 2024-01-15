#include "module_io/parameter_pool.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "module_base/constants.h"
#include "module_base/global_file.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "module_io/input.h"
#include "module_md/md_para.h"

/**
 * @param input_parameters Save all input parameters
 * @param default_parametes_type Save   the names and types of all parameters
 */
namespace ModuleIO
{

std::map<std::string, InputParameter> input_parameters;
std::map<std::string, std::string> default_parametes_type;
std::map<std::string, InputParameter> default_parametes_value;

// Conut how many types of atoms are listed in STRU
int count_ntype(const std::string& fn)
{
    // Only RANK0 core can reach here, because this function is called during Input::Read.
    assert(GlobalV::MY_RANK == 0);

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
            if (temp == "LATTICE_CONSTANT" || temp == "NUMERICAL_ORBITAL" || temp == "NUMERICAL_DESCRIPTOR" || temp == "PAW_FILES"
                || ifa.eof())
            {
                break;
            }
            else if (isalpha(temp[0]))
            {
                ntype_stru += 1;
            }
        }
    }
    return ntype_stru;
}

/**
 * @brief New param init function. First, the default parameter types are read,
 *        then the default parameter values file is read and populated,
 *        and the values in the map are passed to input.h
 *        then do the same to input parameter values
 * @param default_type_path parameter default type file path
 * @param input_value_path parameter default value file path
 * @param input_value_path parameter input value file path
 */
bool Init(const std::string& default_type_path,
          const std::string& default_value_path,
          const std::string& input_value_path)
{
    ModuleBase::timer::tick("Input", "Init");
    default_parametes_reader(default_type_path, default_parametes_type);
    input_parameters_get(default_value_path, default_parametes_value);
    input_parameters_set(default_parametes_value);

    input_parameters_get(input_value_path, input_parameters);
    input_parameters_set(input_parameters);
}

/**
 * @brief Converts the string sa to a minor character and stores it in the string sb
 * @param sa Pointer to the string to be converted
 * @param sb Store a pointer to the converted string, which needs to be pre-allocated with enough space
 */
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
} // namespace Modlevoid strtolower(char*sa,char*sb)

/**
 * @brief Reads the default parameters from the specified file and saves them to the global variable
 *        default_parametes_type
 * @param fn Specifies the path to the file
 * @return true Read successfully
 * @return false Read failure
 */
bool default_parametes_reader(const std::string& fn, std::map<std::string, std::string>& default_parametes_type)
{
    std::ifstream inputFile(fn.c_str());
    if (inputFile.is_open())
    {
        std::string word1, word2;
        while (inputFile >> word1 >> word2)
        {
            // default_parametes_type[word1] = word2.c_str();
            default_parametes_type.insert(std::pair<std::string, std::string>(word1, word2));
        }
        // Close file
        inputFile.close();
    }
    else
    {
        std::cout << "Cannot open file !" << std::endl;
    }
}
/**
 * @brief This function is used to read the input parameter file and store it as a key-value pair
 * @param fn Enter the path to the parameter file
 */
bool input_parameters_get(const std::string& fn, std::map<std::string, InputParameter>& input)
{
    // The module title information is displayed
    ModuleBase::TITLE("Input", "Read");
    // If it is not the primary node, return false
    if (GlobalV::MY_RANK != 0)
        return false;

    // Open the input parameter file
    std::ifstream ifs(fn.c_str(), std::ios::in); // "in_datas/input_parameters"
    // If the opening fails, an error message is printed and false is returned
    if (!ifs)
    {
        std::cout << " Can't find the INPUT file." << std::endl;
        return false;
    }
    ifs.clear();
    ifs.seekg(0);
    char word[80], word1[80];
    int ierr = 0;

    // Read file contents
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> word;
        ifs.ignore(150, '\n'); // Ignore end of line
        if (strcmp(word, "INPUT_PARAMETERS") == 0)
        {
            ierr = 1;
            break;
        }
        ifs.rdstate();
    }
    // If ierr is 0, the word "INPUT_PARAMETERS" is not found, and an error message is printed with false
    if (ierr == 0)
    {
        std::cout << " Error parameter list." << std::endl;
        return false; // return error : false
    }
    ifs.rdstate();

    bool param_bool;
    int param_int;
    double param_double;
    SimpleString param_string;
    SimpleVector<int> param_vector_int;
    SimpleVector<double> param_vector_double;
    void* param_value;

    // Parameter values are read and stored as key-value pairs
    while (ifs.good())
    {
        ifs >> word1;
        if (ifs.eof())
            break;
        strtolower(word1, word);

        //----------------------------------------------------------
        // main parameters
        //----------------------------------------------------------
        // The parameter values are read according to the parameter type and stored as key-value pairs
        std::string param_type = default_parametes_type[word];
        InputParameter input_param;
        if (strcmp(param_type.c_str(), "int") == 0)
        {
            INPUT.read_value(ifs, param_int);
            input_param.type = INT;
            input_param.set((void*)(&param_int));
        }
        else if (strcmp(param_type.c_str(), "bool") == 0)
        {

            INPUT.read_bool(ifs, param_bool);
            input_param.type = BOOL;
            input_param.set((void*)(&param_bool));
        }
        else if (strcmp(param_type.c_str(), "double") == 0)
        {
            INPUT.read_value(ifs, param_double);
            input_param.type = DOUBLE;
            input_param.set((void*)(&param_double));
        }
        else if (strcmp(param_type.c_str(), "string") == 0)
        {
            std::string s;
            INPUT.read_value(ifs, s);
            param_string = s.c_str();
            input_param.type = STRING;
            input_param.set((void*)(&param_string));
        }
        else if (strcmp(param_type.c_str(), "vector_int") == 0)
        {
            int tmp;
            std::string s;
            std::getline(ifs, s);
            std::stringstream ss(s);
            while ((ss >> tmp))
            {
                param_vector_int.push_back(tmp);
            }
            input_param.type = VECTOR_I;
            input_param.set((void*)(&param_vector_int));
        }
        else if (strcmp(param_type.c_str(), "vector_double") == 0)
        {
            double tmp;
            std::string s;
            std::getline(ifs, s);
            std::stringstream ss(s);
            while ((ss >> tmp))
            {
                param_vector_double.push_back(tmp);
            }
            input_param.type = VECTOR_D;
            input_param.set((void*)(&param_vector_double));
        }
        else
        {
            // xiaohui add 2015-09-15
            if (word[0] != '#' && word[0] != '/')
            {
                INPUT.input_error = 1;
                std::cout << " THE PARAMETER NAME '" << word << "' IS NOT USED!" << std::endl;
            }
            // mohan screen this 2012-06-30
            //			std::cout << " THE PARAMETER NAME '" << word
            //			<< "' IS NOT USED!" << std::endl;
            ifs.ignore(150, '\n');
        }
        ifs.rdstate();

        /*if(gamma_only == 1)
        {
           gamma_only_local = 1;      //pengfei 2014-10-15
           gamma_only = 0;
           std::cout << "gamma_only_local = " << gamma_only_local <<std::endl;
        }*/

        if (ifs.eof() != 0)
        {
            break;
        }
        else if (ifs.bad() != 0)
        {
            std::cout << " Bad input parameters. " << std::endl;
            return false;
        }
        else if (ifs.fail() != 0)
        {
            std::cout << " word = " << word << std::endl;
            std::cout << " Fail to read parameters. " << std::endl;
            ifs.clear();
            return false;
        }
        else if (ifs.good() == 0)
        {
            break;
        }
        input[word] = input_param;
    }

    if (INPUT.stru_file == "")
    {
        INPUT.stru_file = "STRU";
    }
    double ntype_stru = count_ntype(INPUT.stru_file);
    if (INPUT.ntype == 0)
    {
        INPUT.ntype = ntype_stru;
        GlobalV::ofs_running << "ntype in INPUT is 0, and it is automatically set to " << INPUT.ntype
                             << " according to STRU" << std::endl;
    }
    else if (INPUT.ntype != ntype_stru)
    {
        ModuleBase::WARNING_QUIT("Input", "The ntype in INPUT is not equal to the ntype counted in STRU, check it.");
    }

    return true;
}

bool input_parameters_set(std::map<std::string, InputParameter> input_parameters)
{
    if (input_parameters.count("nupdown") != 0)
    {
        INPUT.nupdown = *static_cast<double*>(input_parameters["nupdown"].get());
    }
    else if (input_parameters.count("suffix") != 0)
    {
        INPUT.suffix = static_cast<SimpleString*>(input_parameters["suffix"].get())->c_str();
    }
    else if (input_parameters.count("stru_file") != 0)
    {
        INPUT.stru_file = static_cast<SimpleString*>(input_parameters["stru_file"].get())->c_str();
    }
    else if (input_parameters.count("pseudo_dir") != 0)
    {
        INPUT.pseudo_dir = static_cast<SimpleString*>(input_parameters["pseudo_dir"].get())->c_str();
    }
    else if (input_parameters.count("orbital_dir") != 0)
    {
        INPUT.orbital_dir = static_cast<SimpleString*>(input_parameters["orbital_dir"].get())->c_str();
    }
    else if (input_parameters.count("read_file_dir") != 0)
    {
        INPUT.read_file_dir = static_cast<SimpleString*>(input_parameters["read_file_dir"].get())->c_str();
    }
    else if (input_parameters.count("kpoint_file") != 0)
    {
        INPUT.kpoint_file = static_cast<SimpleString*>(input_parameters["kpoint_file"].get())->c_str();
    }
    else if (input_parameters.count("wannier_card") != 0)
    {
        INPUT.wannier_card = static_cast<SimpleString*>(input_parameters["wannier_card"].get())->c_str();
    }
    else if (input_parameters.count("latname") != 0)
    {
        INPUT.latname = static_cast<SimpleString*>(input_parameters["latname"].get())->c_str();
    }
    else if (input_parameters.count("calculation") != 0)
    {
        INPUT.calculation = static_cast<SimpleString*>(input_parameters["calculation"].get())->c_str();
    }
    else if (input_parameters.count("esolver_type") != 0)
    {
        INPUT.esolver_type = static_cast<SimpleString*>(input_parameters["esolver_type"].get())->c_str();
    }
    else if (input_parameters.count("pseudo_rcut") != 0)
    {
        INPUT.pseudo_rcut = *static_cast<double*>(input_parameters["pseudo_rcut"].get());
    }
    else if (input_parameters.count("pseudo_mesh") != 0)
    {
        INPUT.pseudo_mesh = *static_cast<bool*>(input_parameters["pseudo_mesh"].get());
    }
    else if (input_parameters.count("ntype") != 0)
    {
        INPUT.ntype = *static_cast<int*>(input_parameters["ntype"].get());
    }
    else if (input_parameters.count("nbands") != 0)
    {
        INPUT.nbands = *static_cast<int*>(input_parameters["nbands"].get());
    }
    else if (input_parameters.count("nbands_istate") != 0)
    {
        INPUT.nbands_istate = *static_cast<int*>(input_parameters["nbands_istate"].get());
    }
    else if (input_parameters.count("pw_seed") != 0)
    {
        INPUT.pw_seed = *static_cast<int*>(input_parameters["pw_seed"].get());
    }
    else if (input_parameters.count("init_vel") != 0)
    {
        INPUT.init_vel = *static_cast<bool*>(input_parameters["init_vel"].get());
    }
    else if (input_parameters.count("ref_cell_factor") != 0)
    {
        INPUT.ref_cell_factor = *static_cast<double*>(input_parameters["ref_cell_factor"].get());
    }
    else if (input_parameters.count("symmetry") != 0)
    {
        INPUT.symmetry = *static_cast<int*>(input_parameters["symmetry"].get());
    }
    else if (input_parameters.count("symmetry_prec") != 0)
    {
        INPUT.symmetry_prec = *static_cast<double*>(input_parameters["symmetry_prec"].get());
    }
    else if (input_parameters.count("kpar") != 0)
    {
        INPUT.kpar = *static_cast<int*>(input_parameters["kpar"].get());
    }
    else if (input_parameters.count("berry_phase") != 0)
    {
        INPUT.berry_phase = *static_cast<bool*>(input_parameters["berry_phase"].get());
    }
    else if (input_parameters.count("gdir") != 0)
    {
        INPUT.gdir = *static_cast<int*>(input_parameters["gdir"].get());
    }
    else if (input_parameters.count("kspacing") != 0)
    {
        SimpleVector<double> vec_D = *static_cast<SimpleVector<double>*>(input_parameters["kspacing"].get());
        for (int i = 0; i <= vec_D.size(); i++)
        {
            INPUT.kspacing[i] = vec_D[i];
        }
        if (vec_D.size() == 0 || vec_D.size() == 2)
        {
            std::cout << "kspacing can only accept one or three double values." << std::endl;
            // ifs.setstate(std::ios::failbit);
        }
        // if only read one value, set all to kspacing[0]
        if (vec_D.size() == 1)
        {
            INPUT.kspacing[1] = INPUT.kspacing[0];
            INPUT.kspacing[2] = INPUT.kspacing[0];
        }
    }
    else if (input_parameters.count("min_dist_coef") != 0)
    {
        INPUT.min_dist_coef = *static_cast<double*>(input_parameters["min_dist_coef"].get());
    }
    else if (input_parameters.count("towannier90") != 0)
    {
        INPUT.towannier90 = *static_cast<bool*>(input_parameters["towannier90"].get());
    }
    else if (input_parameters.count("nnkpfile") != 0)
    {
        INPUT.nnkpfile = static_cast<SimpleString*>(input_parameters["nnkpfile"].get())->c_str();
    }
    else if (input_parameters.count("wannier_spin") != 0)
    {
        INPUT.wannier_spin = static_cast<SimpleString*>(input_parameters["wannier_spin"].get())->c_str();
    }
    else if (input_parameters.count("wannier_method") != 0)
    {
        INPUT.wannier_method = *static_cast<int*>(input_parameters["wannier_method"].get());
    }
    else if (input_parameters.count("out_wannier_mmn") != 0)
    {
        INPUT.out_wannier_mmn = *static_cast<bool*>(input_parameters["out_wannier_mmn"].get());
    }
    else if (input_parameters.count("out_wannier_amn") != 0)
    {
        INPUT.out_wannier_amn = *static_cast<bool*>(input_parameters["out_wannier_amn"].get());
    }
    else if (input_parameters.count("out_wannier_unk") != 0)
    {
        INPUT.out_wannier_unk = *static_cast<bool*>(input_parameters["out_wannier_unk"].get());
    }
    else if (input_parameters.count("out_wannier_eig") != 0)
    {
        INPUT.out_wannier_eig = *static_cast<bool*>(input_parameters["out_wannier_eig"].get());
    }
    else if (input_parameters.count("out_wannier_wvfn_formatted") != 0)
    {
        INPUT.out_wannier_wvfn_formatted = *static_cast<bool*>(input_parameters["out_wannier_wvfn_formatted"].get());
    }
    else if (input_parameters.count("nche_sto") != 0)
    {
        INPUT.nche_sto = *static_cast<int*>(input_parameters["nche_sto"].get());
    }
    else if (input_parameters.count("nbands_sto") != 0)
    {
        INPUT.nbands_sto = *static_cast<int*>(input_parameters["nbands_sto"].get());
    }
    else if (input_parameters.count("nbndsto_str") != 0)
    {
        INPUT.nbndsto_str = static_cast<SimpleString*>(input_parameters["nbndsto_str"].get())->c_str();
    }
    else if (input_parameters.count("seed_sto") != 0)
    {
        INPUT.seed_sto = *static_cast<int*>(input_parameters["seed_sto"].get());
    }
    else if (input_parameters.count("initsto_ecut") != 0)
    {
        INPUT.initsto_ecut = *static_cast<double*>(input_parameters["initsto_ecut"].get());
    }
    else if (input_parameters.count("emax_sto") != 0)
    {
        INPUT.emax_sto = *static_cast<double*>(input_parameters["emax_sto"].get());
    }
    else if (input_parameters.count("emin_sto") != 0)
    {
        INPUT.emin_sto = *static_cast<double*>(input_parameters["emin_sto"].get());
    }
    else if (input_parameters.count("bndpar") != 0)
    {
        INPUT.bndpar = *static_cast<int*>(input_parameters["bndpar"].get());
    }
    else if (input_parameters.count("initsto_freq") != 0)
    {
        INPUT.initsto_freq = *static_cast<int*>(input_parameters["initsto_freq"].get());
    }
    else if (input_parameters.count("method_sto") != 0)
    {
        INPUT.method_sto = *static_cast<int*>(input_parameters["method_sto"].get());
    }
    else if (input_parameters.count("npart_sto") != 0)
    {
        INPUT.npart_sto = *static_cast<int*>(input_parameters["npart_sto"].get());
    }
    else if (input_parameters.count("cal_cond") != 0)
    {
        INPUT.cal_cond = *static_cast<bool*>(input_parameters["cal_cond"].get());
    }
    else if (input_parameters.count("cond_che_thr") != 0)
    {
        INPUT.cond_che_thr = *static_cast<double*>(input_parameters["cond_che_thr"].get());
    }
    else if (input_parameters.count("cond_dw") != 0)
    {
        INPUT.cond_dw = *static_cast<double*>(input_parameters["cond_dw"].get());
    }
    else if (input_parameters.count("cond_wcut") != 0)
    {
        INPUT.cond_wcut = *static_cast<double*>(input_parameters["cond_wcut"].get());
    }
    else if (input_parameters.count("cond_dt") != 0)
    {
        INPUT.cond_dt = *static_cast<int*>(input_parameters["cond_dt"].get());
    }
    else if (input_parameters.count("cond_dtbatch") != 0)
    {
        INPUT.cond_dtbatch = *static_cast<int*>(input_parameters["cond_dtbatch"].get());
    }
    else if (input_parameters.count("cond_smear") != 0)
    {
        INPUT.cond_smear = *static_cast<int*>(input_parameters["cond_smear"].get());
    }
    else if (input_parameters.count("cond_fwhm") != 0)
    {
        INPUT.cond_fwhm = *static_cast<double*>(input_parameters["cond_fwhm"].get());
    }
    else if (input_parameters.count("cond_nonlocal") != 0)
    {
        INPUT.cond_nonlocal = *static_cast<bool*>(input_parameters["cond_nonlocal"].get());
    }
    else if (input_parameters.count("dft_functional") != 0)
    {
        INPUT.dft_functional = static_cast<SimpleString*>(input_parameters["dft_functional"].get())->c_str();
    }
    else if (input_parameters.count("xc_temperature") != 0)
    {
        INPUT.xc_temperature = *static_cast<double*>(input_parameters["xc_temperature"].get());
    }
    else if (input_parameters.count("nspin") != 0)
    {
        INPUT.nspin = *static_cast<int*>(input_parameters["nspin"].get());
    }
    else if (input_parameters.count("nupdown") != 0)
    {
        INPUT.nupdown = *static_cast<double*>(input_parameters["nupdown"].get());
    }
    else if (input_parameters.count("nelec") != 0)
    {
        INPUT.nelec = *static_cast<double*>(input_parameters["nelec"].get());
    }
    else if (input_parameters.count("lmaxmax") != 0)
    {
        INPUT.lmaxmax = *static_cast<int*>(input_parameters["lmaxmax"].get());
    }
    else if (input_parameters.count("basis_type") != 0)
    {
        INPUT.basis_type = static_cast<SimpleString*>(input_parameters["basis_type"].get())->c_str();
    }
    else if (input_parameters.count("ks_solver") != 0)
    {
        INPUT.ks_solver = static_cast<SimpleString*>(input_parameters["ks_solver"].get())->c_str();
    }
    else if (input_parameters.count("cal_force") != 0)
    {
        INPUT.cal_force = *static_cast<bool*>(input_parameters["cal_force"].get());
    }
    else if (input_parameters.count("force_thr") != 0)
    {
        INPUT.force_thr = *static_cast<double*>(input_parameters["force_thr"].get());
    }
    else if (input_parameters.count("force_thr_ev2") != 0)
    {
        INPUT.force_thr_ev2 = *static_cast<double*>(input_parameters["force_thr_ev2"].get());
    }
    else if (input_parameters.count("stress_thr") != 0)
    {
        INPUT.stress_thr = *static_cast<double*>(input_parameters["stress_thr"].get());
    }
    else if (input_parameters.count("press1") != 0)
    {
        INPUT.press1 = *static_cast<double*>(input_parameters["press1"].get());
    }
    else if (input_parameters.count("press2") != 0)
    {
        INPUT.press2 = *static_cast<double*>(input_parameters["press2"].get());
    }
    else if (input_parameters.count("press3") != 0)
    {
        INPUT.press3 = *static_cast<double*>(input_parameters["press3"].get());
    }
    else if (input_parameters.count("cal_stress") != 0)
    {
        INPUT.cal_stress = *static_cast<bool*>(input_parameters["cal_stress"].get());
    }
    else if (input_parameters.count("fixed_axes") != 0)
    {
        INPUT.fixed_axes = static_cast<SimpleString*>(input_parameters["fixed_axes"].get())->c_str();
    }
    else if (input_parameters.count("fixed_ibrav") != 0)
    {
        INPUT.fixed_ibrav = *static_cast<bool*>(input_parameters["fixed_ibrav"].get());
    }
    else if (input_parameters.count("fixed_atoms") != 0)
    {
        INPUT.fixed_atoms = *static_cast<bool*>(input_parameters["fixed_atoms"].get());
    }
    else if (input_parameters.count("relax_method") != 0)
    {
        INPUT.relax_method = static_cast<SimpleString*>(input_parameters["relax_method"].get())->c_str();
    }
    else if (input_parameters.count("relax_new") != 0)
    {
        INPUT.relax_new = *static_cast<bool*>(input_parameters["relax_new"].get());
    }
    else if (input_parameters.count("relax_cg_thr") != 0)
    {
        INPUT.relax_cg_thr = *static_cast<double*>(input_parameters["relax_cg_thr"].get());
    }
    else if (input_parameters.count("relax_bfgs_w1") != 0)
    {
        INPUT.relax_bfgs_w1 = *static_cast<double*>(input_parameters["relax_bfgs_w1"].get());
    }
    else if (input_parameters.count("relax_bfgs_w2") != 0)
    {
        INPUT.relax_bfgs_w2 = *static_cast<double*>(input_parameters["relax_bfgs_w2"].get());
    }
    else if (input_parameters.count("relax_bfgs_rmax") != 0)
    {
        INPUT.relax_bfgs_rmax = *static_cast<double*>(input_parameters["relax_bfgs_rmax"].get());
    }
    else if (input_parameters.count("relax_bfgs_rmin") != 0)
    {
        INPUT.relax_bfgs_rmin = *static_cast<double*>(input_parameters["relax_bfgs_rmin"].get());
    }
    else if (input_parameters.count("relax_bfgs_init") != 0)
    {
        INPUT.relax_bfgs_init = *static_cast<double*>(input_parameters["relax_bfgs_init"].get());
    }
    else if (input_parameters.count("relax_scale_force") != 0)
    {
        INPUT.relax_scale_force = *static_cast<double*>(input_parameters["relax_scale_force"].get());
    }
    else if (input_parameters.count("gamma_only") != 0)
    {
        INPUT.gamma_only = *static_cast<bool*>(input_parameters["gamma_only"].get());
    }
    else if (input_parameters.count("gamma_only_local") != 0)
    {
        INPUT.gamma_only_local = *static_cast<bool*>(input_parameters["gamma_only_local"].get());
    }
    else if (input_parameters.count("fft_mode") != 0)
    {
        INPUT.fft_mode = *static_cast<int*>(input_parameters["fft_mode"].get());
    }
    else if (input_parameters.count("ecutwfc") != 0)
    {
        INPUT.ecutwfc = *static_cast<double*>(input_parameters["ecutwfc"].get());
    }
    else if (input_parameters.count("ecutrho") != 0)
    {
        INPUT.ecutrho = *static_cast<double*>(input_parameters["ecutrho"].get());
    }
    else if (input_parameters.count("ncx") != 0)
    {
        INPUT.ncx = *static_cast<int*>(input_parameters["ncx"].get());
    }
    else if (input_parameters.count("ncy") != 0)
    {
        INPUT.ncy = *static_cast<int*>(input_parameters["ncy"].get());
    }
    else if (input_parameters.count("ncz") != 0)
    {
        INPUT.ncz = *static_cast<int*>(input_parameters["ncz"].get());
    }
    else if (input_parameters.count("nx") != 0)
    {
        INPUT.nx = *static_cast<int*>(input_parameters["nx"].get());
    }
    else if (input_parameters.count("ny") != 0)
    {
        INPUT.ny = *static_cast<int*>(input_parameters["ny"].get());
    }
    else if (input_parameters.count("nz") != 0)
    {
        INPUT.nz = *static_cast<int*>(input_parameters["nz"].get());
    }
    else if (input_parameters.count("bx") != 0)
    {
        INPUT.bx = *static_cast<int*>(input_parameters["bx"].get());
    }
    else if (input_parameters.count("by") != 0)
    {
        INPUT.by = *static_cast<int*>(input_parameters["by"].get());
    }
    else if (input_parameters.count("bz") != 0)
    {
        INPUT.bz = *static_cast<int*>(input_parameters["bz"].get());
    }
    else if (input_parameters.count("diago_proc") != 0)
    {
        INPUT.diago_proc = *static_cast<int*>(input_parameters["diago_proc"].get());
    }
    else if (input_parameters.count("pw_diag_nmax") != 0)
    {
        INPUT.pw_diag_nmax = *static_cast<int*>(input_parameters["pw_diag_nmax"].get());
    }
    else if (input_parameters.count("diago_cg_prec") != 0)
    {
        INPUT.diago_cg_prec = *static_cast<int*>(input_parameters["diago_cg_prec"].get());
    }
    else if (input_parameters.count("pw_diag_ndim") != 0)
    {
        INPUT.pw_diag_ndim = *static_cast<int*>(input_parameters["pw_diag_ndim"].get());
    }
    else if (input_parameters.count("pw_diag_thr") != 0)
    {
        INPUT.pw_diag_thr = *static_cast<double*>(input_parameters["pw_diag_thr"].get());
    }
    else if (input_parameters.count("nb2d") != 0)
    {
        INPUT.nb2d = *static_cast<int*>(input_parameters["nb2d"].get());
    }
    else if (input_parameters.count("nurse") != 0)
    {
        INPUT.nurse = *static_cast<int*>(input_parameters["nurse"].get());
    }
    else if (input_parameters.count("nbspline") != 0)
    {
        INPUT.nbspline = *static_cast<int*>(input_parameters["nbspline"].get());
    }
    else if (input_parameters.count("colour") != 0)
    {
        INPUT.colour = *static_cast<bool*>(input_parameters["colour"].get());
    }
    else if (input_parameters.count("t_in_h") != 0)
    {
        INPUT.t_in_h = *static_cast<bool*>(input_parameters["t_in_h"].get());
    }
    else if (input_parameters.count("vl_in_h") != 0)
    {
        INPUT.vl_in_h = *static_cast<bool*>(input_parameters["vl_in_h"].get());
    }
    else if (input_parameters.count("vnl_in_h") != 0)
    {
        INPUT.vnl_in_h = *static_cast<bool*>(input_parameters["vnl_in_h"].get());
    }
    else if (input_parameters.count("vh_in_h") != 0)
    {
        INPUT.vh_in_h = *static_cast<bool*>(input_parameters["vh_in_h"].get());
    }
    else if (input_parameters.count("vion_in_h") != 0)
    {
        INPUT.vion_in_h = *static_cast<bool*>(input_parameters["vion_in_h"].get());
    }
    else if (input_parameters.count("test_force") != 0)
    {
        INPUT.test_force = *static_cast<bool*>(input_parameters["test_force"].get());
    }
    else if (input_parameters.count("test_stress") != 0)
    {
        INPUT.test_stress = *static_cast<bool*>(input_parameters["test_stress"].get());
    }
    else if (input_parameters.count("scf_thr") != 0)
    {
        INPUT.scf_thr = *static_cast<double*>(input_parameters["scf_thr"].get());
    }
    else if (input_parameters.count("scf_thr_type") != 0)
    {
        INPUT.scf_thr_type = *static_cast<int*>(input_parameters["scf_thr_type"].get());
    }
    else if (input_parameters.count("scf_nmax") != 0)
    {
        INPUT.scf_nmax = *static_cast<int*>(input_parameters["scf_nmax"].get());
    }
    else if (input_parameters.count("relax_nmax") != 0)
    {
        INPUT.relax_nmax = *static_cast<int*>(input_parameters["relax_nmax"].get());
    }
    else if (input_parameters.count("out_stru") != 0)
    {
        INPUT.out_stru = *static_cast<bool*>(input_parameters["out_stru"].get());
    }
    else if (input_parameters.count("out_level") != 0)
    {
        INPUT.out_level = static_cast<SimpleString*>(input_parameters["out_level"].get())->c_str();
    }
    else if (input_parameters.count("out_md_control") != 0)
    {
        INPUT.out_md_control = *static_cast<bool*>(input_parameters["out_md_control"].get());
    }
    else if (input_parameters.count("occupations") != 0)
    {
        INPUT.occupations = static_cast<SimpleString*>(input_parameters["occupations"].get())->c_str();
    }
    else if (input_parameters.count("smearing_method") != 0)
    {
        INPUT.smearing_method = static_cast<SimpleString*>(input_parameters["smearing_method"].get())->c_str();
    }
    else if (input_parameters.count("smearing_sigma") != 0)
    {
        INPUT.smearing_sigma = *static_cast<double*>(input_parameters["smearing_sigma"].get());
    }
    else if (input_parameters.count("mixing_mode") != 0)
    {
        INPUT.mixing_mode = static_cast<SimpleString*>(input_parameters["mixing_mode"].get())->c_str();
    }
    else if (input_parameters.count("mixing_beta") != 0)
    {
        INPUT.mixing_beta = *static_cast<double*>(input_parameters["mixing_beta"].get());
    }
    else if (input_parameters.count("mixing_ndim") != 0)
    {
        INPUT.mixing_ndim = *static_cast<int*>(input_parameters["mixing_ndim"].get());
    }
    else if (input_parameters.count("mixing_gg0") != 0)
    {
        INPUT.mixing_gg0 = *static_cast<double*>(input_parameters["mixing_gg0"].get());
    }
    else if (input_parameters.count("mixing_beta_mag") != 0)
    {
        INPUT.mixing_beta_mag = *static_cast<double*>(input_parameters["mixing_beta_mag"].get());
    }
    else if (input_parameters.count("mixing_gg0_mag") != 0)
    {
        INPUT.mixing_gg0_mag = *static_cast<double*>(input_parameters["mixing_gg0_mag"].get());
    }
    else if (input_parameters.count("mixing_gg0_min") != 0)
    {
        INPUT.mixing_gg0_min = *static_cast<double*>(input_parameters["mixing_gg0_min"].get());
    }
    else if (input_parameters.count("mixing_angle") != 0)
    {
        INPUT.mixing_angle = *static_cast<double*>(input_parameters["mixing_angle"].get());
    }
    else if (input_parameters.count("mixing_tau") != 0)
    {
        INPUT.mixing_tau = *static_cast<bool*>(input_parameters["mixing_tau"].get());
    }
    else if (input_parameters.count("mixing_dftu") != 0)
    {
        INPUT.mixing_dftu = *static_cast<bool*>(input_parameters["mixing_dftu"].get());
    }
    else if (input_parameters.count("init_wfc") != 0)
    {
        INPUT.init_wfc = static_cast<SimpleString*>(input_parameters["init_wfc"].get())->c_str();
    }
    else if (input_parameters.count("init_chg") != 0)
    {
        INPUT.init_chg = static_cast<SimpleString*>(input_parameters["init_chg"].get())->c_str();
    }
    else if (input_parameters.count("chg_extrap") != 0)
    {
        INPUT.chg_extrap = static_cast<SimpleString*>(input_parameters["chg_extrap"].get())->c_str();
    }
    else if (input_parameters.count("mem_saver") != 0)
    {
        INPUT.mem_saver = *static_cast<int*>(input_parameters["mem_saver"].get());
    }
    else if (input_parameters.count("printe") != 0)
    {
        INPUT.printe = *static_cast<int*>(input_parameters["printe"].get());
    }
    else if (input_parameters.count("out_freq_elec") != 0)
    {
        INPUT.out_freq_elec = *static_cast<int*>(input_parameters["out_freq_elec"].get());
    }
    else if (input_parameters.count("out_freq_ion") != 0)
    {
        INPUT.out_freq_ion = *static_cast<int*>(input_parameters["out_freq_ion"].get());
    }
    else if (input_parameters.count("out_chg") != 0)
    {
        INPUT.out_chg = *static_cast<bool*>(input_parameters["out_chg"].get());
    }
    else if (input_parameters.count("out_dm") != 0)
    {
        INPUT.out_dm = *static_cast<bool*>(input_parameters["out_dm"].get());
    }
    else if (input_parameters.count("out_dm1") != 0)
    {
        INPUT.out_dm1 = *static_cast<bool*>(input_parameters["out_dm1"].get());
    }
    else if (input_parameters.count("out_pot") != 0)
    {
        INPUT.out_pot = *static_cast<int*>(input_parameters["out_pot"].get());
    }
    else if (input_parameters.count("out_wfc_pw") != 0)
    {
        INPUT.out_wfc_pw = *static_cast<int*>(input_parameters["out_wfc_pw"].get());
    }
    else if (input_parameters.count("out_wfc_r") != 0)
    {
        INPUT.out_wfc_r = *static_cast<bool*>(input_parameters["out_wfc_r"].get());
    }
    else if (input_parameters.count("out_dos") != 0)
    {
        INPUT.out_dos = *static_cast<bool*>(input_parameters["out_dos"].get());
    }
    else if (input_parameters.count("out_band") != 0)
    {
        INPUT.out_band = *static_cast<bool*>(input_parameters["out_band"].get());
    }
    else if (input_parameters.count("out_proj_band") != 0)
    {
        INPUT.out_proj_band = *static_cast<bool*>(input_parameters["out_proj_band"].get());
    }
    else if (input_parameters.count("out_mat_hs") != 0)
    {
        INPUT.out_mat_hs = *static_cast<bool*>(input_parameters["out_mat_hs"].get());
    }
    else if (input_parameters.count("out_mat_xc") != 0)
    {
        INPUT.out_mat_xc = *static_cast<bool*>(input_parameters["out_mat_xc"].get());
    }
    else if (input_parameters.count("cal_syns") != 0)
    {
        INPUT.cal_syns = *static_cast<bool*>(input_parameters["cal_syns"].get());
    }
    else if (input_parameters.count("dmax") != 0)
    {
        INPUT.dmax = *static_cast<double*>(input_parameters["dmax"].get());
    }
    else if (input_parameters.count("out_mat_hs2") != 0)
    {
        INPUT.out_mat_hs2 = *static_cast<bool*>(input_parameters["out_mat_hs2"].get());
    }
    else if (input_parameters.count("out_mat_dh") != 0)
    {
        INPUT.out_mat_dh = *static_cast<bool*>(input_parameters["out_mat_dh"].get());
    }
    else if (input_parameters.count("out_interval") != 0)
    {
        INPUT.out_interval = *static_cast<int*>(input_parameters["out_interval"].get());
    }
    else if (input_parameters.count("out_app_flag") != 0)
    {
        INPUT.out_app_flag = *static_cast<bool*>(input_parameters["out_app_flag"].get());
    }
    else if (input_parameters.count("out_mat_t") != 0)
    {
        INPUT.out_mat_t = *static_cast<bool*>(input_parameters["out_mat_t"].get());
    }
    else if (input_parameters.count("out_mat_r") != 0)
    {
        INPUT.out_mat_r = *static_cast<bool*>(input_parameters["out_mat_r"].get());
    }
    else if (input_parameters.count("out_wfc_lcao") != 0)
    {
        INPUT.out_wfc_lcao = *static_cast<bool*>(input_parameters["out_wfc_lcao"].get());
    }
    else if (input_parameters.count("out_alllog") != 0)
    {
        INPUT.out_alllog = *static_cast<bool*>(input_parameters["out_alllog"].get());
    }
    else if (input_parameters.count("out_element_info") != 0)
    {
        INPUT.out_element_info = *static_cast<bool*>(input_parameters["out_element_info"].get());
    }
    else if (input_parameters.count("out_bandgap") != 0)
    {
        INPUT.out_bandgap = *static_cast<bool*>(input_parameters["out_bandgap"].get());
    }
    else if (input_parameters.count("dos_emin_ev") != 0)
    {
        INPUT.dos_emin_ev = *static_cast<double*>(input_parameters["dos_emin_ev"].get());
    }
    else if (input_parameters.count("dos_emax_ev") != 0)
    {
        INPUT.dos_emax_ev = *static_cast<double*>(input_parameters["dos_emax_ev"].get());
    }
    if (input_parameters.count("dos_edelta_ev") != 0)
    {
        INPUT.dos_edelta_ev = *static_cast<double*>(input_parameters["dos_edelta_ev"].get());
    }
    else if (input_parameters.count("dos_scale") != 0)
    {
        INPUT.dos_scale = *static_cast<double*>(input_parameters["dos_scale"].get());
    }
    else if (input_parameters.count("dos_nche") != 0)
    {
        INPUT.dos_nche = *static_cast<int*>(input_parameters["dos_nche"].get());
    }
    else if (input_parameters.count("dos_setemin") != 0)
    {
        INPUT.dos_setemin = *static_cast<bool*>(input_parameters["dos_setemin"].get());
    }
    else if (input_parameters.count("dos_setemax") != 0)
    {
        INPUT.dos_setemax = *static_cast<bool*>(input_parameters["dos_setemax"].get());
    }
    else if (input_parameters.count("dos_sigma") != 0)
    {
        INPUT.dos_sigma = *static_cast<double*>(input_parameters["dos_sigma"].get());
    }
    else if (input_parameters.count("lcao_ecut") != 0)
    {
        INPUT.lcao_ecut = *static_cast<double*>(input_parameters["lcao_ecut"].get());
    }
    else if (input_parameters.count("lcao_dk") != 0)
    {
        INPUT.lcao_dk = *static_cast<double*>(input_parameters["lcao_dk"].get());
    }
    else if (input_parameters.count("lcao_dr") != 0)
    {
        INPUT.lcao_dr = *static_cast<double*>(input_parameters["lcao_dr"].get());
    }
    else if (input_parameters.count("lcao_rmax") != 0)
    {
        INPUT.lcao_rmax = *static_cast<double*>(input_parameters["lcao_rmax"].get());
    }
    else if (input_parameters.count("search_radius") != 0)
    {
        INPUT.search_radius = *static_cast<double*>(input_parameters["search_radius"].get());
    }
    else if (input_parameters.count("search_pbc") != 0)
    {
        INPUT.search_pbc = *static_cast<bool*>(input_parameters["search_pbc"].get());
    }
    else if (input_parameters.count("mdp") != 0)
    {
        // INPUT.mdp = static_cast<MD_para>(input_parameters["mdp"].get());
    }
    else if (input_parameters.count("efield_flag") != 0)
    {
        INPUT.efield_flag = *static_cast<bool*>(input_parameters["efield_flag"].get());
    }
    else if (input_parameters.count("dip_cor_flag") != 0)
    {
        INPUT.dip_cor_flag = *static_cast<bool*>(input_parameters["dip_cor_flag"].get());
    }
    else if (input_parameters.count("efield_dir") != 0)
    {
        INPUT.efield_dir = *static_cast<int*>(input_parameters["efield_dir"].get());
    }
    else if (input_parameters.count("efield_pos_max") != 0)
    {
        INPUT.efield_pos_max = *static_cast<double*>(input_parameters["efield_pos_max"].get());
    }
    else if (input_parameters.count("efield_pos_dec") != 0)
    {
        INPUT.efield_pos_dec = *static_cast<double*>(input_parameters["efield_pos_dec"].get());
    }
    else if (input_parameters.count("efield_amp") != 0)
    {
        INPUT.efield_amp = *static_cast<double*>(input_parameters["efield_amp"].get());
    }
    else if (input_parameters.count("gate_flag") != 0)
    {
        INPUT.gate_flag = *static_cast<bool*>(input_parameters["gate_flag"].get());
    }
    else if (input_parameters.count("zgate") != 0)
    {
        INPUT.zgate = *static_cast<double*>(input_parameters["zgate"].get());
    }
    else if (input_parameters.count("relax") != 0)
    {
        INPUT.relax = *static_cast<bool*>(input_parameters["relax"].get());
    }
    else if (input_parameters.count("block") != 0)
    {
        INPUT.block = *static_cast<bool*>(input_parameters["block"].get());
    }
    else if (input_parameters.count("block_down") != 0)
    {
        INPUT.block_down = *static_cast<double*>(input_parameters["block_down"].get());
    }
    else if (input_parameters.count("block_up") != 0)
    {
        INPUT.block_up = *static_cast<double*>(input_parameters["block_up"].get());
    }
    else if (input_parameters.count("block_height") != 0)
    {
        INPUT.block_height = *static_cast<double*>(input_parameters["block_height"].get());
    }
    else if (input_parameters.count("vdw_method") != 0)
    {
        INPUT.vdw_method = static_cast<SimpleString*>(input_parameters["vdw_method"].get())->c_str();
    }
    else if (input_parameters.count("vdw_s6") != 0)
    {
        INPUT.vdw_s6 = static_cast<SimpleString*>(input_parameters["vdw_s6"].get())->c_str();
    }
    else if (input_parameters.count("vdw_s8") != 0)
    {
        INPUT.vdw_s8 = static_cast<SimpleString*>(input_parameters["vdw_s8"].get())->c_str();
    }
    else if (input_parameters.count("vdw_a1") != 0)
    {
        INPUT.vdw_a1 = static_cast<SimpleString*>(input_parameters["vdw_a1"].get())->c_str();
    }
    else if (input_parameters.count("vdw_a2") != 0)
    {
        INPUT.vdw_a2 = static_cast<SimpleString*>(input_parameters["vdw_a2"].get())->c_str();
    }
    else if (input_parameters.count("vdw_d") != 0)
    {
        INPUT.vdw_d = *static_cast<double*>(input_parameters["vdw_d"].get());
    }
    else if (input_parameters.count("vdw_abc") != 0)
    {
        INPUT.vdw_abc = *static_cast<bool*>(input_parameters["vdw_abc"].get());
    }
    else if (input_parameters.count("vdw_cutoff_radius") != 0)
    {
        INPUT.vdw_cutoff_radius = static_cast<SimpleString*>(input_parameters["vdw_cutoff_radius"].get())->c_str();
    }
    else if (input_parameters.count("vdw_radius_unit") != 0)
    {
        INPUT.vdw_radius_unit = static_cast<SimpleString*>(input_parameters["vdw_radius_unit"].get())->c_str();
    }
    else if (input_parameters.count("vdw_cn_thr") != 0)
    {
        INPUT.vdw_cn_thr = *static_cast<double*>(input_parameters["vdw_cn_thr"].get());
    }
    else if (input_parameters.count("vdw_cn_thr_unit") != 0)
    {
        INPUT.vdw_cn_thr_unit = static_cast<SimpleString*>(input_parameters["vdw_cn_thr_unit"].get())->c_str();
    }
    else if (input_parameters.count("vdw_C6_file") != 0)
    {
        INPUT.vdw_C6_file = static_cast<SimpleString*>(input_parameters["vdw_C6_file"].get())->c_str();
    }
    else if (input_parameters.count("vdw_C6_unit") != 0)
    {
        INPUT.vdw_C6_unit = static_cast<SimpleString*>(input_parameters["vdw_C6_unit"].get())->c_str();
    }
    else if (input_parameters.count("vdw_R0_file") != 0)
    {
        INPUT.vdw_R0_file = static_cast<SimpleString*>(input_parameters["vdw_R0_file"].get())->c_str();
    }
    else if (input_parameters.count("vdw_R0_unit") != 0)
    {
        INPUT.vdw_R0_unit = static_cast<SimpleString*>(input_parameters["vdw_R0_unit"].get())->c_str();
    }
    else if (input_parameters.count("vdw_cutoff_type") != 0)
    {
        INPUT.vdw_cutoff_type = static_cast<SimpleString*>(input_parameters["vdw_cutoff_type"].get())->c_str();
    }
    else if (input_parameters.count("vdw_cutoff_period") != 0)
    {
        // INPUT.vdw_cutoff_period = static_cast<ModuleBase::Vector3<int>>(input_parameters["vdw_cutoff_period"].get());
    }
    else if (input_parameters.count("ocp") != 0)
    {
        INPUT.ocp = *static_cast<bool*>(input_parameters["ocp"].get());
    }
    else if (input_parameters.count("ocp_set") != 0)
    {
        INPUT.ocp_set = static_cast<SimpleString*>(input_parameters["ocp_set"].get())->c_str();
    }
    else if (input_parameters.count("out_mul") != 0)
    {
        INPUT.out_mul = *static_cast<bool*>(input_parameters["out_mul"].get());
    }
    else if (input_parameters.count("noncolin") != 0)
    {
        INPUT.noncolin = *static_cast<bool*>(input_parameters["noncolin"].get());
    }
    else if (input_parameters.count("lspinorb") != 0)
    {
        INPUT.lspinorb = *static_cast<bool*>(input_parameters["lspinorb"].get());
    }
    else if (input_parameters.count("soc_lambda") != 0)
    {
        INPUT.soc_lambda = *static_cast<double*>(input_parameters["soc_lambda"].get());
    }
    else if (input_parameters.count("exx_hybrid_alpha") != 0)
    {
        INPUT.exx_hybrid_alpha = static_cast<SimpleString*>(input_parameters["exx_hybrid_alpha"].get())->c_str();
    }
    else if (input_parameters.count("exx_hse_omega") != 0)
    {
        INPUT.exx_hse_omega = *static_cast<double*>(input_parameters["exx_hse_omega"].get());
    }
    else if (input_parameters.count("exx_separate_loop") != 0)
    {
        INPUT.exx_separate_loop = *static_cast<bool*>(input_parameters["exx_separate_loop"].get());
    }
    else if (input_parameters.count("exx_hybrid_step") != 0)
    {
        INPUT.exx_hybrid_step = *static_cast<int*>(input_parameters["exx_hybrid_step"].get());
    }
    else if (input_parameters.count("exx_mixing_beta") != 0)
    {
        INPUT.exx_mixing_beta = *static_cast<double*>(input_parameters["exx_mixing_beta"].get());
    }
    else if (input_parameters.count("exx_lambda") != 0)
    {
        INPUT.exx_lambda = *static_cast<double*>(input_parameters["exx_lambda"].get());
    }
    else if (input_parameters.count("exx_real_number") != 0)
    {
        INPUT.exx_real_number = static_cast<SimpleString*>(input_parameters["exx_real_number"].get())->c_str();
    }
    else if (input_parameters.count("exx_pca_threshold") != 0)
    {
        INPUT.exx_pca_threshold = *static_cast<double*>(input_parameters["exx_pca_threshold"].get());
    }
    else if (input_parameters.count("exx_c_threshold") != 0)
    {
        INPUT.exx_c_threshold = *static_cast<double*>(input_parameters["exx_c_threshold"].get());
    }
    else if (input_parameters.count("exx_v_threshold") != 0)
    {
        INPUT.exx_v_threshold = *static_cast<double*>(input_parameters["exx_v_threshold"].get());
    }
    else if (input_parameters.count("exx_dm_threshold") != 0)
    {
        INPUT.exx_dm_threshold = *static_cast<double*>(input_parameters["exx_dm_threshold"].get());
    }
    else if (input_parameters.count("exx_schwarz_threshold") != 0)
    {
        INPUT.exx_schwarz_threshold = *static_cast<double*>(input_parameters["exx_schwarz_threshold"].get());
    }
    else if (input_parameters.count("exx_cauchy_threshold") != 0)
    {
        INPUT.exx_cauchy_threshold = *static_cast<double*>(input_parameters["exx_cauchy_threshold"].get());
    }
    else if (input_parameters.count("exx_c_grad_threshold") != 0)
    {
        INPUT.exx_c_grad_threshold = *static_cast<double*>(input_parameters["exx_c_grad_threshold"].get());
    }
    else if (input_parameters.count("exx_v_grad_threshold") != 0)
    {
        INPUT.exx_v_grad_threshold = *static_cast<double*>(input_parameters["exx_v_grad_threshold"].get());
    }
    else if (input_parameters.count("exx_cauchy_force_threshold") != 0)
    {
        INPUT.exx_cauchy_force_threshold = *static_cast<double*>(input_parameters["exx_cauchy_force_threshold"].get());
    }
    else if (input_parameters.count("exx_cauchy_stress_threshold") != 0)
    {
        INPUT.exx_cauchy_stress_threshold
            = *static_cast<double*>(input_parameters["exx_cauchy_stress_threshold"].get());
    }
    else if (input_parameters.count("exx_ccp_threshold") != 0)
    {
        INPUT.exx_ccp_threshold = *static_cast<double*>(input_parameters["exx_ccp_threshold"].get());
    }
    else if (input_parameters.count("exx_ccp_rmesh_times") != 0)
    {
        INPUT.exx_ccp_rmesh_times = static_cast<SimpleString*>(input_parameters["exx_ccp_rmesh_times"].get())->c_str();
    }
    else if (input_parameters.count("exx_distribute_type") != 0)
    {
        INPUT.exx_distribute_type = static_cast<SimpleString*>(input_parameters["exx_distribute_type"].get())->c_str();
    }
    else if (input_parameters.count("exx_opt_orb_lmax") != 0)
    {
        INPUT.exx_opt_orb_lmax = *static_cast<int*>(input_parameters["exx_opt_orb_lmax"].get());
    }
    else if (input_parameters.count("exx_opt_orb_ecut") != 0)
    {
        INPUT.exx_opt_orb_ecut = *static_cast<double*>(input_parameters["exx_opt_orb_ecut"].get());
    }
    else if (input_parameters.count("exx_opt_orb_tolerence") != 0)
    {
        INPUT.exx_opt_orb_tolerence = *static_cast<double*>(input_parameters["exx_opt_orb_tolerence"].get());
    }
    else if (input_parameters.count("td_force_dt") != 0)
    {
        INPUT.td_force_dt = *static_cast<double*>(input_parameters["td_force_dt"].get());
    }
    else if (input_parameters.count("td_vext") != 0)
    {
        INPUT.td_vext = *static_cast<bool*>(input_parameters["td_vext"].get());
    }
    else if (input_parameters.count("td_vext_dire") != 0)
    {
        INPUT.td_vext_dire = static_cast<SimpleString*>(input_parameters["td_vext_dire"].get())->c_str();
    }
    else if (input_parameters.count("out_dipole") != 0)
    {
        INPUT.out_dipole = *static_cast<bool*>(input_parameters["out_dipole"].get());
    }
    else if (input_parameters.count("out_efield") != 0)
    {
        INPUT.out_efield = *static_cast<bool*>(input_parameters["out_efield"].get());
    }
    else if (input_parameters.count("td_print_eij") != 0)
    {
        INPUT.td_print_eij = *static_cast<double*>(input_parameters["td_print_eij"].get());
    }
    else if (input_parameters.count("td_edm") != 0)
    {
        INPUT.td_edm = *static_cast<int*>(input_parameters["td_edm"].get());
    }
    else if (input_parameters.count("propagator") != 0)
    {
        INPUT.propagator = *static_cast<int*>(input_parameters["propagator"].get());
    }
    else if (input_parameters.count("td_stype") != 0)
    {
        INPUT.td_stype = *static_cast<int*>(input_parameters["td_stype"].get());
    }
    else if (input_parameters.count("td_ttype") != 0)
    {
        INPUT.td_ttype = static_cast<SimpleString*>(input_parameters["td_ttype"].get())->c_str();
    }
    else if (input_parameters.count("td_tstart") != 0)
    {
        INPUT.td_tstart = *static_cast<int*>(input_parameters["td_tstart"].get());
    }
    else if (input_parameters.count("td_tend") != 0)
    {
        INPUT.td_tend = *static_cast<int*>(input_parameters["td_tend"].get());
    }
    else if (input_parameters.count("td_lcut1") != 0)
    {
        INPUT.td_lcut1 = *static_cast<double*>(input_parameters["td_lcut1"].get());
    }
    else if (input_parameters.count("td_lcut2") != 0)
    {
        INPUT.td_lcut2 = *static_cast<double*>(input_parameters["td_lcut2"].get());
    }
    else if (input_parameters.count("td_gauss_freq") != 0)
    {
        INPUT.td_gauss_freq = static_cast<SimpleString*>(input_parameters["td_gauss_freq"].get())->c_str();
    }
    else if (input_parameters.count("td_gauss_phase") != 0)
    {
        INPUT.td_gauss_phase = static_cast<SimpleString*>(input_parameters["td_gauss_phase"].get())->c_str();
    }
    else if (input_parameters.count("td_gauss_sigma") != 0)
    {
        INPUT.td_gauss_sigma = static_cast<SimpleString*>(input_parameters["td_gauss_sigma"].get())->c_str();
    }
    else if (input_parameters.count("td_gauss_t0") != 0)
    {
        INPUT.td_gauss_t0 = static_cast<SimpleString*>(input_parameters["td_gauss_t0"].get())->c_str();
    }
    else if (input_parameters.count("td_gauss_amp") != 0)
    {
        INPUT.td_gauss_amp = static_cast<SimpleString*>(input_parameters["td_gauss_amp"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_freq") != 0)
    {
        INPUT.td_trape_freq = static_cast<SimpleString*>(input_parameters["td_trape_freq"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_phase") != 0)
    {
        INPUT.td_trape_phase = static_cast<SimpleString*>(input_parameters["td_trape_phase"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_t1") != 0)
    {
        INPUT.td_trape_t1 = static_cast<SimpleString*>(input_parameters["td_trape_t1"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_t2") != 0)
    {
        INPUT.td_trape_t2 = static_cast<SimpleString*>(input_parameters["td_trape_t2"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_t3") != 0)
    {
        INPUT.td_trape_t3 = static_cast<SimpleString*>(input_parameters["td_trape_t3"].get())->c_str();
    }
    else if (input_parameters.count("td_trape_amp") != 0)
    {
        INPUT.td_trape_amp = static_cast<SimpleString*>(input_parameters["td_trape_amp"].get())->c_str();
    }
    else if (input_parameters.count("td_trigo_freq1") != 0)
    {
        INPUT.td_trigo_freq1 = static_cast<SimpleString*>(input_parameters["td_trigo_freq1"].get())->c_str();
    }
    else if (input_parameters.count("td_trigo_freq2") != 0)
    {
        INPUT.td_trigo_freq2 = static_cast<SimpleString*>(input_parameters["td_trigo_freq2"].get())->c_str();
    }
    else if (input_parameters.count("td_trigo_phase1") != 0)
    {
        INPUT.td_trigo_phase1 = static_cast<SimpleString*>(input_parameters["td_trigo_phase1"].get())->c_str();
    }
    else if (input_parameters.count("td_trigo_phase2") != 0)
    {
        INPUT.td_trigo_phase2 = static_cast<SimpleString*>(input_parameters["td_trigo_phase2"].get())->c_str();
    }
    else if (input_parameters.count("td_trigo_amp") != 0)
    {
        INPUT.td_trigo_amp = static_cast<SimpleString*>(input_parameters["td_trigo_amp"].get())->c_str();
    }
    else if (input_parameters.count("td_heavi_t0") != 0)
    {
        INPUT.td_heavi_t0 = static_cast<SimpleString*>(input_parameters["td_heavi_t0"].get())->c_str();
    }
    else if (input_parameters.count("td_heavi_amp") != 0)
    {
        INPUT.td_heavi_amp = static_cast<SimpleString*>(input_parameters["td_heavi_amp"].get())->c_str();
    }
    else if (input_parameters.count("restart_save") != 0)
    {
        INPUT.restart_save = *static_cast<bool*>(input_parameters["restart_save"].get());
    }
    else if (input_parameters.count("restart_load") != 0)
    {
        INPUT.restart_load = *static_cast<bool*>(input_parameters["restart_load"].get());
    }
    else if (input_parameters.count("input_error") != 0)
    {
        INPUT.input_error = *static_cast<bool*>(input_parameters["input_error"].get());
    }
    else if (input_parameters.count("cell_factor") != 0)
    {
        INPUT.cell_factor = *static_cast<double*>(input_parameters["cell_factor"].get());
    }
    else if (input_parameters.count("dft_plus_u") != 0)
    {
        INPUT.dft_plus_u = *static_cast<bool*>(input_parameters["dft_plus_u"].get());
    }
    else if (input_parameters.count("orbital_corr") != 0)
    {
        INPUT.dft_plus_u = 0;
        bool dmft_flag = false;
        SimpleVector<double> vec_D = *static_cast<SimpleVector<double>*>(input_parameters["orbital_corr"].get());
        for (int i = 0; i < INPUT.ntype; i++)
        {
            INPUT.orbital_corr[i] = vec_D[i];
            if ((INPUT.orbital_corr[i] != -1) && (INPUT.orbital_corr[i] != 0) && (INPUT.orbital_corr[i] != 1)
                && (INPUT.orbital_corr[i] != 2) && (INPUT.orbital_corr[i] != 3))
            {
                std::cout << " WRONG ARGUMENTS OF orbital_corr " << std::endl;
                exit(0);
            }
            if (INPUT.orbital_corr[i] != -1)
            {
                dmft_flag = true;
                INPUT.dft_plus_u = 1;
            }
        }
        if (!dmft_flag)
        {
            std::cout << "No atoms are correlated!!!" << std::endl;
            exit(0);
        }

        if (strcmp("lcao", INPUT.basis_type.c_str()) != 0)
        {
            std::cout << " WRONG ARGUMENTS OF basis_type, only lcao is support " << std::endl;
            exit(0);
        }
    }
    else if (input_parameters.count("hubbard_u") != 0)
    {
        SimpleVector<double> vec_D = *static_cast<SimpleVector<double>*>(input_parameters["hubbard_u"].get());
        for (int i = 0; i < INPUT.ntype; i++)
        {
            INPUT.hubbard_u[i] = vec_D[i] / ModuleBase::Ry_to_eV;
            if (INPUT.hubbard_u[i] < -1.0e-3)
            {
                std::cout << " WRONG ARGUMENTS OF hubbard_u " << std::endl;
                exit(0);
            }
        }
    }
    else if (input_parameters.count("omc") != 0)
    {
        INPUT.omc = *static_cast<int*>(input_parameters["omc"].get());
    }
    else if (input_parameters.count("yukawa_potential") != 0)
    {
        INPUT.yukawa_potential = *static_cast<bool*>(input_parameters["yukawa_potential"].get());
    }
    else if (input_parameters.count("yukawa_lambda") != 0)
    {
        INPUT.yukawa_lambda = *static_cast<double*>(input_parameters["yukawa_lambda"].get());
    }
    else if (input_parameters.count("dft_plus_dmft") != 0)
    {
        INPUT.dft_plus_dmft = *static_cast<bool*>(input_parameters["dft_plus_dmft"].get());
    }
    else if (input_parameters.count("rpa") != 0)
    {
        INPUT.rpa = *static_cast<bool*>(input_parameters["rpa"].get());
    }
    else if (input_parameters.count("coulomb_type") != 0)
    {
        INPUT.coulomb_type = static_cast<SimpleString*>(input_parameters["coulomb_type"].get())->c_str();
    }
    else if (input_parameters.count("deepks_out_labels") != 0)
    {
        INPUT.deepks_out_labels = *static_cast<bool*>(input_parameters["deepks_out_labels"].get());
    }
    else if (input_parameters.count("deepks_scf") != 0)
    {
        INPUT.deepks_scf = *static_cast<bool*>(input_parameters["deepks_scf"].get());
    }
    else if (input_parameters.count("deepks_bandgap") != 0)
    {
        INPUT.deepks_bandgap = *static_cast<bool*>(input_parameters["deepks_bandgap"].get());
    }
    else if (input_parameters.count("deepks_out_unittest") != 0)
    {
        INPUT.deepks_out_unittest = *static_cast<bool*>(input_parameters["deepks_out_unittest"].get());
    }
    else if (input_parameters.count("deepks_model") != 0)
    {
        INPUT.deepks_model = static_cast<SimpleString*>(input_parameters["deepks_model"].get())->c_str();
    }
    else if (input_parameters.count("imp_sol") != 0)
    {
        INPUT.imp_sol = *static_cast<bool*>(input_parameters["imp_sol"].get());
    }
    else if (input_parameters.count("eb_k") != 0)
    {
        INPUT.eb_k = *static_cast<double*>(input_parameters["eb_k"].get());
    }
    else if (input_parameters.count("tau") != 0)
    {
        INPUT.tau = *static_cast<double*>(input_parameters["tau"].get());
    }
    else if (input_parameters.count("sigma_k") != 0)
    {
        INPUT.sigma_k = *static_cast<double*>(input_parameters["sigma_k"].get());
    }
    else if (input_parameters.count("nc_k") != 0)
    {
        INPUT.nc_k = *static_cast<double*>(input_parameters["nc_k"].get());
    }
    else if (input_parameters.count("of_kinetic") != 0)
    {
        INPUT.of_kinetic = static_cast<SimpleString*>(input_parameters["of_kinetic"].get())->c_str();
    }
    else if (input_parameters.count("of_method") != 0)
    {
        INPUT.of_method = static_cast<SimpleString*>(input_parameters["of_method"].get())->c_str();
    }
    else if (input_parameters.count("of_conv") != 0)
    {
        INPUT.of_conv = static_cast<SimpleString*>(input_parameters["of_conv"].get())->c_str();
    }
    else if (input_parameters.count("of_tole") != 0)
    {
        INPUT.of_tole = *static_cast<double*>(input_parameters["of_tole"].get());
    }
    else if (input_parameters.count("of_tolp") != 0)
    {
        INPUT.of_tolp = *static_cast<double*>(input_parameters["of_tolp"].get());
    }
    else if (input_parameters.count("of_tf_weight") != 0)
    {
        INPUT.of_tf_weight = *static_cast<double*>(input_parameters["of_tf_weight"].get());
    }
    else if (input_parameters.count("of_vw_weight") != 0)
    {
        INPUT.of_vw_weight = *static_cast<double*>(input_parameters["of_vw_weight"].get());
    }
    else if (input_parameters.count("of_wt_alpha") != 0)
    {
        INPUT.of_wt_alpha = *static_cast<double*>(input_parameters["of_wt_alpha"].get());
    }
    else if (input_parameters.count("of_wt_beta") != 0)
    {
        INPUT.of_wt_beta = *static_cast<double*>(input_parameters["of_wt_beta"].get());
    }
    else if (input_parameters.count("of_wt_rho0") != 0)
    {
        INPUT.of_wt_rho0 = *static_cast<double*>(input_parameters["of_wt_rho0"].get());
    }
    else if (input_parameters.count("of_hold_rho0") != 0)
    {
        INPUT.of_hold_rho0 = *static_cast<bool*>(input_parameters["of_hold_rho0"].get());
    }
    else if (input_parameters.count("of_lkt_a") != 0)
    {
        INPUT.of_lkt_a = *static_cast<double*>(input_parameters["of_lkt_a"].get());
    }
    else if (input_parameters.count("of_full_pw") != 0)
    {
        INPUT.of_full_pw = *static_cast<bool*>(input_parameters["of_full_pw"].get());
    }
    else if (input_parameters.count("of_full_pw_dim") != 0)
    {
        INPUT.of_full_pw_dim = *static_cast<int*>(input_parameters["of_full_pw_dim"].get());
    }
    else if (input_parameters.count("of_read_kernel") != 0)
    {
        INPUT.of_read_kernel = *static_cast<bool*>(input_parameters["of_read_kernel"].get());
    }
    else if (input_parameters.count("of_kernel_file") != 0)
    {
        INPUT.of_kernel_file = static_cast<SimpleString*>(input_parameters["of_kernel_file"].get())->c_str();
    }
    else if (input_parameters.count("bessel_nao_smooth") != 0)
    {
        INPUT.bessel_nao_smooth = *static_cast<bool*>(input_parameters["bessel_nao_smooth"].get());
    }
    else if (input_parameters.count("bessel_nao_sigma") != 0)
    {
        INPUT.bessel_nao_sigma = *static_cast<double*>(input_parameters["bessel_nao_sigma"].get());
    }
    else if (input_parameters.count("bessel_nao_ecut") != 0)
    {
        INPUT.bessel_nao_ecut = static_cast<SimpleString*>(input_parameters["bessel_nao_ecut"].get())->c_str();
    }
    else if (input_parameters.count("bessel_nao_rcut") != 0)
    {
        INPUT.bessel_nao_rcut = *static_cast<double*>(input_parameters["bessel_nao_rcut"].get());
    }
    else if (input_parameters.count("bessel_nao_tolerence") != 0)
    {
        INPUT.bessel_nao_tolerence = *static_cast<double*>(input_parameters["bessel_nao_tolerence"].get());
    }
    else if (input_parameters.count("bessel_descriptor_lmax") != 0)
    {
        INPUT.bessel_descriptor_lmax = *static_cast<int*>(input_parameters["bessel_descriptor_lmax"].get());
    }
    else if (input_parameters.count("bessel_descriptor_smooth") != 0)
    {
        INPUT.bessel_descriptor_smooth = *static_cast<bool*>(input_parameters["bessel_descriptor_smooth"].get());
    }
    else if (input_parameters.count("bessel_descriptor_sigma") != 0)
    {
        INPUT.bessel_descriptor_sigma = *static_cast<double*>(input_parameters["bessel_descriptor_sigma"].get());
    }
    else if (input_parameters.count("bessel_descriptor_ecut") != 0)
    {
        INPUT.bessel_descriptor_ecut
            = static_cast<SimpleString*>(input_parameters["bessel_descriptor_ecut"].get())->c_str();
    }
    else if (input_parameters.count("bessel_descriptor_rcut") != 0)
    {
        INPUT.bessel_descriptor_rcut = *static_cast<double*>(input_parameters["bessel_descriptor_rcut"].get());
    }
    else if (input_parameters.count("bessel_descriptor_tolerence") != 0)
    {
        INPUT.bessel_descriptor_tolerence
            = *static_cast<double*>(input_parameters["bessel_descriptor_tolerence"].get());
    }
    else if (input_parameters.count("device") != 0)
    {
        INPUT.device = static_cast<SimpleString*>(input_parameters["device"].get())->c_str();
    }
    else if (input_parameters.count("precision") != 0)
    {
        INPUT.precision = static_cast<SimpleString*>(input_parameters["precision"].get())->c_str();
    }
    else if (input_parameters.count("test_skip_ewald") != 0)
    {
        INPUT.test_skip_ewald = *static_cast<bool*>(input_parameters["test_skip_ewald"].get());
    }
}

// namespace ModuleIO
} // namespace ModuleIO
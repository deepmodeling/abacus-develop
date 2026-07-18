#include "deepks_test.h"
#include "source_base/global_variable.h"
#include "source_estate/read_pseudo.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_io/module_parameter/parameter.h"
#include "../../LCAO_nonlocal_info.h"

#include <gtest/gtest.h>

namespace
{
Input_para& mutable_input_for_deepks_unit()
{
    return const_cast<Input_para&>(PARAM.inp);
}

System_para& mutable_system_for_deepks_unit()
{
    return const_cast<System_para&>(PARAM.globalv);
}
} // namespace

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}
namespace GlobalC
{
Exx_Info exx_info;
}

template <typename T>
void test_deepks<T>::preparation()
{
    this->count_ntype();
    this->set_parameters();
    if (testing::Test::HasFatalFailure())
    {
        return;
    }

    this->setup_cell();

    this->setup_kpt();

    this->set_ekcut();
    this->set_orbs();
    this->prep_neighbour();

    this->ParaO.set_serial(this->nlocal, this->nlocal);
    this->ParaO.nrow_bands = this->nlocal;
    this->ParaO.ncol_bands = this->nbands;
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06

    this->ParaO.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, this->nlocal);
}

template <typename T>
void test_deepks<T>::set_parameters()
{
    Input_para& input = mutable_input_for_deepks_unit();
    System_para& system = mutable_system_for_deepks_unit();

    input.basis_type = "lcao";
    input.kpoint_file = "KPT";
    input.pseudo_rcut = 15.0;
    input.cal_force = this->cal_force;
    input.gamma_only = this->gamma_only_local;
    input.nspin = this->nspin;
    input.orbital_dir = this->orbital_dir;
    input.out_element_info = this->out_element_info;
    system.global_out_dir = "./";
    GlobalV::ofs_warning.open("warning.log");
    GlobalV::ofs_running.open("running.log");
    system.deepks_setorb = this->deepks_setorb;

    std::ifstream ifs("INPUT");
    ASSERT_TRUE(ifs.is_open()) << "Cannot open DeePKS unit-test INPUT";
    char word[80];
    ASSERT_TRUE(ifs >> word);
    ASSERT_STREQ(word, "gamma_only_local");
    ASSERT_TRUE(ifs >> this->gamma_only_local);
    ifs.close();

    input.gamma_only = this->gamma_only_local;
    system.gamma_only_local = this->gamma_only_local;
    system.npol = this->npol;
    GlobalV::KPAR = 1;
    GlobalV::MY_POOL = 0;
    GlobalV::RANK_IN_POOL = 0;
    GlobalV::NPROC_IN_POOL = 1;

    ucell.latName = "user_defined_lattice";
    ucell.ntype = ntype;
    return;
}

template <typename T>
void test_deepks<T>::count_ntype()
{
    GlobalV::ofs_running << "count number of atom types" << std::endl;
    std::ifstream ifs("STRU", std::ios::in);

    if (!ifs)
    {
        GlobalV::ofs_running << "ERROR : file STRU does not exist" << std::endl;
        exit(1);
    }

    ModuleBase::GlobalFunc::SCAN_BEGIN(ifs, "ATOMIC_SPECIES");

    ntype = 0;

    std::string x;
    ifs.rdstate();
    while (ifs.good())
    {
        // read a line
        std::getline(ifs, x);

        // trim white space
        const char* typeOfWhitespaces = " \t\n\r\f\v";
        x.erase(x.find_last_not_of(typeOfWhitespaces) + 1);
        x.erase(0, x.find_first_not_of(typeOfWhitespaces));

        if (x == "LATTICE_CONSTANT" || x == "NUMERICAL_ORBITAL" || x == "LATTICE_VECTORS" || x == "ATOMIC_POSITIONS"
            || x == "NUMERICAL_DESCRIPTOR")
        {
            break;
        }

        std::string tmpid = x.substr(0, 1);
        if (!x.empty() && tmpid != "#")
        {
            ntype++;
        }
    }

    GlobalV::ofs_running << "ntype : " << ntype << std::endl;
    ifs.close();

    return;
}

template <typename T>
void test_deepks<T>::set_ekcut()
{
    GlobalV::ofs_running << "set lcao_ecut from LCAO files" << std::endl;
    // set as max of ekcut from every element

    lcao_ecut = 0.0;
    std::ifstream in_ao;

    for (int it = 0; it < ntype; it++)
    {
        double ek_current;

        in_ao.open(ucell.orbital_fn[it].c_str());
        if (!in_ao)
        {
            GlobalV::ofs_running << "error : cannot find LCAO file : " << ucell.orbital_fn[it] << std::endl;
        }

        std::string word;
        while (in_ao.good())
        {
            in_ao >> word;
            if (word == "Cutoff(Ry)")
            {
                break;
            }
        }
        in_ao >> ek_current;
        lcao_ecut = std::max(lcao_ecut, ek_current);

        in_ao.close();
    }

    ORB.ecutwfc = lcao_ecut;
    GlobalV::ofs_running << "lcao_ecut : " << lcao_ecut << std::endl;

    return;
}

template <typename T>
void test_deepks<T>::setup_cell()
{
    ucell.setup_cell("STRU", GlobalV::ofs_running, PARAM.inp.symmetry_prec, PARAM.inp.dfthalf_type, PARAM.inp.pseudo_dir, PARAM.inp.nspin,
        PARAM.inp.basis_type, PARAM.inp.orbital_dir, PARAM.inp.init_wfc,
        PARAM.inp.onsite_radius, PARAM.globalv.deepks_setorb, PARAM.inp.rpa,
        PARAM.inp.fixed_atoms, PARAM.inp.noncolin, PARAM.inp.calculation, PARAM.inp.esolver_type);
    const std::string pseudo_dir = PARAM.inp.pseudo_dir;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool out_element_info = PARAM.inp.out_element_info;
    const std::string dft_functional = PARAM.inp.dft_functional;
    const bool lspinorb = PARAM.inp.lspinorb;
    const double pseudo_rcut = PARAM.inp.pseudo_rcut;
    const double soc_lambda = PARAM.inp.soc_lambda;
    const int nspin = PARAM.inp.nspin;
    const int nlocal = PARAM.globalv.nlocal;
    const int npol = PARAM.globalv.npol;
    const std::string basis_type = PARAM.inp.basis_type;
    const std::string esolver_type = PARAM.inp.esolver_type;
    const std::string init_wfc = PARAM.inp.init_wfc;
    const int nbands = PARAM.inp.nbands;
    elecstate::read_pseudo(GlobalV::ofs_running, ucell, pseudo_dir, global_out_dir, out_element_info, dft_functional, lspinorb, pseudo_rcut, soc_lambda, nspin, nlocal, npol, basis_type, esolver_type, init_wfc, nbands);
    this->nlocal = PARAM.globalv.nlocal;
    this->nbands = PARAM.inp.nbands;
    this->npol = PARAM.globalv.npol;

    return;
}

template <typename T>
void test_deepks<T>::prep_neighbour()
{
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                   this->out_level,
                                                   ORB.get_rcutmax_Phi(),
                                                   ucell.infoNL->get_rcutmax_Beta(),
                                                   this->gamma_only_local);

    atom_arrange::search(this->search_pbc,
                         GlobalV::ofs_running,
                         Test_Deepks::GridD,
                         ucell,
                         search_radius,
                         this->test_atom_input);
}

template <typename T>
void test_deepks<T>::set_orbs()
{
    ORB.init(GlobalV::ofs_running,
             ucell.ntype,
             this->orbital_dir,
             ucell.orbital_fn.data(),
             ucell.descriptor_file,
             ucell.lmax,
             lcao_ecut,
             lcao_dk,
             lcao_dr,
             lcao_rmax,
             this->deepks_setorb,
             out_mat_r,
             this->out_element_info,
             this->cal_force,
             my_rank);

    auto* lcao_nl = new LCAONonlocalInfo();
    lcao_nl->setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, ORB);
    ucell.infoNL = lcao_nl;

    orb_.build(ntype, ucell.orbital_fn.data());

    std::string file_alpha = this->orbital_dir + ucell.descriptor_file;
    alpha_.build(1, &file_alpha);

    double rmax = std::max(orb_.rcut_max(), alpha_.rcut_max());
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / lcao_dr) + 1;

    orb_.set_uniform_grid(true, nr, cutoff, 'i', true);
    alpha_.set_uniform_grid(true, nr, cutoff, 'i', true);

    overlap_orb_alpha_.tabulate(orb_, alpha_, 'S', nr, cutoff);

    return;
}

template <typename T>
void test_deepks<T>::setup_kpt()
{
    ModuleSymmetry::Symmetry::symm_flag = -1;
    const bool use_ibz = false;
    const std::string global_out_dir = PARAM.globalv.global_out_dir;
    const bool gamma_only_local = PARAM.globalv.gamma_only_local;
    const double kspacing[3] = {PARAM.inp.kspacing[0], PARAM.inp.kspacing[1], PARAM.inp.kspacing[2]};
    const std::string kmesh_type = PARAM.inp.kmesh_type;
    const double koffset[3] = {PARAM.inp.koffset[0], PARAM.inp.koffset[1], PARAM.inp.koffset[2]};
    this->kv.set(ucell,
                 ucell.symm,
                 PARAM.inp.kpoint_file,
                 this->nspin,
                 ucell.G,
                 ucell.latvec,
                 GlobalV::ofs_running,
                 use_ibz,
                 global_out_dir,
                 gamma_only_local,
                 kspacing,
                 kmesh_type,
                 koffset);
}

template class test_deepks<double>;
template class test_deepks<std::complex<double>>;

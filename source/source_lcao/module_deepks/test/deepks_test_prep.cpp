#include "deepks_test.h"
#include "source_base/global_variable.h"
#include "source_estate/read_pseudo.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "../../LCAO_nonlocal_info.h"

#include <gtest/gtest.h>

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

    this->ParaO.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, this->nlocal);
}

template <typename T>
void test_deepks<T>::set_parameters()
{
    std::ifstream ifs("INPUT");
    ASSERT_TRUE(ifs.is_open()) << "Cannot open DeePKS unit-test INPUT";
    char word[80];
    ASSERT_TRUE(ifs >> word);
    ASSERT_STREQ(word, "gamma_only_local");
    ASSERT_TRUE(ifs >> this->gamma_only_local);
    ifs.close();

    GlobalV::ofs_warning.open("warning.log");
    GlobalV::ofs_running.open("running.log");
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
        std::getline(ifs, x);

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
    const std::string basis_type = "lcao";
    const std::string orbital_dir = this->orbital_dir;
    const std::string init_wfc = "atomic";
    const double onsite_radius = 0.0;
    const bool deepks_setorb = this->deepks_setorb;
    const bool rpa = false;
    const bool fixed_atoms = false;
    const bool noncolin = false;
    const std::string calculation = "scf";
    const std::string esolver_type = "cg";
    const double symmetry_prec = 1e-5;
    const int dfthalf_type = 0;
    const std::string pseudo_dir = "";
    const int nspin = this->nspin;

    ucell.setup_cell("STRU", GlobalV::ofs_running, symmetry_prec, dfthalf_type, pseudo_dir, nspin,
        basis_type, orbital_dir, init_wfc,
        onsite_radius, deepks_setorb, rpa,
        fixed_atoms, noncolin, calculation, esolver_type);

    const std::string global_out_dir = "./";
    const bool out_element_info = this->out_element_info;
    const std::string dft_functional = "default";
    const bool lspinorb = false;
    const double pseudo_rcut = 15.0;
    const double soc_lambda = 0.0;
    const int npol = this->npol;
    const int nbands = this->nbands;
    const bool two_fermi = false;
    const double nelec_delta = 0.0;
    const std::string smearing_method = "gaussian";
    const std::string ks_solver = "cg";
    const int bndpar = 1;

    auto atoms_info = elecstate::read_pseudo(GlobalV::ofs_running, ucell, pseudo_dir, global_out_dir, out_element_info, dft_functional, lspinorb, pseudo_rcut, soc_lambda, nspin, npol, basis_type, esolver_type, init_wfc, nbands, two_fermi, nelec_delta, smearing_method, ks_solver, bndpar);

    this->nlocal = atoms_info.nlocal;
    this->nbands = atoms_info.nbands;

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

    const std::string basis_type = "lcao";
    const bool out_element_info = this->out_element_info;
    const bool lspinorb = false;
    const int nspin = this->nspin;

    auto* lcao_nl = new LCAONonlocalInfo();
    lcao_nl->setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, ORB,
                           basis_type, out_element_info, lspinorb, nspin);
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
    const std::string global_out_dir = "./";
    const bool gamma_only_local = this->gamma_only_local;
    const double kspacing[3] = {0.0, 0.0, 0.0};
    const std::string kmesh_type = "gamma";
    const double koffset[3] = {0.0, 0.0, 0.0};
    const std::string kpoint_file = "KPT";

    this->kv.set(ucell,
                 ucell.symm,
                 kpoint_file,
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
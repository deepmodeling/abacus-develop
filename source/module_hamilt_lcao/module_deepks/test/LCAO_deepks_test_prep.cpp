#include "LCAO_deepks_test.h"
#include "module_base/global_variable.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#include "module_elecstate/read_pseudo.h"
#include "module_hamilt_general/module_xc/exx_info.h"

Magnetism::Magnetism() {
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism() { delete[] this->start_magnetization; }
namespace GlobalC
{
	Exx_Info exx_info;
}

void test_deepks::preparation()
{
    this->count_ntype();
    this->set_parameters();

    this->setup_cell();

    this->setup_kpt();

    this->set_ekcut();
    this->set_orbs();
    this->prep_neighbour();

    this->ParaO.set_serial(PARAM.globalv.nlocal, PARAM.globalv.nlocal);
    this->ParaO.nrow_bands = PARAM.globalv.nlocal;
    this->ParaO.ncol_bands = PARAM.inp.nbands;
    // Zhang Xiaoyang enable the serial version of LCAO and recovered this function usage. 2024-07-06

    this->ParaO.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, PARAM.globalv.nlocal);
}

void test_deepks::set_parameters()
{
    PARAM.input.basis_type = "lcao";
    // GlobalV::global_pseudo_type= "auto";
    PARAM.input.pseudo_rcut = 15.0;
    PARAM.sys.global_out_dir = "./";
    GlobalV::ofs_warning.open("warning.log");
    GlobalV::ofs_running.open("running.log");
    PARAM.sys.deepks_setorb = true;
    PARAM.input.cal_force = 1;

    std::ifstream ifs("INPUT");
    char word[80];
    ifs >> word;
    ifs >> PARAM.sys.gamma_only_local;
    ifs.close();

    ucell.latName = "none";
    ucell.ntype = ntype;
    return;
}

void test_deepks::count_ntype()
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
            || x == "NUMERICAL_DESCRIPTOR") {
            break;
}

        std::string tmpid = x.substr(0, 1);
        if (!x.empty() && tmpid != "#") {
            ntype++;
}
    }

    GlobalV::ofs_running << "ntype : " << ntype << std::endl;
    ifs.close();

    return;
}

void test_deepks::set_ekcut()
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
            if (word == "Cutoff(Ry)") {
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

void test_deepks::setup_cell()
{
    ucell.setup_cell("STRU", GlobalV::ofs_running);
    elecstate::read_pseudo(GlobalV::ofs_running, ucell);

    return;
}

void test_deepks::prep_neighbour()
{
    double search_radius = atom_arrange::set_sr_NL(GlobalV::ofs_running,
                                                   PARAM.input.out_level,
                                                   ORB.get_rcutmax_Phi(),
                                                   ucell.infoNL.get_rcutmax_Beta(),
                                                   PARAM.sys.gamma_only_local);

    atom_arrange::search(PARAM.inp.search_pbc,
                         GlobalV::ofs_running,
                         Test_Deepks::GridD,
                         ucell,
                         search_radius,
                         PARAM.inp.test_atom_input);
}

void test_deepks::set_orbs()
{
    ORB.init(GlobalV::ofs_running,
                        ucell.ntype,
                        PARAM.inp.orbital_dir,
                        ucell.orbital_fn,
                        ucell.descriptor_file,
                        ucell.lmax,
                        lcao_ecut,
                        lcao_dk,
                        lcao_dr,
                        lcao_rmax,
                        PARAM.sys.deepks_setorb,
                        out_mat_r,
                        PARAM.input.cal_force,
                        my_rank);

    ucell.infoNL.setupNonlocal(ucell.ntype, ucell.atoms, GlobalV::ofs_running, ORB);

    std::vector<std::string> file_orb(ntype);
    std::transform(ucell.orbital_fn, ucell.orbital_fn + ntype, file_orb.begin(), [](const std::string& file) {
        return PARAM.inp.orbital_dir + file;
    });
    orb_.build(ntype, file_orb.data());

    std::string file_alpha = PARAM.inp.orbital_dir + ucell.descriptor_file;
    alpha_.build(1, &file_alpha);

    double rmax = std::max(orb_.rcut_max(), alpha_.rcut_max());
    double cutoff = 2.0 * rmax;
    int nr = static_cast<int>(rmax / lcao_dr) + 1;

    orb_.set_uniform_grid(true,nr,cutoff,'i',true);
    alpha_.set_uniform_grid(true,nr,cutoff,'i',true);

    overlap_orb_alpha_.tabulate(orb_, alpha_, 'S', nr, cutoff);

    return;
}

void test_deepks::setup_kpt()
{
    this->kv.set("KPT",
                 PARAM.input.nspin,
                 ucell.G,
                 ucell.latvec,
                 PARAM.sys.gamma_only_local,
                 GlobalV::ofs_running,
                 GlobalV::ofs_warning);
}

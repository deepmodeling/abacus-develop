#include "klist.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_basis/module_ao/ORB_read.h"
#include "module_cell/module_neighbor/sltk_atom_arrange.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
//#include "parallel_orbitals.h"

#include "../LCAO_deepks.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Test_Deepks
{
extern Grid_Driver GridD;
}

namespace GlobalC
{
  extern LCAO_Deepks ld;
}

class test_deepks
{

  public:
    test_deepks();
    ~test_deepks();

    LCAO_Orbitals ORB;

    RadialCollection orb_;
    RadialCollection alpha_;
    TwoCenterIntegrator overlap_orb_alpha_;

    UnitCell ucell;

    Parallel_Orbitals ParaO;
    Test_Deepks::K_Vectors kv;
    // LCAO_Deepks ld;

    int failed_check = 0;
    int total_check = 0;

    int my_rank = 0;

    double lcao_ecut = 0; // (Ry)
    double lcao_dk = 0.01;
    double lcao_dr = 0.01;
    double lcao_rmax = 30; // (a.u.)

    int out_mat_r = 0;

    int lmax = 2;
    int ntype = 0;

    std::vector<ModuleBase::matrix> dm;
    std::vector<ModuleBase::ComplexMatrix> dm_k;

    std::vector<std::vector<double>> dm_new;
    std::vector<std::vector<std::complex<double>>> dm_k_new;

    elecstate::DensityMatrix<double, double>* p_elec_DM;
    elecstate::DensityMatrix<std::complex<double>, double>* p_elec_DM_k;

    // preparation
    void preparation();
    void set_parameters(); // set some global variables
    void setup_cell();

    void count_ntype(); // from STRU, count types of elements
    void set_ekcut();   // from LCAO files, read and set ekcut

    void prep_neighbour();
    void setup_kpt();
    void set_orbs();

    // tranfer Matrix into vector<T>
    void set_dm_new();
    void set_dm_k_new();

    // tranfer vector<T> into DensityMatrix
    void set_p_elec_DM();
    void set_p_elec_DM_k();

    // checking
    void check_dstable();
    void check_phialpha();

    void read_dm();
    void read_dm_k(const int nks);

    void check_pdm();
    void check_gdmx();

    void check_descriptor();
    void check_gvx();

    void check_edelta();

    // calculate H_V_delta
    void cal_H_V_delta();
    void cal_H_V_delta_k();

    void check_e_deltabands();
    void check_f_delta_and_stress_delta();

    // compares numbers stored in two files
    void compare_with_ref(const std::string f1, const std::string f2);
};

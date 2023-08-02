#include "gtest/gtest.h"
#include "../overlap_new.h"

// mock of OperatorLCAO

template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::Operator(){}

template<typename FPTYPE, typename Device>
hamilt::Operator<FPTYPE, Device>::~Operator(){}

template<typename FPTYPE, typename Device>
typename hamilt::Operator<FPTYPE, Device>::hpsi_info hamilt::Operator<FPTYPE, Device>::hPsi(hpsi_info&) const 
{
    return hpsi_info(nullptr, 0, nullptr);
}

template<typename FPTYPE, typename Device>
void hamilt::Operator<FPTYPE, Device>::init(const int ik_in) 
{
    return;
}

template<typename FPTYPE, typename Device>
void hamilt::Operator<FPTYPE, Device>::add(Operator* next) 
{
    return;
}

template<typename FPTYPE, typename Device>
FPTYPE* hamilt::Operator<FPTYPE, Device>::get_hpsi(const hpsi_info& info) const
{
    return nullptr;
}

template class hamilt::Operator<double, psi::DEVICE_CPU>;
template class hamilt::Operator<std::complex<double>, psi::DEVICE_CPU>;

// mock of OperatorLCAO
template<typename TK>
void hamilt::OperatorLCAO<TK>::init(const int ik_in)
{
    return;
}
template<typename TK>
void hamilt::OperatorLCAO<TK>::get_hs_pointers()
{
    return;
}
template class hamilt::OperatorLCAO<double>;
template class hamilt::OperatorLCAO<std::complex<double>>;

// mock of ORB_gen_tables and LCAO_Orbitals
#include "module_basis/module_ao/ORB_gen_tables.h"
const ORB_gen_tables& ORB_gen_tables::get_const_instance()
{
    static ORB_gen_tables instance;
    return instance;
}
ORB_gen_tables::ORB_gen_tables() {}
ORB_gen_tables::~ORB_gen_tables() {}
ORB_gaunt_table::ORB_gaunt_table() {}
ORB_gaunt_table::~ORB_gaunt_table() {}
ORB_table_phi::ORB_table_phi() {}
ORB_table_phi::~ORB_table_phi() {}
ORB_table_alpha::ORB_table_alpha() {}
ORB_table_alpha::~ORB_table_alpha() {}
ORB_table_beta::ORB_table_beta() {}
ORB_table_beta::~ORB_table_beta() {}
// mock of snap_psipsi
void ORB_gen_tables::snap_psipsi(
    const LCAO_Orbitals &orb,
    double olm[],
    const int &job, ///<[in]0 for matrix element of either S or T, 1 for its derivatives
    const char &dtype, ///<[in] derivative type, 'S' for overlap, 'T' for kinetic energy, 'D' for descriptor in deepks
    const ModuleBase::Vector3<double> &R1,
    const int &I1,
    const int &l1,
    const int &m1,
    const int &n1,
    const ModuleBase::Vector3<double> &R2,
    const int &I2,
    const int &l2,
    const int &m2,
    const int &n2,
    bool cal_syns,
    double dmax)const
{
    olm[0] = 1.0;
}

#include "module_basis/module_ao/ORB_read.h"
const LCAO_Orbitals& LCAO_Orbitals::get_const_instance()
{
    static LCAO_Orbitals instance;
    return instance;
}
LCAO_Orbitals::LCAO_Orbitals() {}
LCAO_Orbitals::~LCAO_Orbitals() {}

// mock find_atom() function
void Grid_Driver::Find_atom(const UnitCell& ucell,
                            const ModuleBase::Vector3<double>& tau,
                            const int& T,
                            const int& I,
                            AdjacentAtomInfo* adjs)
{
    adjs->adj_num = ucell.nat - 1 ;
    adjs->adjacent_tau.resize(ucell.nat);
    adjs->ntype.resize(ucell.nat, 0);
    adjs->natom.resize(ucell.nat);
    adjs->box.resize(ucell.nat);
    for(int iat = 0;iat<ucell.nat;iat++)
    {
        adjs->natom[iat] = iat;
        adjs->box[iat].x = 0;
        adjs->box[iat].y = 0;
        adjs->box[iat].z = 0;
        adjs->adjacent_tau[iat] = ucell.get_tau(iat);
    }
}
Grid::Grid(const int &test_grid_in):test_grid(test_grid_in)
{}
Grid::~Grid(){}
Grid_Driver::Grid_Driver(const int &test_d_in, 
		const int &test_gd_in, 
		const int &test_grid_in) :Grid(test_grid_in), test_deconstructor(test_d_in), test_grid_driver(test_gd_in) {}
Grid_Driver::~Grid_Driver() {}
// simple test for OverlapNew class
// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 10;
int test_nw = 10;
class OverlapNewTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        // MPI parallel settings
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        // set up a unitcell, with one element and test_size atoms, each atom has test_nw orbitals
        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
        ucell.atoms[0].tau = new ModuleBase::Vector3<double>[ucell.nat];
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        ucell.iat2iwt.resize(ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
            ucell.iat2iwt[iat] = iat * test_nw;
        }
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l = new int[test_nw];
        ucell.atoms[0].iw2m = new int[test_nw];
        ucell.atoms[0].iw2n = new int[test_nw];
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 0;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        init_parav();
        // set up a HContainer with ucell
        SR = new hamilt::HContainer<double>(paraV);
    }

    void TearDown() override
    {
        delete SR;
        delete paraV;
        delete[] ucell.atoms[0].tau;
        delete[] ucell.atoms[0].iw2l;
        delete[] ucell.atoms[0].iw2m;
        delete[] ucell.atoms[0].iw2n;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
        delete[] ucell.iat2ia;
    }

#ifdef __MPI
    void init_parav()
    {
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->set_block_size(2/* nb_2d set to be 2*/);
        paraV->set_proc_dim(dsize, 0);
        paraV->mpi_create_cart(MPI_COMM_WORLD);
        paraV->set_local2global(global_row, global_col, ofs_running, ofs_running);
        int lr = paraV->get_row_size();
        int lc = paraV->get_col_size();
        paraV->set_desc(global_row, global_col, lr);
        paraV->set_global2local(global_row, global_col, true, ofs_running);
        paraV->set_atomic_trace(ucell.iat2iwt.data(), test_size, global_row);
    }
#else
    void init_parav()
    {}
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* SR;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

// using TEST_F to test OverlapNew
TEST_F(OverlapNewTest, constructHR)
{
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    std::vector<double> hk(paraV->get_row_size() * paraV->get_col_size(), 0.0);
    Grid_Driver gd(0,0,0);
    hamilt::OverlapNew<hamilt::OperatorLCAO<double>, double> op(
        nullptr, 
        kvec_d_in, 
        SR, 
        hk.data(), 
        &ucell, 
        &gd,
        paraV
    );
    op.contributeHR();
    // check the value of SR
    for (int iap = 0; iap < SR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = SR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        auto indexes1 = paraV->get_indexes_row(iat1);
        auto indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_EQ(tmp.get_HR_values(0, 0, 0).get_pointer()[i], 1.0);
        }
    }
    // calculate SK
    op.contributeHk(0);
    // check the value of SK
    for (int i = 0; i < paraV->get_row_size() * paraV->get_col_size(); ++i)
    {
        EXPECT_EQ(hk[i], 1.0);
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
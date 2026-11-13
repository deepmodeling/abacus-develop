#include <gtest/gtest.h>
#include "source_base/atom_basis_view.h"
#include "source_base/basis_flattened_indexing.h"
#include "source_base/element_basis_index.h"
#include "source_hamilt/module_hcontainer/transfer.h"
#include "source_hamilt/module_hcontainer/hcontainer.h"
#include <chrono>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#endif

// test_size is the number of atoms in the unitcell
// modify test_size to test different size of unitcell
int test_size = 2;
int test_nw = 10;

class TransferTest : public ::testing::Test
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
        ucell.iat2it.resize(ucell.nat);
        ucell.iat2ia.resize(ucell.nat);
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.set_iat2iwt(1);
        init_parav();
        // set up a HContainer with ucell
        HR_para = new hamilt::HContainer<double>(ucell, paraV);
    }

    void TearDown() override
    {
        delete HR_para;
        delete paraV;
        delete[] ucell.atoms;
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 2;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->para_2d().init(global_row, global_col, nb, MPI_COMM_WORLD, false);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {}
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR_para = nullptr;
    Parallel_Orbitals *paraV;

    int dsize;
    int my_rank = 0;
};

TEST_F(TransferTest, serialToPara)
{
// get rank of process
#ifdef __MPI

    hamilt::HContainer<double>* HR_serial = nullptr;

// initialize HR_serial
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // if the master process, calculate the value of HR_serial and send to other processes
    if (my_rank == 0)
    {
        HR_serial = new hamilt::HContainer<double>(ucell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < HR_serial->size_atom_pairs(); i++)
        {
            hamilt::AtomPair<double>& atom_pair = HR_serial->get_atom_pair(i);
            int atom_i = atom_pair.get_atom_i();
            int atom_j = atom_pair.get_atom_j();
            //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
            auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
            double* data = atom_pair.get_pointer(0);
            for(int k = 0; k < test_nw; k++)
            {
                for(int l = 0; l < test_nw; l++)
                {
                    *data = value(k, l);
                    ++data;
                }
            }
        }
    }
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();
    hamilt::transferSerial2Parallels(*HR_serial, HR_para, 0);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    // check data in HR_para
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                EXPECT_NEAR(*data, value(row_indexes[k], col_indexes[l]), 1e-10);
                ++data;
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    if(my_rank == 0)
    {
        delete HR_serial;
    }

    std::cout <<"rank: "<<my_rank<< " HR init time: " << elapsed_time0.count()<<" transfer_s2p time: "<<elapsed_time1.count()<<" check time: "<<elapsed_time2.count()<<" seconds." << std::endl;
#endif
}

#ifdef __MPI
TEST_F(TransferTest, BasisFlattenedLayoutMatchesParallelOrbitals)
{
    ModuleBase::Element_Basis_Index::Range range(1);
    range[0].resize(2);
    range[0][0].N = 1;
    range[0][0].M = 1;
    range[0][1].N = 3;
    range[0][1].M = 3;
    const ModuleBase::Element_Basis_Index::IndexLNM hierarchy =
        ModuleBase::Element_Basis_Index::construct_index(range);
    const std::vector<int> atom_types(test_size, 0);
    const BasisFlattenedIndexing basis(hierarchy, atom_types, 1);

    hamilt::HContainer<double> from_basis(basis, paraV->para_2d());
    for (int iat = 0; iat < test_size; ++iat)
    {
        const AtomBasisView basis_rows = AtomBasisView::rows(basis, paraV->para_2d(), iat);
        const AtomBasisView basis_cols = AtomBasisView::cols(basis, paraV->para_2d(), iat);
        const std::vector<int> indexes_row = paraV->get_indexes_row(iat);
        const std::vector<int> indexes_col = paraV->get_indexes_col(iat);
        const int begin_row_atom = paraV->atom_begin_row[iat];
        const int begin_col_atom = paraV->atom_begin_col[iat];
        const int nrow_atom = paraV->get_nrow_atom(iat);
        const int ncol_atom = paraV->get_ncol_atom(iat);
        const int indexes_row_size = static_cast<int>(indexes_row.size());
        const int indexes_col_size = static_cast<int>(indexes_col.size());
        EXPECT_EQ(basis_rows.local0(), begin_row_atom);
        EXPECT_EQ(basis_cols.local0(), begin_col_atom);
        EXPECT_EQ(basis_rows.nloc(), nrow_atom);
        EXPECT_EQ(basis_cols.nloc(), ncol_atom);
        EXPECT_EQ(indexes_row_size, basis_rows.nloc());
        EXPECT_EQ(indexes_col_size, basis_cols.nloc());
        for (int k = 0; k < basis_rows.nloc(); ++k)
        {
            EXPECT_EQ(basis_rows[k], indexes_row[k]);
        }
        for (int k = 0; k < basis_cols.nloc(); ++k)
        {
            EXPECT_EQ(basis_cols[k], indexes_col[k]);
        }
    }

    ASSERT_EQ(from_basis.size_atom_pairs(), HR_para->size_atom_pairs());
    EXPECT_EQ(from_basis.get_nbasis(), HR_para->get_nbasis());
    for (int iap = 0; iap < from_basis.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& basis_pair = from_basis.get_atom_pair(iap);
        const hamilt::AtomPair<double>* orbital_pair
            = HR_para->find_pair(basis_pair.get_atom_i(), basis_pair.get_atom_j());
        ASSERT_NE(orbital_pair, nullptr);
        EXPECT_EQ(basis_pair.get_begin_row(), orbital_pair->get_begin_row());
        EXPECT_EQ(basis_pair.get_begin_col(), orbital_pair->get_begin_col());
        EXPECT_EQ(basis_pair.get_row_size(), orbital_pair->get_row_size());
        EXPECT_EQ(basis_pair.get_col_size(), orbital_pair->get_col_size());
    }

    hamilt::HContainer<double> inserted(basis, paraV->para_2d(), nullptr);
    const ModuleBase::Vector3<int> R_index(0, 0, 0);
    for (int iat1 = 0; iat1 < test_size; ++iat1)
    {
        for (int iat2 = 0; iat2 < test_size; ++iat2)
        {
            inserted.insert_pair(iat1, iat2, R_index);
        }
    }
    EXPECT_EQ(inserted.size_atom_pairs(), HR_para->size_atom_pairs());

    const std::vector<int> ijr_info = HR_para->get_ijr_info();
    hamilt::HContainer<double> from_sparse(basis, paraV->para_2d(), nullptr, &ijr_info);
    ASSERT_EQ(from_sparse.size_atom_pairs(), HR_para->size_atom_pairs());
    EXPECT_EQ(from_sparse.get_nbasis(), HR_para->get_nbasis());
    for (int iap = 0; iap < from_sparse.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& basis_pair = from_sparse.get_atom_pair(iap);
        const hamilt::AtomPair<double>* orbital_pair
            = HR_para->find_pair(basis_pair.get_atom_i(), basis_pair.get_atom_j());
        ASSERT_NE(orbital_pair, nullptr);
        EXPECT_EQ(basis_pair.get_begin_row(), orbital_pair->get_begin_row());
        EXPECT_EQ(basis_pair.get_begin_col(), orbital_pair->get_begin_col());
        EXPECT_EQ(basis_pair.get_row_size(), orbital_pair->get_row_size());
        EXPECT_EQ(basis_pair.get_col_size(), orbital_pair->get_col_size());
    }

    hamilt::HContainer<double> from_basis_inserted(basis, paraV->para_2d());
    hamilt::HContainer<double> from_orbital_inserted(paraV);
    for (int iap = 0; iap < HR_para->size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& source_pair = HR_para->get_atom_pair(iap);
        const int iat1 = source_pair.get_atom_i();
        const int iat2 = source_pair.get_atom_j();
        EXPECT_TRUE(from_basis_inserted.owns_atom_pair(iat1, iat2));
        EXPECT_TRUE(from_orbital_inserted.owns_atom_pair(iat1, iat2));
        for (int ir = 0; ir < source_pair.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> R_index = source_pair.get_R_index(ir);
            from_basis_inserted.insert_pair(iat1, iat2, R_index);
            from_orbital_inserted.insert_pair(iat1, iat2, R_index);
        }
    }
    EXPECT_FALSE(from_basis_inserted.owns_atom_pair(-1, 0));
    EXPECT_FALSE(from_orbital_inserted.owns_atom_pair(0, test_size));
    ASSERT_EQ(from_basis_inserted.size_atom_pairs(), HR_para->size_atom_pairs());
    ASSERT_EQ(from_orbital_inserted.size_atom_pairs(), HR_para->size_atom_pairs());
    for (int iap = 0; iap < from_basis_inserted.size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<double>& basis_pair = from_basis_inserted.get_atom_pair(iap);
        const hamilt::AtomPair<double>* orbital_pair
            = from_orbital_inserted.find_pair(basis_pair.get_atom_i(), basis_pair.get_atom_j());
        ASSERT_NE(orbital_pair, nullptr);
        EXPECT_EQ(basis_pair.get_begin_row(), orbital_pair->get_begin_row());
        EXPECT_EQ(basis_pair.get_begin_col(), orbital_pair->get_begin_col());
        EXPECT_EQ(basis_pair.get_row_size(), orbital_pair->get_row_size());
        EXPECT_EQ(basis_pair.get_col_size(), orbital_pair->get_col_size());
        EXPECT_EQ(basis_pair.get_R_size(), orbital_pair->get_R_size());
    }
}

#endif

TEST_F(TransferTest, paraToSerial)
{
    // get rank of process
#ifdef __MPI

    hamilt::HContainer<double>* HR_serial = nullptr;

    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // initialize HR_para
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                *data = value(row_indexes[k], col_indexes[l]);
                ++data;
            }
        }
    }

    // initialize HR_serial
    // if the master process, calculate the value of HR_serial and send to other processes
    if (my_rank == 0)
    {
        HR_serial = new hamilt::HContainer<double>(ucell);
    }
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();
    hamilt::transferParallels2Serial(*HR_para, HR_serial, 0);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    start_time = std::chrono::high_resolution_clock::now();
    // check data in HR_serial
    if(my_rank == 0)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < HR_serial->size_atom_pairs(); i++)
        {
            hamilt::AtomPair<double>& atom_pair = HR_serial->get_atom_pair(i);
            int atom_i = atom_pair.get_atom_i();
            int atom_j = atom_pair.get_atom_j();
            //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
            auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
            double* data = atom_pair.get_pointer(0);
            for(int k = 0; k < test_nw; k++)
            {
                for(int l = 0; l < test_nw; l++)
                {
                    EXPECT_NEAR(*data , value(k, l), 1e-10);
                    ++data;
                }
            }
        }
        delete HR_serial;
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout <<"rank: "<<my_rank<< " HR init time: " << elapsed_time0.count()<<" transfer_p2s time: "<<elapsed_time1.count()<<" check time: "<<elapsed_time2.count()<<" seconds." << std::endl;

    // now assume the HR_serial is empty and do the same test
    Parallel_Orbitals serialV;
    std::ofstream ofs("test.log");
    serialV.para_2d().set_serial(test_size*test_nw, test_size*test_nw);
    serialV.set_atomic_trace(ucell.get_iat2iwt(), test_size, test_size*test_nw);
    hamilt::HContainer<double>* HR_serial2;
    if(my_rank == 0)
    {
        HR_serial2 = new hamilt::HContainer<double>(&serialV);
    }
    hamilt::gatherParallels(*HR_para, HR_serial2, 0);
    if(my_rank == 0)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i = 0; i < HR_serial2->size_atom_pairs(); i++)
        {
            hamilt::AtomPair<double>& atom_pair = HR_serial2->get_atom_pair(i);
            int atom_i = atom_pair.get_atom_i();
            int atom_j = atom_pair.get_atom_j();
            //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
            auto value = [&](int k, int l) -> double {return (((atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
            double* data = atom_pair.get_pointer(0);
            for(int k = 0; k < test_nw; k++)
            {
                for(int l = 0; l < test_nw; l++)
                {
                    EXPECT_NEAR(*data , value(k, l), 1e-10);
                    ++data;
                }
            }
        }
        delete HR_serial2;
    }
#endif
}

TEST_F(TransferTest, serialAllToParaAll)
{
// get rank of process
#ifdef __MPI

// initialize HR_serial
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // for each process, calculate the value of HR_serial and send to other processes
    hamilt::HContainer<double> HR_serial(ucell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_serial.size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_serial.get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        for(int k = 0; k < test_nw; k++)
        {
            for(int l = 0; l < test_nw; l++)
            {
                *data = value(k, l);
                ++data;
            }
        }
    }
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    start_time = std::chrono::high_resolution_clock::now();
    hamilt::transferSerials2Parallels(HR_serial, HR_para);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    // check data in HR_para
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l)*dsize;};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                EXPECT_NEAR(*data, value(row_indexes[k], col_indexes[l]), 1e-10);
                ++data;
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    std::cout <<"rank: "<<my_rank<< " HR init time: " << elapsed_time0.count()<<" transfer_sa2pa time: "<<elapsed_time1.count()<<" check time: "<<elapsed_time2.count()<<" seconds." << std::endl;
#endif
}

TEST_F(TransferTest, paraAllToserialAll)
{
// get rank of process
#ifdef __MPI
    std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
    // for each process, calculate the value of HR_para and send to other processes
    // initialize HR_para
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_para->size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_para->get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        auto row_indexes = paraV->get_indexes_row(atom_i);
        auto col_indexes = paraV->get_indexes_col(atom_j);
        for(int k = 0; k < row_indexes.size(); k++)
        {
            for(int l = 0; l < col_indexes.size(); l++)
            {
                *data = value(row_indexes[k], col_indexes[l]);
                ++data;
            }
        }
    }
    // prepare memory for HR_serial
    hamilt::HContainer<double> HR_serial(ucell);
    std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time0 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    hamilt::transferParallels2Serials(*HR_para, &HR_serial);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time1 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    // check data in HR_serial
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i = 0; i < HR_serial.size_atom_pairs(); i++)
    {
        hamilt::AtomPair<double>& atom_pair = HR_serial.get_atom_pair(i);
        int atom_i = atom_pair.get_atom_i();
        int atom_j = atom_pair.get_atom_j();
        //lambda function to calculate value of array: (atom_i*test_size+atom_j+k)*test_nw + l
        auto value = [&](int k, int l) -> double {return ((double(atom_i*test_nw+k)*test_size+atom_j)*test_nw + l);};
        double* data = atom_pair.get_pointer(0);
        for(int k = 0; k < test_nw; k++)
        {
            for(int l = 0; l < test_nw; l++)
            {
                EXPECT_NEAR(*data , value(k, l), 1e-10);
                ++data;
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time2 = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

    std::cout <<"rank: "<<my_rank<< " HR init time: " << elapsed_time0.count()<<" transfer_pa2sa time: "<<elapsed_time1.count()<<" check time: "<<elapsed_time2.count()<<" seconds." << std::endl;
#endif
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

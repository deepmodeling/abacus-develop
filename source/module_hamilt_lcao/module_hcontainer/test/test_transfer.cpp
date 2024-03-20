#include "gtest/gtest.h"
#include "../transfer.h"
#include "../hcontainer.h"
#include <chrono>
#ifdef __MPI
#include <mpi.h>
#include "../hcontainer_funcs.h"
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
        ucell.iat2it = new int[ucell.nat];
        ucell.iat2ia = new int[ucell.nat];
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
    serialV.set_global2local(test_size*test_nw, test_size*test_nw, false, ofs);
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
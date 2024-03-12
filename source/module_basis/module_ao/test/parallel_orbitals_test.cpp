#include "gtest/gtest.h"
#include "../parallel_orbitals.h"

//--------------------------------------------------------------
//      unit test of class "Parallel_Orbitals"
//--------------------------------------------------------------
/**
 * Test functions:
 * - set_atomic_trace
 * - get_col_size
 * - get_row_size
 * - get_col_size(iat)
 * - get_row_size(iat)
 * - get_indexes_row
 * - get_indexes_col
 * - get_indexes_row(iat)
 * - get_indexes_col(iat)
 * 
 * the test framework is based on parallel_2d_test.cpp
*/
class TestParaO : public testing::Test
{
protected:
    int dsize;
    int my_rank = 0;
    std::vector<std::pair<int, int>> sizes{ {50, 50} , {60, 60}};
    std::vector<int> nat{ 10, 5};
    std::vector<int> nbs{ 1,2,3 };
    std::ofstream ofs_running;
#ifdef __MPI
    void SetUp() override
    {
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        this->ofs_running.open("log" + std::to_string(my_rank) + ".txt");
        ofs_running << "dsize(nproc) = " << dsize << std::endl;
        ofs_running << "my_rank = " << my_rank << std::endl;
    }
    void TearDown() override
    {
        ofs_running.close();
    }
#endif
};

#ifdef __MPI
TEST_F(TestParaO, Divide2D)
{
    for (auto& size : sizes)
    {
        int gr = size.first;
        int gc = size.second;
        for (auto nb : nbs)
        {
            Parallel_Orbitals po;
            po.set_block_size(nb);
            EXPECT_EQ(po.get_block_size(), nb);

            for (auto mode : { 0,1 })
            {
                //1. set dim0 and dim1
                po.set_proc_dim(dsize, mode);
                EXPECT_EQ(po.dim0 * po.dim1, dsize);
                if (mode)EXPECT_LE(po.dim1, po.dim0);
                else EXPECT_LE(po.dim0, po.dim1);

                //2. mpi_create_cart
                po.mpi_create_cart(MPI_COMM_WORLD);
                EXPECT_NE(po.comm_2D, MPI_COMM_NULL);

                //3. set_local2global and local sizes
                po.set_local2global(gr, gc, ofs_running, ofs_running);
                int lr = po.get_row_size();
                int lc = po.get_col_size();
                EXPECT_EQ(lr * lc, po.get_local_size());
                auto cal_lsize = [](const int& gsize, const int& nb, const int& np, const int& pcoord) -> int
                    {
                        int nblock = gsize / nb;
                        return nblock / np * nb + static_cast<int>(nblock % np > pcoord) * nb //full blocks' contribution
                            + static_cast<int>(nblock % np == pcoord) * (gsize % nb);   // the last block's contribution
                    };
                EXPECT_EQ(lr, cal_lsize(gr, nb, po.dim0, po.coord[0]));
                EXPECT_EQ(lc, cal_lsize(gc, nb, po.dim1, po.coord[1]));

                //4. set_desc
                po.set_desc(gr, gc, lr);
                EXPECT_EQ(po.desc[0], 1);
                EXPECT_EQ(po.desc[1], po.blacs_ctxt);
                EXPECT_EQ(po.desc[2], gr);
                EXPECT_EQ(po.desc[3], gc);
                EXPECT_EQ(po.desc[4], po.get_block_size());
                EXPECT_EQ(po.desc[5], po.get_block_size());
                EXPECT_EQ(po.desc[6], 0);
                EXPECT_EQ(po.desc[7], 0);
                EXPECT_EQ(po.desc[8], lr);

                //5. set_global2local
                po.set_global2local(gr, gc, true, ofs_running);
                auto sum_array = [&po](const int& gr, const int& gc) -> std::pair<int, int>
                {
                    int sum_row = 0; int sum_col = 0;
                    for (int i = 0; i < gr; ++i)
                        sum_row += po.global2local_row(i);
                    for (int i = 0; i < gc; ++i)
                        sum_col += po.global2local_col(i);
                    return { sum_row, sum_col };
                };
                std::pair<int, int> sumrc = sum_array(gr, gc);
                EXPECT_EQ(std::get<0>(sumrc), lr * (lr - 1) / 2 - (gr - lr));
                EXPECT_EQ(std::get<1>(sumrc), lc * (lc - 1) / 2 - (gc - lc));
                for (int i = 0;i < lr;++i)
                    for (int j = 0;j < lc;++j)
                        EXPECT_TRUE(po.in_this_processor(po.local2global_row(i), po.local2global_col(j)));
                
                //6. set_atomic_trace
                for(auto nat0 : nat)
                {
                    EXPECT_EQ(gr, gc);
                    std::vector<int> iat2iwt(nat0);
                    int nw = gr / nat0;
                    for (int i = 0; i < nat0; ++i)
                    {
                        iat2iwt[i] = i * nw;
                    }
                    po.set_atomic_trace(iat2iwt.data(), nat0, gr);
                    auto global_row_array = po.get_indexes_row();
                    auto global_col_array = po.get_indexes_col();
                    int local_index_trace_row = 0;
                    int local_index_trace_col = 0;
                    // check get_col_size(iat) and get_row_size(iat)
                    for (int i = 0; i < nat0; ++i)
                    {
                        auto atomic_row_array = po.get_indexes_row(i);
                        auto atomic_col_array = po.get_indexes_col(i);
                        EXPECT_EQ(po.get_col_size(i), atomic_col_array.size());
                        EXPECT_EQ(po.get_row_size(i), atomic_row_array.size());
                        for (int j = 0; j < atomic_row_array.size(); ++j)
                        {
                            //check global_index == global_index
                            EXPECT_EQ(atomic_row_array[j]+iat2iwt[i], global_row_array[local_index_trace_row]);
                            //check local_index == local_index
                            EXPECT_EQ(local_index_trace_row, po.global2local_row(atomic_row_array[j]+iat2iwt[i]));
                            local_index_trace_row++;
                        }
                        for (int j = 0; j < atomic_col_array.size(); ++j)
                        {
                            //check global_index == global_index
                            EXPECT_EQ(atomic_col_array[j]+iat2iwt[i], global_col_array[local_index_trace_col]);
                            //check local_index == local_index
                            EXPECT_EQ(local_index_trace_col, po.global2local_col(atomic_col_array[j]+iat2iwt[i]));
                            local_index_trace_col++;
                        }
                        
                    }
                }
            }
        }
    }
}
#else
TEST_F(TestParaO, Serial)
{
    for (auto& size : sizes)
    {
        int gr = size.first;
        int gc = size.second;

        Parallel_Orbitals po;

        //1. set dim0 and dim1
        po.set_proc_dim(1);
        EXPECT_EQ(po.dim0 * po.dim1, 1);

        //2. set_serial
        po.set_serial(gr, gc);
        EXPECT_EQ(po.get_row_size(), gr);
        EXPECT_EQ(po.get_col_size(), gc);
        EXPECT_EQ(po.get_local_size(), gr * gc);

        //3. set_global2local
        po.set_global2local(gr, gc, false, ofs_running);
        for (int i = 0;i < gr;++i)
            EXPECT_EQ(po.global2local_row(i), i);
        for (int i = 0;i < gc;++i)
            EXPECT_EQ(po.global2local_col(i), i);
        //6. set_atomic_trace
        for(auto nat0 : nat)
        {
            EXPECT_EQ(gr, gc);
            std::vector<int> iat2iwt(nat0);
            int nw = gr / nat0;
            for (int i = 0; i < nat0; ++i)
            {
                iat2iwt[i] = i * nw;
            }
            po.set_atomic_trace(iat2iwt.data(), nat0, gr);
            EXPECT_EQ(po.get_col_size(), gr);
            EXPECT_EQ(po.get_row_size(), gr);
            // check get_col_size(iat) and get_row_size(iat)
            for (int i = 0; i < nat0; ++i)
            {
                std::cout<<__FILE__<<__LINE__<<" i = "<<i<<" size = "<<po.get_row_size(i)<<" "<<po.get_col_size(i)<<std::endl;
                //EXPECT_EQ(po.get_col_size(i), nw);
                //EXPECT_EQ(po.get_row_size(i), nw);
            }
        }
    }
}
#endif

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


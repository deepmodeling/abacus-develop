#include "gtest/gtest.h"
#include "../parallel_orbitals.h"

//--------------------------------------------------------------
//      unit test of class "Parallel_Orbitals"
//--------------------------------------------------------------
/**
 * Test functions:
 * - set_atomic_trace
 * - ncloc
 * - nrloc
 * - get_ncol_atom(iat)
 * - get_nrow_atom(iat)
 * - get_indexes_row
 * - get_indexes_col
 * - get_indexes_row(iat)
 * - get_indexes_col(iat)
 * 
 * the test framework follows the block-cyclic layout tests
*/
class TestParaO : public testing::Test
{
protected:
    int dsize;
    int my_rank = 0;
    std::vector<std::pair<int, int>> sizes{ {50, 50} , {60, 60}};
    std::vector<int> nat{ 10, 5};
    std::vector<int> nbs{ 1,2,3 };
#ifdef __MPI
    void SetUp() override
    {
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    }
    void TearDown() override
    {
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

            for (auto mode : { 0,1 })
            {
                po.para_2d().init(gr, gc, nb, MPI_COMM_WORLD, mode);
                EXPECT_EQ(po.para_2d().nb(), nb);

                //1. dim0 and dim1
                EXPECT_EQ(po.para_2d().nrproc() * po.para_2d().ncproc(), dsize);
                if (mode)EXPECT_LE(po.para_2d().ncproc(), po.para_2d().nrproc());
                else EXPECT_LE(po.para_2d().nrproc(), po.para_2d().ncproc());

                //2. comm_2D
                //EXPECT_NE(po.comm_2D, MPI_COMM_NULL);

                //3. local2global and local sizes
                int lr = po.para_2d().nrloc();
                int lc = po.para_2d().ncloc();
                EXPECT_EQ(lr * lc, po.para_2d().size());
                auto cal_lsize = [](const int& gsize, const int& nb, const int& np, const int& pcoord) -> int
                    {
                        int nblock = gsize / nb;
                        return nblock / np * nb + static_cast<int>(nblock % np > pcoord) * nb //full blocks' contribution
                            + static_cast<int>(nblock % np == pcoord) * (gsize % nb);   // the last block's contribution
                    };
                EXPECT_EQ(lr, cal_lsize(gr, nb, po.para_2d().nrproc(), po.para_2d().irproc()));
                EXPECT_EQ(lc, cal_lsize(gc, nb, po.para_2d().ncproc(), po.para_2d().icproc()));

                //4. ScaLAPACK descriptor
                EXPECT_EQ(po.para_2d().desc()[0], 1);
                EXPECT_EQ(po.para_2d().desc()[1], po.para_2d().ictxt());
                EXPECT_EQ(po.para_2d().desc()[2], gr);
                EXPECT_EQ(po.para_2d().desc()[3], gc);
                EXPECT_EQ(po.para_2d().desc()[4], po.para_2d().nb());
                EXPECT_EQ(po.para_2d().desc()[5], po.para_2d().nb());
                EXPECT_EQ(po.para_2d().desc()[6], 0);
                EXPECT_EQ(po.para_2d().desc()[7], 0);
                EXPECT_EQ(po.para_2d().desc()[8], lr);

                //5. global2local
                auto sum_array = [&po](const int& gr, const int& gc) -> std::pair<int, int>
                {
                    int sum_row = 0; int sum_col = 0;
                    for (int i = 0; i < gr; ++i)
                        sum_row += po.para_2d().irloc(i);
                    for (int i = 0; i < gc; ++i)
                        sum_col += po.para_2d().icloc(i);
                    return { sum_row, sum_col };
                };
                std::pair<int, int> sumrc = sum_array(gr, gc);
                EXPECT_EQ(std::get<0>(sumrc), lr * (lr - 1) / 2 - (gr - lr));
                EXPECT_EQ(std::get<1>(sumrc), lc * (lc - 1) / 2 - (gc - lc));
                for (int i = 0;i < lr;++i)
                    for (int j = 0;j < lc;++j)
                        EXPECT_TRUE(po.para_2d().in_this_processor(po.para_2d().irglb(i), po.para_2d().icglb(j)));
                
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
                    // check get_ncol_atom(iat) and get_nrow_atom(iat)
                    for (int i = 0; i < nat0; ++i)
                    {
                        auto atomic_row_array = po.get_indexes_row(i);
                        auto atomic_col_array = po.get_indexes_col(i);
                        EXPECT_EQ(po.get_ncol_atom(i), atomic_col_array.size());
                        EXPECT_EQ(po.get_nrow_atom(i), atomic_row_array.size());
                        for (int j = 0; j < atomic_row_array.size(); ++j)
                        {
                            //check global_index == global_index
                            EXPECT_EQ(atomic_row_array[j]+iat2iwt[i], global_row_array[local_index_trace_row]);
                            //check local_index == local_index
                            EXPECT_EQ(local_index_trace_row, po.para_2d().irloc(atomic_row_array[j]+iat2iwt[i]));
                            local_index_trace_row++;
                        }
                        for (int j = 0; j < atomic_col_array.size(); ++j)
                        {
                            //check global_index == global_index
                            EXPECT_EQ(atomic_col_array[j]+iat2iwt[i], global_col_array[local_index_trace_col]);
                            //check local_index == local_index
                            EXPECT_EQ(local_index_trace_col, po.para_2d().icloc(atomic_col_array[j]+iat2iwt[i]));
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
        //2. set_serial
        po.para_2d().set_serial(gr, gc);
        EXPECT_EQ(po.para_2d().nrproc() * po.para_2d().ncproc(), 1);
        EXPECT_EQ(po.para_2d().nrloc(), gr);
        EXPECT_EQ(po.para_2d().ncloc(), gc);
        EXPECT_EQ(po.para_2d().size(), gr * gc);

        //3. global2local
        for (int i = 0;i < gr;++i)
            EXPECT_EQ(po.para_2d().irloc(i), i);
        for (int i = 0;i < gc;++i)
            EXPECT_EQ(po.para_2d().icloc(i), i);
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
            EXPECT_EQ(po.para_2d().ncloc(), gr);
            EXPECT_EQ(po.para_2d().nrloc(), gr);
            // check get_ncol_atom(iat) and get_nrow_atom(iat)
            for (int i = 0; i < nat0; ++i)
            {
                std::cout<<__FILE__<<__LINE__<<" i = "<<i<<" size = "<<po.get_nrow_atom(i)<<" "<<po.get_ncol_atom(i)<<std::endl;
                //EXPECT_EQ(po.get_ncol_atom(i), nw);
                //EXPECT_EQ(po.get_nrow_atom(i), nw);
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


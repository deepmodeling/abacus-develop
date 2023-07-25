#include "gtest/gtest.h"
#include "../parallel_2d.h"
/***********************************************************
 *      unit test of class "Parallel_2D"
 ***********************************************************/

 /* Tested functions (in order):
  *
  * - set_proc_dim
  *   set the 2D-structure of processors in each dimension.
  *
  * - mpi_create_cart (parallel)
  *   create a 2D MPI_Cart structure.
  *
  * - set_local2global (parallel)
  *   set the map from local index to global index and the local sizes.
  *
  * - set_desc (parallel)
  *   set the desc[9] of the 2D-block-cyclic distribution.
  *
  * - set_global2local
  *   set the map from global index to local index (init, reuse).
  *
  * - set_serial (serial)
  *   set the local(=global) sizes.
  *
  * - some getters:
  *  - get_row_size, get_col_size, get_local_size, get_block_size
  *  - in_this_processor
  *
  * Result check:
  * - local sizes
  * - index maps
  * - desc[9]
  ***********************************************************/
class test_para2d : public testing::Test
{
protected:
    int dsize;
    int my_rank = 0;
    std::vector<std::pair<int, int>> sizes{ {30, 35}, {49, 94}, {57, 57} };
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
TEST_F(test_para2d, Divide2D)
{
    for (auto& size : sizes)
    {
        int gr = size.first;
        int gc = size.second;
        for (auto nb : nbs)
        {
            Parallel_2D p2d;
            p2d.set_block_size(nb);
            EXPECT_EQ(p2d.get_block_size(), nb);

            for (auto mode : { 0,1 })
            {
                //1. set dim0 and dim1
                p2d.set_proc_dim(dsize, mode);
                EXPECT_EQ(p2d.dim0 * p2d.dim1, dsize);
                if (mode)EXPECT_LE(p2d.dim1, p2d.dim0);
                else EXPECT_LE(p2d.dim0, p2d.dim1);

                //2. mpi_create_cart
                p2d.mpi_create_cart(MPI_COMM_WORLD);
                EXPECT_NE(p2d.comm_2D, MPI_COMM_NULL);

                //3. set_local2global and local sizes
                p2d.set_local2global(gr, gc, ofs_running, ofs_running);
                int lr = p2d.get_row_size();
                int lc = p2d.get_col_size();
                EXPECT_EQ(lr * lc, p2d.get_local_size());
                auto cal_lsize = [](const int& gsize, const int& nb, const int& np, const int& pcoord) -> int
                    {
                        int nblock = gsize / nb;
                        return nblock / np * nb + static_cast<int>(nblock % np > pcoord) * nb //full blocks' contribution
                            + static_cast<int>(nblock % np == pcoord) * (gsize % nb);   // the last block's contribution
                    };
                EXPECT_EQ(lr, cal_lsize(gr, nb, p2d.dim0, p2d.coord[0]));
                EXPECT_EQ(lc, cal_lsize(gc, nb, p2d.dim1, p2d.coord[1]));

                //4. set_desc
                p2d.set_desc(gr, gc, lr);
                EXPECT_EQ(p2d.desc[0], 1);
                EXPECT_EQ(p2d.desc[1], p2d.blacs_ctxt);
                EXPECT_EQ(p2d.desc[2], gr);
                EXPECT_EQ(p2d.desc[3], gc);
                EXPECT_EQ(p2d.desc[4], p2d.get_block_size());
                EXPECT_EQ(p2d.desc[5], p2d.get_block_size());
                EXPECT_EQ(p2d.desc[6], 0);
                EXPECT_EQ(p2d.desc[7], 0);
                EXPECT_EQ(p2d.desc[8], lr);

                //5. set_global2local
                p2d.set_global2local(gr, gc, true, ofs_running);
                auto sum_array = [&p2d](const int& gr, const int& gc) -> std::pair<int, int>
                    {
                        int sum_row = 0; int sum_col = 0;
                        for (int i = 0; i < gr; ++i)
                            sum_row += p2d.global2local_row(i);
                        for (int i = 0; i < gc; ++i)
                            sum_col += p2d.global2local_col(i);
                        return { sum_row, sum_col };
                    };
                std::pair<int, int> sumrc = sum_array(gr, gc);
                EXPECT_EQ(std::get<0>(sumrc), lr * (lr - 1) / 2 - (gr - lr));
                EXPECT_EQ(std::get<1>(sumrc), lc * (lc - 1) / 2 - (gc - lc));
                for (int i = 0;i < lr;++i)
                    for (int j = 0;j < lc;++j)
                        EXPECT_TRUE(p2d.in_this_processor(p2d.local2global_row(i), p2d.local2global_col(j)));
            }
        }
    }
}

TEST_F(test_para2d, DescReuseCtxt)
{
    for (auto nb : nbs)
    {
        Parallel_2D p1;
        p1.set_block_size(nb);
        p1.set_proc_dim(dsize);
        p1.mpi_create_cart(MPI_COMM_WORLD);
        p1.set_local2global(sizes[0].first, sizes[0].second, ofs_running, ofs_running);
        p1.set_desc(sizes[0].first, sizes[0].second, p1.get_row_size());

        Parallel_2D p2;    // use 2 different sizes, but they can share the same ctxt
        p2.set_block_size(nb);
        p2.set_proc_dim(dsize);
        p2.comm_2D = p1.comm_2D;
        p2.blacs_ctxt = p1.blacs_ctxt;
        p2.set_local2global(sizes[1].first, sizes[1].second, ofs_running, ofs_running);
        p2.set_desc(sizes[1].first, sizes[1].second, p2.get_row_size(), false);

        EXPECT_EQ(p1.desc[1], p2.desc[1]);

        Parallel_2D p3;    // using default `set_desc`, p3 can't share the same ctxt with p1
        p3.set_block_size(nb);
        p3.set_proc_dim(dsize);
        p3.comm_2D = p1.comm_2D;
        p3.blacs_ctxt = p1.blacs_ctxt;
        p3.set_local2global(sizes[2].first, sizes[2].second, ofs_running, ofs_running);
        p3.set_desc(sizes[2].first, sizes[2].second, p3.get_row_size());

        EXPECT_NE(p1.desc[1], p3.desc[1]);
    }
}
#else
TEST_F(test_para2d, Serial)
{
    for (auto& size : sizes)
    {
        int gr = size.first;
        int gc = size.second;

        Parallel_2D p2d;

        //1. set dim0 and dim1
        p2d.set_proc_dim(1);
        EXPECT_EQ(p2d.dim0 * p2d.dim1, 1);

        //2. set_serial
        p2d.set_serial(gr, gc);
        EXPECT_EQ(p2d.get_row_size(), gr);
        EXPECT_EQ(p2d.get_col_size(), gc);
        EXPECT_EQ(p2d.get_local_size(), gr * gc);

        //3. set_global2local
        p2d.set_global2local(gr, gc, false, ofs_running);
        for (int i = 0;i < gr;++i)
            EXPECT_EQ(p2d.global2local_row(i), i);
        for (int i = 0;i < gc;++i)
            EXPECT_EQ(p2d.global2local_col(i), i);
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


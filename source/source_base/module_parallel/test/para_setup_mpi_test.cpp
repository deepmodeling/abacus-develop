#include "gtest/gtest.h"

#include "../para_setup.h"
#include "../para_collection.h"
#include "../para_tag.h"
#include "../para_world.h"

#include <mpi.h>

// These tests run under mpirun -np 4 (see para_setup_mpi_test.sh).

namespace
{
int world_rank = -1;
int world_size = -1;
}

// Single image: the esolver domain wraps the whole world and the
// cross-image domain is absent (same convention as KP_WORLD at kpar == 1).
TEST(ParaSetupMpiTest, SingleImage)
{
    Parallel::ParaCollection worlds
        = Parallel::setup_para_worlds(world_size, world_rank, /*nimage=*/1,
                                      /*bndpar=*/1, /*kpar=*/1, /*diag_np=*/1);

    const Parallel::ParaWorld& esolver = worlds.find(Parallel::ParaTag::esolver);
    EXPECT_TRUE(esolver.valid());
    EXPECT_EQ(esolver.size(), world_size);
    EXPECT_EQ(esolver.rank(), world_rank);

    const Parallel::ParaWorld& images = worlds.find(Parallel::ParaTag::images);
    EXPECT_FALSE(images.valid());

    // All solver domains span the full world in a single-image run.
    EXPECT_EQ(worlds.find(Parallel::ParaTag::pw).size(), world_size);
    EXPECT_FALSE(worlds.find(Parallel::ParaTag::kmesh).valid());
    EXPECT_EQ(worlds.find(Parallel::ParaTag::diag).size(), world_size);
    EXPECT_EQ(worlds.find(Parallel::ParaTag::rgrid).size(), world_size);
}

// Two images on 4 ranks: each esolver owns 2 ranks; the images domain
// connects corresponding ranks (rank_in_esolver) across the two images.
TEST(ParaSetupMpiTest, TwoImages)
{
    if (world_size < 4)
    {
        GTEST_SKIP() << "requires 4 MPI ranks (run via para_setup_mpi_test.sh)";
    }
    const int nimage = 2;
    int image_id = 0;
    int rank_in_esolver = 0;
    int esolver_size = 0;
    Parallel::ParaWorld esolver_world = Parallel::ParaWorld::make_mpi("esolver", MPI_COMM_NULL);
    Parallel::ParaWorld images_world = Parallel::ParaWorld::make_mpi("images", MPI_COMM_NULL);

    Parallel::split_images(world_size, world_rank, nimage,
                           image_id, rank_in_esolver, esolver_size,
                           esolver_world, images_world);

    EXPECT_EQ(esolver_size, world_size / nimage);
    EXPECT_EQ(image_id, world_rank / (world_size / nimage));
    EXPECT_EQ(rank_in_esolver, world_rank % (world_size / nimage));

    EXPECT_TRUE(esolver_world.valid());
    EXPECT_EQ(esolver_world.size(), world_size / nimage);
    EXPECT_EQ(esolver_world.rank(), rank_in_esolver);

    // Even split: the inter-image domain exists and contains one rank
    // per image, i.e. its size equals nimage.
    EXPECT_TRUE(images_world.valid());
    EXPECT_EQ(images_world.size(), nimage);
    EXPECT_EQ(images_world.rank(), image_id);
}

// Full hierarchy with two images: every solver domain must be derived
// from the esolver domain, so its size never exceeds esolver_size.
TEST(ParaSetupMpiTest, TwoImagesFullHierarchy)
{
    if (world_size < 4)
    {
        GTEST_SKIP() << "requires 4 MPI ranks (run via para_setup_mpi_test.sh)";
    }
    Parallel::ParaCollection worlds
        = Parallel::setup_para_worlds(world_size, world_rank, /*nimage=*/2,
                                      /*bndpar=*/1, /*kpar=*/1, /*diag_np=*/1);

    const Parallel::ParaWorld& esolver = worlds.find(Parallel::ParaTag::esolver);
    EXPECT_EQ(esolver.size(), world_size / 2);

    const Parallel::ParaWorld& images = worlds.find(Parallel::ParaTag::images);
    EXPECT_TRUE(images.valid());
    EXPECT_EQ(images.size(), 2);

    // Domains inside one esolver never see ranks of the other image.
    EXPECT_EQ(worlds.find(Parallel::ParaTag::pw).size(), world_size / 2);
    EXPECT_EQ(worlds.find(Parallel::ParaTag::bsame_kdiff).size(), world_size / 2);
    EXPECT_EQ(worlds.find(Parallel::ParaTag::diag).size(), world_size / 2);
    EXPECT_EQ(worlds.find(Parallel::ParaTag::rgrid).size(), world_size / 2);
}

// k-parallelism inside an image: with 2 images and kpar=2 each pool
// contains one rank, and the inter-pool domain has size kpar = 2.
TEST(ParaSetupMpiTest, TwoImagesWithKpar)
{
    if (world_size < 4)
    {
        GTEST_SKIP() << "requires 4 MPI ranks (run via para_setup_mpi_test.sh)";
    }
    Parallel::ParaCollection worlds
        = Parallel::setup_para_worlds(world_size, world_rank, /*nimage=*/2,
                                      /*bndpar=*/1, /*kpar=*/2, /*diag_np=*/1);

    const Parallel::ParaWorld& esolver = worlds.find(Parallel::ParaTag::esolver);
    EXPECT_EQ(esolver.size(), 2);

    const Parallel::ParaWorld& pw = worlds.find(Parallel::ParaTag::pw);
    EXPECT_TRUE(pw.valid());
    EXPECT_EQ(pw.size(), 1); // 2 ranks per image / kpar 2

    const Parallel::ParaWorld& kmesh = worlds.find(Parallel::ParaTag::kmesh);
    EXPECT_TRUE(kmesh.valid());
    EXPECT_EQ(kmesh.size(), 2); // one corresponding rank per pool
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}

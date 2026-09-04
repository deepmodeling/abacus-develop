#include "gtest/gtest.h"

#include "../para_collection.h"
#include "../para_kmesh_world.h"
#include "../para_tag.h"

TEST(ParaCollectionMpiTest, AssembleAndFind)
{
    Parallel::ParaCollection coll;
    coll.add(std::unique_ptr<Parallel::ParaWorld>(
        new Parallel::ParaKmeshWorld(MPI_COMM_WORLD, 1, 0, 1, 4, 1)));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));
    EXPECT_EQ(coll.size(), 2u);

    const Parallel::ParaWorld& kmesh = coll.find(Parallel::ParaTag::kmesh);
    EXPECT_TRUE(kmesh.valid());
    EXPECT_EQ(kmesh.tag(), "kmesh");

    const Parallel::ParaWorld& pw = coll.find(Parallel::ParaTag::pw);
    EXPECT_TRUE(pw.valid());
    EXPECT_EQ(pw.size(), 1);
}

TEST(ParaCollectionMpiTest, FindMissingReturnsInvalid)
{
    Parallel::ParaCollection coll;
    coll.add(std::unique_ptr<Parallel::ParaWorld>(
        new Parallel::ParaKmeshWorld(MPI_COMM_WORLD, 1, 0, 1, 4, 1)));

    const Parallel::ParaWorld& missing = coll.find("nonexistent");
    EXPECT_FALSE(missing.valid());
}

TEST(ParaCollectionMpiTest, FindAsSubclass)
{
    Parallel::ParaCollection coll;
    coll.add(std::unique_ptr<Parallel::ParaWorld>(
        new Parallel::ParaKmeshWorld(MPI_COMM_WORLD, 1, 0, 1, 8, 1)));

    const Parallel::ParaKmeshWorld* kmesh = coll.find_as<Parallel::ParaKmeshWorld>(Parallel::ParaTag::kmesh);
    ASSERT_NE(kmesh, nullptr);
    EXPECT_EQ(kmesh->nkstot(), 8);
    EXPECT_EQ(kmesh->nks_local(), 8);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}

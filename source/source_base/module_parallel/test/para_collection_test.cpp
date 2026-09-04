#include "gtest/gtest.h"

#include "../para_collection.h"
#include "../para_kmesh_world.h"
#include "../para_tag.h"

TEST(ParaCollectionTest, DefaultIsEmpty)
{
    const Parallel::ParaCollection coll;
    EXPECT_EQ(coll.size(), 0u);
}

TEST(ParaCollectionTest, AddAndFind)
{
    Parallel::ParaCollection coll;
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));
    coll.add(std::unique_ptr<Parallel::ParaWorld>(new Parallel::ParaKmeshWorld(4, 1)));
    EXPECT_EQ(coll.size(), 2u);

    const Parallel::ParaWorld& pw = coll.find(Parallel::ParaTag::pw);
    EXPECT_EQ(pw.tag(), "pw");
    EXPECT_TRUE(pw.valid());

    const Parallel::ParaWorld& kmesh = coll.find(Parallel::ParaTag::kmesh);
    EXPECT_EQ(kmesh.tag(), "kmesh");
    EXPECT_TRUE(kmesh.valid());
}

TEST(ParaCollectionTest, FindMissingReturnsEmpty)
{
    Parallel::ParaCollection coll;
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));

    const Parallel::ParaWorld& missing = coll.find("nonexistent");
    EXPECT_TRUE(missing.tag().empty());
    EXPECT_FALSE(missing.valid());
}

TEST(ParaCollectionTest, DuplicateTagRejected)
{
    Parallel::ParaCollection coll;
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));
    EXPECT_EQ(coll.size(), 1u);
}

TEST(ParaCollectionTest, FindAsSubclass)
{
    Parallel::ParaCollection coll;
    coll.add(std::unique_ptr<Parallel::ParaWorld>(new Parallel::ParaKmeshWorld(6, 1)));

    const Parallel::ParaKmeshWorld* kmesh = coll.find_as<Parallel::ParaKmeshWorld>(Parallel::ParaTag::kmesh);
    ASSERT_NE(kmesh, nullptr);
    EXPECT_EQ(kmesh->nkstot(), 6);
    EXPECT_EQ(kmesh->kpar(), 1);
}

TEST(ParaCollectionTest, FindAllEightDomains)
{
    Parallel::ParaCollection coll;
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::pw));
    coll.add(std::unique_ptr<Parallel::ParaWorld>(new Parallel::ParaKmeshWorld(4, 1)));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::bsame_kdiff));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::bdiff_ksame));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::rgrid));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::diag));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::matrix));
    coll.add(Parallel::ParaWorld::make_serial(Parallel::ParaTag::atom));
    EXPECT_EQ(coll.size(), 8u);

    EXPECT_TRUE(coll.find(Parallel::ParaTag::pw).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::kmesh).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::bsame_kdiff).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::bdiff_ksame).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::rgrid).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::diag).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::matrix).valid());
    EXPECT_TRUE(coll.find(Parallel::ParaTag::atom).valid());
}

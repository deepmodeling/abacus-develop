#include "gtest/gtest.h"
#include "source_lcao/module_hcontainer/hcontainer.h"

// Tests for add_value_intersection and add_value_union.
// HContainer is built via Parallel_Orbitals (serial):
//   nat=2 atoms, atom0: 2 orbitals, atom1: 3 orbitals, nlocal=5
//   All (i,j) pairs at R=(0,0,0) are pre-inserted and zero-allocated.
//   Values are written directly via find_matrix/get_atom_pair.
//
// For extra R vectors, register via get_atom_pair(i,j).get_HR_values(rx,ry,rz),
// then call hc->allocate(nullptr, true) once to reallocate (only in that case).

class AddValueTest : public ::testing::Test
{
  protected:
    Parallel_Orbitals paraV;
    int iat2iwt[2] = {0, 2};

    void SetUp() override
    {
        paraV.set_serial(5, 5);
        paraV.set_atomic_trace(iat2iwt, 2, 5);
    }

    // Insert all 2*2 pairs at R=(0,0,0) and allocate memory (zero-initialized).
    void insert_all_pairs(hamilt::HContainer<double>* hc)
    {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
            {
                hamilt::AtomPair<double> ap(i, j, 0, 0, 0, &paraV);
                hc->insert_pair(ap);
            }
        hc->allocate(nullptr, true);
    }

    // Build an HContainer with only R=(0,0,0), writing the given values into each pair.
    // fill: { {i, j, values}, ... }  where values.size() == nw(i) * nw(j)
    hamilt::HContainer<double>* make_hc(
        const std::vector<std::tuple<int, int, std::vector<double>>>& fill)
    {
        auto* hc = new hamilt::HContainer<double>(&paraV);
        insert_all_pairs(hc);
        for (auto& [i, j, vals] : fill)
        {
            double* ptr = hc->find_matrix(i, j, 0, 0, 0)->get_pointer();
            for (int k = 0; k < (int)vals.size(); k++)
                ptr[k] = vals[k];
        }
        return hc;
    }

    // Build an HContainer that also has R=(rx,ry,rz) for pair (i,j), with values written.
    // Used only for multi-R tests; calls allocate a second time to include the extra R.
    hamilt::HContainer<double>* make_hc_multiR(
        int i, int j, int rx, int ry, int rz,
        const std::vector<double>& vals_000,
        const std::vector<double>& vals_R)
    {
        auto* hc = new hamilt::HContainer<double>(&paraV);
        insert_all_pairs(hc);
        // Register extra R vector
        hc->get_atom_pair(i, j).get_HR_values(rx, ry, rz);
        // Reallocate so the extra R is included in the wrapper
        hc->allocate(nullptr, true);
        // Write R=(0,0,0) values
        double* p0 = hc->find_matrix(i, j, 0, 0, 0)->get_pointer();
        for (int k = 0; k < (int)vals_000.size(); k++)
            p0[k] = vals_000[k];
        // Write extra-R values
        double* pR = hc->find_matrix(i, j, rx, ry, rz)->get_pointer();
        for (int k = 0; k < (int)vals_R.size(); k++)
            pR[k] = vals_R[k];
        return hc;
    }
};

// ═══════════════════════════════════════════════════════════════════
// add_value_intersection tests
// ═══════════════════════════════════════════════════════════════════

// 1. Identical sparsity pattern: this += 1.0 * other, each element is the sum of both
TEST_F(AddValueTest, intersection_same_sparsity)
{
    // pair(0,1): 2×3=6 elements; pair(1,0): 3×2=6 elements
    auto* dst = make_hc({
        {0, 1, {1, 2, 3, 4, 5, 6}},
        {1, 0, {7, 8, 9, 10, 11, 12}},
    });
    auto* src = make_hc({
        {0, 1, {10, 20, 30, 40, 50, 60}},
        {1, 0, {70, 80, 90, 100, 110, 120}},
    });

    dst->add_value_intersection(*src, 1.0);

    double* p01 = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(p01[0], 11.0);
    EXPECT_DOUBLE_EQ(p01[1], 22.0);
    EXPECT_DOUBLE_EQ(p01[5], 66.0);

    double* p10 = dst->find_matrix(1, 0, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(p10[0], 77.0);
    EXPECT_DOUBLE_EQ(p10[5], 132.0);

    delete dst;
    delete src;
}

// 2. factor parameter: this += 2.0 * other
TEST_F(AddValueTest, intersection_factor)
{
    auto* dst = make_hc({{0, 1, {1, 0, 0, 0, 0, 0}}});
    auto* src = make_hc({{0, 1, {3, 0, 0, 0, 0, 0}}});

    dst->add_value_intersection(*src, 2.0);

    double* ptr = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 7.0);    // 1 + 2*3 = 7
    EXPECT_DOUBLE_EQ(ptr[1], 0.0);

    delete dst;
    delete src;
}

// 3. Negative factor: this += -1.0 * other
TEST_F(AddValueTest, intersection_negative_factor)
{
    auto* dst = make_hc({{0, 1, {10, 20, 30, 40, 50, 60}}});
    auto* src = make_hc({{0, 1, {1, 2, 3, 4, 5, 6}}});

    dst->add_value_intersection(*src, -1.0);

    double* ptr = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 9.0);
    EXPECT_DOUBLE_EQ(ptr[1], 18.0);
    EXPECT_DOUBLE_EQ(ptr[5], 54.0);

    delete dst;
    delete src;
}

// 4. other has pairs absent from dst: only the intersection is added;
//    values in other for pairs not in dst do not affect dst
TEST_F(AddValueTest, intersection_partial_overlap)
{
    // dst only fills (0,1); (1,0) stays 0 (zero-allocated by make_hc)
    auto* dst = make_hc({{0, 1, {1, 2, 3, 4, 5, 6}}});
    // src fills both (0,1) and (1,0)
    auto* src = make_hc({
        {0, 1, {10, 20, 30, 40, 50, 60}},
        {1, 0, {99, 99, 99, 99, 99, 99}},
    });

    dst->add_value_intersection(*src, 1.0);

    // (0,1) should be correctly accumulated
    double* p01 = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(p01[0], 11.0);
    EXPECT_DOUBLE_EQ(p01[5], 66.0);

    // (1,0) exists in dst (zero-initialized); intersection iterates this's pairs,
    // so it finds (1,0) in dst and adds src's (1,0) values: 0+99=99
    double* p10 = dst->find_matrix(1, 0, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(p10[0], 99.0);

    delete dst;
    delete src;
}

// 5. Multiple R vectors: only matching (i,j,R) triples are added; R present in dst but not src is unchanged
TEST_F(AddValueTest, intersection_multi_R)
{
    // dst: (0,1) has R=(0,0,0) and R=(1,0,0)
    auto* dst = make_hc_multiR(0, 1, 1, 0, 0,
                                {1, 0, 0, 0, 0, 0},
                                {2, 0, 0, 0, 0, 0});
    // src: (0,1) has only R=(0,0,0)
    auto* src = make_hc({{0, 1, {10, 0, 0, 0, 0, 0}}});

    dst->add_value_intersection(*src, 1.0);

    // R=(0,0,0): accumulated
    EXPECT_DOUBLE_EQ(dst->find_matrix(0, 1, 0, 0, 0)->get_pointer()[0], 11.0);
    // R=(1,0,0): absent in src, so dst value is unchanged
    EXPECT_DOUBLE_EQ(dst->find_matrix(0, 1, 1, 0, 0)->get_pointer()[0], 2.0);

    delete dst;
    delete src;
}

// ═══════════════════════════════════════════════════════════════════
// add_value_union tests
// ═══════════════════════════════════════════════════════════════════

// 6. Basic correctness: same sparsity as intersection case should give the same result
//    (old code produced all zeros here; correct after bug fix)
TEST_F(AddValueTest, union_basic_sum)
{
    auto* dst = make_hc({{0, 1, {1, 2, 3, 4, 5, 6}}});
    auto* src = make_hc({{0, 1, {10, 20, 30, 40, 50, 60}}});

    dst->add_value_union(*src, 1.0);

    double* ptr = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 11.0);
    EXPECT_DOUBLE_EQ(ptr[1], 22.0);
    EXPECT_DOUBLE_EQ(ptr[5], 66.0);

    delete dst;
    delete src;
}

// 7. Core bug regression: other's data must not be corrupted after the call
//    (old add_value_union zeroed out other; this test specifically checks the fix)
TEST_F(AddValueTest, union_does_not_corrupt_other)
{
    auto* dst = make_hc({});  // all pairs default to 0
    auto* src = make_hc({{0, 1, {1, 2, 3, 4, 5, 6}}});

    dst->add_value_union(*src, 1.0);

    // src data must be fully preserved
    double* src_ptr = src->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(src_ptr[0], 1.0);
    EXPECT_DOUBLE_EQ(src_ptr[1], 2.0);
    EXPECT_DOUBLE_EQ(src_ptr[5], 6.0);

    // dst should be correctly accumulated
    double* dst_ptr = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(dst_ptr[0], 1.0);
    EXPECT_DOUBLE_EQ(dst_ptr[5], 6.0);

    delete dst;
    delete src;
}

// 8. factor parameter: this += 2.0 * other
TEST_F(AddValueTest, union_factor)
{
    auto* dst = make_hc({{0, 1, {1, 0, 0, 0, 0, 0}}});
    auto* src = make_hc({{0, 1, {3, 0, 0, 0, 0, 0}}});

    dst->add_value_union(*src, 2.0);

    double* ptr = dst->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 7.0);   // 1 + 2*3 = 7

    delete dst;
    delete src;
}

// 9. Three successive union accumulations: simulates accumulating dH terms as in write_dH_sum;
//    also verifies that each term's data is not corrupted after each call
TEST_F(AddValueTest, union_accumulate_three_terms)
{
    auto* sum = make_hc({});

    auto* t1 = make_hc({{0, 1, {1, 0, 0, 0, 0, 0}}});
    auto* t2 = make_hc({{0, 1, {0, 2, 0, 0, 0, 0}}});
    auto* t3 = make_hc({{0, 1, {0, 0, 3, 0, 0, 0}}});

    sum->add_value_union(*t1, 1.0);
    sum->add_value_union(*t2, 1.0);
    sum->add_value_union(*t3, 1.0);

    double* ptr = sum->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 1.0);
    EXPECT_DOUBLE_EQ(ptr[1], 2.0);
    EXPECT_DOUBLE_EQ(ptr[2], 3.0);

    // Each term's data must remain intact
    EXPECT_DOUBLE_EQ(t1->find_matrix(0, 1, 0, 0, 0)->get_pointer()[0], 1.0);
    EXPECT_DOUBLE_EQ(t2->find_matrix(0, 1, 0, 0, 0)->get_pointer()[1], 2.0);
    EXPECT_DOUBLE_EQ(t3->find_matrix(0, 1, 0, 0, 0)->get_pointer()[2], 3.0);

    delete sum;
    delete t1;
    delete t2;
    delete t3;
}

// 10. union introduces a new R vector: other has an R absent from dst, which should be inserted and assigned
TEST_F(AddValueTest, union_new_R_from_other)
{
    // dst: (0,1) has only R=(0,0,0), value 1
    auto* dst = make_hc({{0, 1, {1, 0, 0, 0, 0, 0}}});
    // src: (0,1) has R=(0,0,0) and R=(1,0,0)
    auto* src = make_hc_multiR(0, 1, 1, 0, 0,
                                {10, 0, 0, 0, 0, 0},
                                {99, 0, 0, 0, 0, 0});

    dst->add_value_union(*src, 1.0);

    // R=(0,0,0): correctly accumulated
    EXPECT_DOUBLE_EQ(dst->find_matrix(0, 1, 0, 0, 0)->get_pointer()[0], 11.0);

    // R=(1,0,0): newly inserted and assigned (0 + 99 = 99)
    double* pR = dst->find_matrix(0, 1, 1, 0, 0)->get_pointer();
    EXPECT_NE(pR, nullptr);
    EXPECT_DOUBLE_EQ(pR[0], 99.0);

    // src data must not be corrupted
    EXPECT_DOUBLE_EQ(src->find_matrix(0, 1, 0, 0, 0)->get_pointer()[0], 10.0);
    EXPECT_DOUBLE_EQ(src->find_matrix(0, 1, 1, 0, 0)->get_pointer()[0], 99.0);

    delete dst;
    delete src;
}

// 11. Multiple unions followed by intersection: verifies the two operations compose correctly
TEST_F(AddValueTest, union_then_intersection)
{
    auto* acc = make_hc({});
    auto* t1  = make_hc({{0, 1, {3, 0, 0, 0, 0, 0}}});
    auto* t2  = make_hc({{0, 1, {0, 5, 0, 0, 0, 0}}});
    auto* ref = make_hc({{0, 1, {2, 2, 0, 0, 0, 0}}});

    acc->add_value_union(*t1, 1.0);
    acc->add_value_union(*t2, 1.0);

    // acc = {3, 5, 0, ...}; subtract ref={2,2,...} via intersection
    acc->add_value_intersection(*ref, -1.0);

    double* ptr = acc->find_matrix(0, 1, 0, 0, 0)->get_pointer();
    EXPECT_DOUBLE_EQ(ptr[0], 1.0);   // 3 - 2
    EXPECT_DOUBLE_EQ(ptr[1], 3.0);   // 5 - 2

    delete acc;
    delete t1;
    delete t2;
    delete ref;
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

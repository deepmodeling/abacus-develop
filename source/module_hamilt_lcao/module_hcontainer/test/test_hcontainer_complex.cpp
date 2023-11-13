#include "gtest/gtest.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"

/**
 * Unit test of HContainer
 * HContainer is a container of hamilt::AtomPair, in this test, we test the following functions:
 * 1. insert_pair
 * 2. find_pair
 * 3. get_atom_pair
 * 4. fix_R and unfix_R
 * 5. fix_gamma
 * 6. loop_R
 * 7. size_atom_pairs
 * 8. data
 *
 */

// Unit test of HContainer with gtest framework
class HContainerTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        // set up a unitcell, with one element and three atoms, each atom has 2 orbitals
        ucell.ntype = 1;
        ucell.nat = 3;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it = new int[ucell.nat];
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
        }
        ucell.atoms[0].nw = 2;

        // set up a HContainer with ucell
        HR = new hamilt::HContainer<std::complex<double>>(ucell);
    }

    void TearDown() override
    {
        delete HR;
        delete[] ucell.atoms;
        delete[] ucell.iat2it;
    }

    UnitCell ucell;
    hamilt::HContainer<std::complex<double>>* HR;
};

// using TEST_F to test HContainer::insert_pair
TEST_F(HContainerTest, insert_pair)
{
    // check HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_col_size(), 2);
    // set up a hamilt::AtomPair
    hamilt::AtomPair<std::complex<double>> atom_ij(0, 3);
    atom_ij.set_size(2, 2);
    atom_ij.allocate(nullptr, true);
    // insert atom_ij into HR
    HR->insert_pair(atom_ij);
    // check if atom_ij is inserted into HR
    EXPECT_EQ(HR->size_atom_pairs(), 10);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_col_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_atom_j(), 3);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_col_size(), 2);
    // insert atom_ij again
    HR->insert_pair(atom_ij);
    // check if atom_ij is inserted into HR
    EXPECT_EQ(HR->size_atom_pairs(), 10);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_atom_j(), 3);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0, 3).get_col_size(), 2);
    // set up another hamilt::AtomPair
    hamilt::AtomPair<std::complex<double>> atom_kl(1, 0);
    atom_kl.set_size(2, 2);
    atom_kl.allocate(nullptr, true);
    std::complex<double> tmp_array[4] = {1, 2, 3, 4};
    atom_kl.get_HR_values(1, 0, 0).add_array(&tmp_array[0]);
    // insert atom_kl into HR
    HR->insert_pair(atom_kl);
    // check if atom_kl is inserted into HR
    EXPECT_EQ(HR->size_atom_pairs(), 10);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_i(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_atom_j(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(2, 2).get_col_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(1, 0).get_atom_i(), 1);
    EXPECT_EQ(HR->get_atom_pair(1, 0).get_atom_j(), 0);
    EXPECT_EQ(HR->get_atom_pair(1, 0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(1, 0).get_col_size(), 2);
    // check if data is correct
    std::complex<double>* data_ptr = HR->get_atom_pair(1, 0).get_HR_values(1, 0, 0).get_pointer();
    EXPECT_EQ(data_ptr[0], std::complex<double>(1));
    EXPECT_EQ(data_ptr[1], std::complex<double>(2));
    EXPECT_EQ(data_ptr[2], std::complex<double>(3));
    EXPECT_EQ(data_ptr[3], std::complex<double>(4));

    // test copy constructor for HContainer
    hamilt::HContainer<std::complex<double>> HR_copy(*HR);
    // test move constructor for HContainer
    hamilt::HContainer<std::complex<double>> HR_move(std::move(HR_copy));
}

// using TEST_F to test HContainer::find_pair
TEST_F(HContainerTest, find_pair)
{
    // check HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_j(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0).get_col_size(), 2);
    // find atom_ij in HR
    hamilt::AtomPair<std::complex<double>>* atom_ij_ptr = HR->find_pair(0, 1);
    // check if atom_ij is found
    EXPECT_EQ(atom_ij_ptr->get_atom_i(), 0);
    EXPECT_EQ(atom_ij_ptr->get_atom_j(), 1);
    EXPECT_EQ(atom_ij_ptr->get_row_size(), 2);
    EXPECT_EQ(atom_ij_ptr->get_col_size(), 2);
    // find atom_kl in HR
    hamilt::AtomPair<std::complex<double>>* atom_kl_ptr = HR->find_pair(1, 0);
    // check if atom_kl is found
    EXPECT_EQ(atom_kl_ptr->get_atom_i(), 1);
    EXPECT_EQ(atom_kl_ptr->get_atom_j(), 0);
    EXPECT_EQ(atom_kl_ptr->get_row_size(), 2);
    EXPECT_EQ(atom_kl_ptr->get_col_size(), 2);
    // find atom_ij not in HR
    hamilt::AtomPair<std::complex<double>>* atom_ij_ptr2 = HR->find_pair(0, 3);
    // check if atom_ij is found
    EXPECT_EQ(atom_ij_ptr2, nullptr);
}

// using TEST_F to test HContainer::get_atom_pair, both with atom_i, atom_j and with index
TEST_F(HContainerTest, get_atom_pair)
{
    // check  HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_j(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0).get_col_size(), 2);
    // get atom_ij in HR
    hamilt::AtomPair<std::complex<double>>& atom_ij_ref = HR->get_atom_pair(0, 1);
    // check if atom_ij is found
    EXPECT_EQ(atom_ij_ref.get_atom_i(), 0);
    EXPECT_EQ(atom_ij_ref.get_atom_j(), 1);
    EXPECT_EQ(atom_ij_ref.get_row_size(), 2);
    EXPECT_EQ(atom_ij_ref.get_col_size(), 2);
    // get atom_kl in HR
    hamilt::AtomPair<std::complex<double>>& atom_kl_ref = HR->get_atom_pair(1, 0);
    // check if atom_kl is found
    EXPECT_EQ(atom_kl_ref.get_atom_i(), 1);
    EXPECT_EQ(atom_kl_ref.get_atom_j(), 0);
    EXPECT_EQ(atom_kl_ref.get_row_size(), 2);
    EXPECT_EQ(atom_kl_ref.get_col_size(), 2);
    // get atom_ij in HR with index
    hamilt::AtomPair<std::complex<double>>& atom_ij_ref2 = HR->get_atom_pair(0);
    // check if atom_ij is found
    EXPECT_EQ(atom_ij_ref2.get_atom_i(), 0);
    EXPECT_EQ(atom_ij_ref2.get_atom_j(), 0);
    EXPECT_EQ(atom_ij_ref2.get_row_size(), 2);
    EXPECT_EQ(atom_ij_ref2.get_col_size(), 2);
    // get atom_kl in HR with index
    hamilt::AtomPair<std::complex<double>>& atom_kl_ref2 = HR->get_atom_pair(8);
    // check if atom_kl is found
    EXPECT_EQ(atom_kl_ref2.get_atom_i(), 2);
    EXPECT_EQ(atom_kl_ref2.get_atom_j(), 2);
    EXPECT_EQ(atom_kl_ref2.get_row_size(), 2);
    EXPECT_EQ(atom_kl_ref2.get_col_size(), 2);
}

// using TEST_F to test HContainer::fix_R and unfix_R
TEST_F(HContainerTest, fix_R)
{
    // check HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_j(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0).get_col_size(), 2);
    EXPECT_EQ(HR->size_R_loop(), 1);
    // fix R
    EXPECT_EQ(HR->fix_R(0, 0, 0), true);
    // check if R is fixed
    EXPECT_EQ(HR->get_current_R(), 0);
    // fix R again
    EXPECT_EQ(HR->fix_R(0, 0, 0), true);
    // check if R is fixed
    EXPECT_EQ(HR->get_current_R(), 0);
    // fix another R
    EXPECT_EQ(HR->fix_R(1, 0, 0), false);
    // check if R is fixed
    EXPECT_EQ(HR->get_current_R(), -1);
    // unfix R
    HR->unfix_R();
    // check if R is unfixed
    EXPECT_EQ(HR->get_current_R(), -1);
}

// using TEST_F to test HContainer::fix_gamma
TEST_F(HContainerTest, fix_gamma)
{
    // check HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_j(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0).get_col_size(), 2);
    EXPECT_EQ(HR->size_R_loop(), 1);
    hamilt::AtomPair<std::complex<double>> atom_ij(0, 1);
    atom_ij.set_size(2, 2);
    hamilt::BaseMatrix<std::complex<double>>& tmp = atom_ij.get_HR_values(1, 0, 0);
    std::complex<double> tmp_array[4] = {1, 2, 3, 4};
    tmp.add_array(tmp_array);
    // insert atom_ij into HR
    HR->insert_pair(atom_ij);
    EXPECT_EQ(HR->size_R_loop(), 2);
    // fix gamma
    EXPECT_EQ(HR->is_gamma_only(), false);
    HR->fix_gamma();
    // check if gamma is fixed
    EXPECT_EQ(HR->is_gamma_only(), true);
    EXPECT_EQ(HR->size_R_loop(), 1);
    // fix gamma again
    HR->fix_gamma();
    // check if gamma is fixed
    EXPECT_EQ(HR->is_gamma_only(), true);
}

/**
 * using TEST_F to test HContainer::loop_R,
 * step: 1. size_R_loop(), 2. for-loop, loop_R(), 3. fix_R(), 4. do something
 */
TEST_F(HContainerTest, loop_R)
{
    // 1. size_R_loop()
    int size_for_loop_R = HR->size_R_loop();
    EXPECT_EQ(size_for_loop_R, 1);
    // 2. for-loop, loop_R()
    int rx, ry, rz;
    for (int i = 0; i < size_for_loop_R; i++)
    {
        HR->loop_R(i, rx, ry, rz);
        EXPECT_EQ(rx, 0);
        EXPECT_EQ(ry, 0);
        EXPECT_EQ(rz, 0);
        HR->fix_R(rx, ry, rz);
        // check if R is fixed
        EXPECT_EQ(HR->get_current_R(), i);
        // 4. do something
    }
    HR->unfix_R();
    // check if R is unfixed
    EXPECT_EQ(HR->get_current_R(), -1);
}

// using TEST_F to test HContainer::size_atom_pairs
// 1. test with R fixed
// 2. test with R unfixed
TEST_F(HContainerTest, size_atom_pairs)
{
    // get size_R_loop
    int size_R_loop = HR->size_R_loop();
    EXPECT_EQ(size_R_loop, 1);
    // 1. test with R fixed
    // fix R
    EXPECT_EQ(HR->get_current_R(), -1);
    bool ok = HR->fix_R(0, 0, 0);
    // check if R is fixed
    EXPECT_EQ(ok, true);
    EXPECT_EQ(HR->get_current_R(), 0);
    // get AP size
    int AP_size = HR->size_atom_pairs();
    EXPECT_EQ(AP_size, 9);
    // fix to another R
    hamilt::AtomPair<std::complex<double>> atom_ij(0, 1);
    atom_ij.set_size(2, 2);
    atom_ij.allocate(nullptr, true);
    hamilt::BaseMatrix<std::complex<double>>& tmp = atom_ij.get_HR_values(1, 0, 0);
    std::complex<double> tmp_array[4] = {1, 2, 3, 4};
    tmp.add_array(tmp_array);
    // insert atom_ij into HR
    HR->insert_pair(atom_ij);
    // get size_R_loop again, it should be 2
    size_R_loop = HR->size_R_loop();
    // if R is fixed, size_R_loop will not be reset, so it should be 1
    EXPECT_EQ(size_R_loop, 1);
    // unfix R
    HR->unfix_R();
    // check size_R_loop after R index unfixed
    size_R_loop = HR->size_R_loop();
    EXPECT_EQ(size_R_loop, 2);
    ok = HR->fix_R(1, 0, 0);
    // check if R is fixed
    EXPECT_EQ(ok, true);
    EXPECT_EQ(HR->get_current_R(), 1);
    EXPECT_EQ(HR->size_atom_pairs(), 1);
    // check if tmp_atom_pairs is correct
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_atom_j(), 1);
    EXPECT_EQ(HR->get_atom_pair(0).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0).get_col_size(), 2);
    int* R_ptr = HR->get_atom_pair(0).get_R_index();
    EXPECT_EQ(R_ptr[0], 1);
    EXPECT_EQ(R_ptr[1], 0);
    EXPECT_EQ(R_ptr[2], 0);
    EXPECT_EQ(HR->get_atom_pair(0).get_R_index(5), nullptr);
    // check if data is correct
    std::complex<double>* data_ptr = HR->get_atom_pair(0).get_pointer();
    EXPECT_EQ(data_ptr[0], std::complex<double>(1));
    EXPECT_EQ(data_ptr[1], std::complex<double>(2));
    EXPECT_EQ(data_ptr[2], std::complex<double>(3));
    EXPECT_EQ(data_ptr[3], std::complex<double>(4));
    // 2. test with R unfixed
    // unfix R
    HR->unfix_R();
    // check if R is unfixed
    EXPECT_EQ(HR->get_current_R(), -1);
    // get AP size
    AP_size = HR->size_atom_pairs();
    EXPECT_EQ(AP_size, 9);
    // fix to another R with no AP
    ok = HR->fix_R(2, 0, 0);
    // check if R is fixed
    EXPECT_EQ(ok, false);
    EXPECT_EQ(HR->get_current_R(), -1);
    EXPECT_EQ(HR->size_atom_pairs(), 9);
}

// using TEST_F to test HContainer::data()
TEST_F(HContainerTest, data)
{
    // set up a hamilt::AtomPair
    hamilt::AtomPair<std::complex<double>> atom_ij(0, 1);
    atom_ij.set_size(2, 2);
    atom_ij.allocate(nullptr, true);
    hamilt::BaseMatrix<std::complex<double>>& tmp = atom_ij.get_HR_values(0, 0, 0);
    std::complex<double> tmp_array[4] = {1, 2, 3, 4};
    tmp.add_array(tmp_array);
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    // insert atom_ij into HR
    HR->insert_pair(atom_ij);
    // check if atom_ij is inserted into HR
    EXPECT_EQ(HR->size_atom_pairs(), 9);
    EXPECT_EQ(HR->get_atom_pair(0, 1).get_atom_i(), 0);
    EXPECT_EQ(HR->get_atom_pair(0, 1).get_atom_j(), 1);
    EXPECT_EQ(HR->get_atom_pair(0, 1).get_row_size(), 2);
    EXPECT_EQ(HR->get_atom_pair(0, 1).get_col_size(), 2);
    // get data pointer
    std::complex<double>* data_ptr = HR->data(0, 1);
    // check if data pointer is correct
    EXPECT_EQ(data_ptr, HR->get_atom_pair(0, 1).get_pointer());
    EXPECT_EQ(data_ptr, HR->get_atom_pair(0, 1).get_pointer(0));
    int r_index[3] = {0, 0, 0};
    EXPECT_EQ(data_ptr, HR->data(0, 1, r_index));
    EXPECT_EQ(HR->data(0, 10), nullptr);
    EXPECT_EQ(HR->data(0, 10, r_index), nullptr);
    // check if data is correct
    EXPECT_EQ(data_ptr[0], std::complex<double>(1));
    EXPECT_EQ(data_ptr[1], std::complex<double>(2));
    EXPECT_EQ(data_ptr[2], std::complex<double>(3));
    EXPECT_EQ(data_ptr[3], std::complex<double>(4));
}

// using TEST_F to test functions in BaseMatrix
// 1. test constructor with existed data
// 4. test add_element
// 5. test get_value
TEST_F(HContainerTest, basematrix_funcs)
{
    // 1. test constructor with existed data
    std::complex<double> data_ptr[4] = {1, 2, 3, 4};
    hamilt::BaseMatrix<std::complex<double>> BM(2, 2, &data_ptr[0]);
    // check if data is correct
    EXPECT_EQ(BM.get_value(0, 0), std::complex<double>(1));
    EXPECT_EQ(BM.get_value(0, 1), std::complex<double>(2));
    EXPECT_EQ(BM.get_value(1, 0), std::complex<double>(3));
    EXPECT_EQ(BM.get_value(1, 1), std::complex<double>(4));
    // copy BM to check copy constructor
    hamilt::BaseMatrix<std::complex<double>> BM_copy(BM);
    BM_copy = BM;
    BM_copy = hamilt::BaseMatrix<std::complex<double>>(BM);
    // check if data is correct
    EXPECT_EQ(BM_copy.get_value(0, 0), std::complex<double>(1));
    EXPECT_EQ(BM_copy.get_value(0, 1), std::complex<double>(2));
    EXPECT_EQ(BM_copy.get_value(1, 0), std::complex<double>(3));
    EXPECT_EQ(BM_copy.get_value(1, 1), std::complex<double>(4));
} 

// using TEST_F to test functions in AtomPair
// 1. constructor
// 2. copy assignment
// 3. move assignment
// 4. identify
// 5. add_to_matrix
// 6. add_to_array
// 7. get_matrix_value
// 8. get_R_index with out of range
// 9. get_value
// 10. get_value_size
TEST_F(HContainerTest, atompair_funcs)
{
    // 1. constructor
    Parallel_Orbitals PO;
    PO.atom_begin_row.resize(3); // natom = 2, size should be natom + 1
    PO.atom_begin_col.resize(3);
    for(int i=0;i<3;i++)
    {
        PO.atom_begin_row[i] = i*2; // nw = 2, value should be i*nw
        PO.atom_begin_col[i] = i*2;
    }
    PO.nrow = 4;
    PO.ncol = 4;
    hamilt::AtomPair<std::complex<double>> atom_ij(0, 0, &PO, nullptr);
    hamilt::AtomPair<std::complex<double>> atom_ij2(0, 1, 1, 1, 1, &PO, nullptr);
    hamilt::AtomPair<std::complex<double>> atom_ij3(1, 0, PO.atom_begin_row.data(), PO.atom_begin_col.data(), 2, nullptr);
    hamilt::AtomPair<std::complex<double>> atom_ij33(1, 1, 1, 1, 1, PO.atom_begin_row.data(), PO.atom_begin_col.data(), 2, nullptr);
    EXPECT_EQ(atom_ij<atom_ij2, true);
    EXPECT_EQ(atom_ij2<atom_ij, false);
    EXPECT_EQ(atom_ij3<atom_ij2, false);
    EXPECT_EQ(atom_ij2<atom_ij3, true);
    atom_ij3 = atom_ij;
    atom_ij33 = hamilt::AtomPair<std::complex<double>>(atom_ij2);
    EXPECT_EQ(atom_ij.identify(atom_ij3), true);
    EXPECT_EQ(atom_ij.identify(atom_ij2), false);
    EXPECT_EQ(atom_ij2.identify(atom_ij33.get_atom_i(), atom_ij33.get_atom_j()), true);
    EXPECT_EQ(atom_ij2.identify(atom_ij3.get_atom_i(), atom_ij3.get_atom_j()), false);
    // 5. add_to_matrix
    // row major case
    std::complex<double> hk_data2[16] = {1, 2, 3, 4,
                                         5, 6, 7, 8,
                                         9, 10, 11, 12,
                                         13, 14, 15, 16};
    // colomn major case
    std::complex<double> hk_data3[16] = {1, 5, 9, 13,
                                         2, 6, 10, 14,
                                         3, 7, 11, 15,
                                         4, 8, 12, 16};
    hamilt::HContainer<std::complex<double>> HR(2);
    for(int atom_i = 0;atom_i<2;++atom_i)
    {
        for(int atom_j = 0; atom_j<2; ++atom_j)
        {
            hamilt::AtomPair<std::complex<double>> tmp(atom_i, atom_j, 0, 0, 0, PO.atom_begin_row.data(), PO.atom_begin_col.data(), 2, nullptr);
            tmp.allocate(nullptr, false);
            std::complex<double>* tmp_data = tmp.get_HR_values(0, 0, 0).get_pointer();
            for(int i=0;i<4;++i)
            {
                tmp_data[i] = atom_i*2 + atom_j*4 + i + 1;
            }
            HR.insert_pair(tmp);
        }
    }
    for(int ir=0;ir<HR.size_R_loop();++ir)
    {
        int rx, ry, rz;
        HR.loop_R(ir, rx, ry, rz);
        HR.fix_R(rx, ry, rz);
        for(int iap = 0;iap<HR.size_atom_pairs();++iap)
        {
            hamilt::AtomPair<std::complex<double>>& tmp = HR.get_atom_pair(iap);
            // row major case
            tmp.add_to_matrix(&hk_data2[0], 4, std::complex<double>(1.0, 0.5), 0);
            // colomn major case
            tmp.add_to_matrix(&hk_data3[0], 4, std::complex<double>(1.0, 0.5), 1);
        }
    }
    HR.unfix_R();
    // check hk_data and hk_data2 are correct
    std::complex<double> hk_data2_correct[16] = {std::complex<double>(2,0.5), std::complex<double>(4,1), std::complex<double>(8,2.5), std::complex<double>(10,3),
                                                 std::complex<double>(8,1.5), std::complex<double>(10,2), std::complex<double>(14,3.5), std::complex<double>(16,4),
                                                 std::complex<double>(12,1.5), std::complex<double>(14,2), std::complex<double>(18,3.5), std::complex<double>(20,4),
                                                 std::complex<double>(18,2.5), std::complex<double>(20,3), std::complex<double>(24,4.5), std::complex<double>(26,5)};
    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            EXPECT_EQ(hk_data2[i*4+j], hk_data2_correct[i*4+j]);
            EXPECT_EQ(hk_data3[j*4+i], hk_data2_correct[i*4+j]);
        }
    }

    // 6. add_to_array
    std::vector<std::complex<double>> hr_array(16, 0.0);
    for(int ir=0;ir<HR.size_R_loop();++ir)
    {
        int rx, ry, rz;
        HR.loop_R(ir, rx, ry, rz);
        HR.fix_R(rx, ry, rz);
        std::complex<double>* ptr1 = hr_array.data();
        for(int iap = 0;iap<HR.size_atom_pairs();++iap)
        {
            auto tmp = HR.get_atom_pair(iap);
            tmp.add_to_array(ptr1, std::complex<double>(1.0, 0.0));
            ptr1 += tmp.get_size();
        }
    }
    HR.unfix_R();
    // check hr_array and hr_array2 are correct
    std::complex<double> correct1;
    std::complex<double> correct_array[16] = {
        1, 2, 3, 4, 
        5, 6, 7, 8, 
        3, 4, 5, 6, 
        7, 8, 9, 10};
    std::complex<double> test_array[16] = {1, 2, 5, 6, 3, 4, 7, 8, 3, 4, 7, 8, 5, 6, 9, 10};

    for(int i=0;i<4;++i)
    {
        for(int j=0;j<4;++j)
        {
            correct1 = correct_array[i*4+j];
            EXPECT_EQ(hr_array[i*4+j], correct1);
        }
    }

    // construct AtomPair from existed matrix
    hamilt::AtomPair<std::complex<double>> atom_ij4(0, 0, &PO, test_array);
    EXPECT_EQ(atom_ij4.get_value(0, 0), correct_array[0]);
    EXPECT_EQ(atom_ij4.get_value(1, 1), correct_array[5]);
    EXPECT_EQ(atom_ij4.get_value(0), correct_array[0]);
    hamilt::AtomPair<std::complex<double>> atom_ij5(0, 1, 1, 1, 1, &PO, &test_array[4]);
    hamilt::AtomPair<std::complex<double>> atom_ij6(1, 0, PO.atom_begin_row.data(), PO.atom_begin_col.data(), 2, &test_array[8]);
    hamilt::AtomPair<std::complex<double>> atom_ij7(1, 1, 1, 1, 1, PO.atom_begin_row.data(), PO.atom_begin_col.data(), 2, &test_array[12]);
    // get_matrix_value will use global2local_row and global2local_col in Parallel_Orbitals
    // so we need to set them
    std::ofstream ofs("test_hcontainer_complex.log");
    PO.set_global2local(4, 4, false, ofs);
    auto checkdata = [&](hamilt::AtomPair<std::complex<double>>& ap_in) {
        auto data_ij4 = ap_in.get_matrix_values();
        int* tmp_index = std::get<0>(data_ij4).data();
        std::complex<double>* tmp_data = std::get<1>(data_ij4);
        double sum_error = 0.0;
        for(int irow = tmp_index[0]; irow < tmp_index[0] + tmp_index[1]; ++irow)
        {
            for(int icol = tmp_index[2]; icol < tmp_index[2] + tmp_index[3]; ++icol)
            {
                sum_error += std::abs((*tmp_data).real() - correct_array[irow*4+icol].real());
                tmp_data++;
            }
        }
        return sum_error;
    };
    EXPECT_EQ(checkdata(atom_ij4), 0.0);
    EXPECT_EQ(checkdata(atom_ij5), 0.0);
    EXPECT_EQ(checkdata(atom_ij6), 0.0);
    EXPECT_EQ(checkdata(atom_ij7), 0.0);
}


int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    return result;
}
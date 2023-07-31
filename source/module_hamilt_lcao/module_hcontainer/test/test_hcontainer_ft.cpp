#include <fstream>

#include "../hcontainer.h"
#include "../output_hcontainer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/complexmatrix.h"
#include "module_base/lapack_connector.h"
#include "module_io/csr_reader.h"
#include "prepare_unitcell.h"

// mock functions
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}
// mocke functions

/************************************************
 *  unit test of read and output hcontainer
 ***********************************************/

/**
 * This unit test read sparse matrices of SR
 * from a file SR.csr, construct SR using
 * hamilt::HContainer<double> and output the
 * matrices of SR to a file SR.out.
 */

class HContainerFTTest : public testing::Test
{
  protected:
    UnitCell* ucell;
    UcellTestPrepare utp = UcellTestLib["Si"];
    // nw is the number of orbitals of each atom
    // it should container ucell.nat elements
    std::vector<int> nw = {13};
    int nlocal;
    void SetUp() override
    {
        ucell = utp.SetUcellInfo(nw, nlocal);
    }
};

TEST_F(HContainerFTTest, ReadAndOutputHContainer)
{
    // read SR
    std::string filename = "./support/data-SR-sparse_SPIN0.csr";
    ModuleIO::csrFileReader<double> csr(filename);
    // std::cout << "csr.getStep " << csr.getStep() << std::endl;
    // std::cout << "csr.getMatrixDimension " << csr.getMatrixDimension() << std::endl;
    // std::cout << "csr.getNumberOfR " << csr.getNumberOfR() << std::endl;
    // std::cout << "nlocal " << nlocal << std::endl;
    //
    // construct paraV
    Parallel_Orbitals paraV;
    std::ofstream ofs("test.log");
    paraV.set_global2local(nlocal, nlocal, false, ofs);
    ofs.close();
    remove("test.log");
    paraV.set_atomic_trace(ucell->iat2iwt.data(), ucell->nat, nlocal);
    // std::cout << paraV.atom_begin_col[0] << " " << paraV.atom_begin_col[1] << std::endl;
    // std::cout << paraV.atom_begin_row[0] << " " << paraV.atom_begin_row[1] << std::endl;
    //
    // construct SR
    hamilt::HContainer<double> SR(&paraV);
    int numberofR = csr.getNumberOfR();
    for (int i = 0; i < numberofR; i++)
    {
        std::vector<int> RCoord = csr.getRCoordinate(i);
        ModuleIO::SparseMatrix<double> sparse_matrix = csr.getMatrix(i);
        for (int iat = 0; iat < ucell->nat; iat++)
        {
            int begin_row = paraV.atom_begin_row[iat];
            int end_row = (iat == ucell->nat - 1) ? paraV.get_row_size() : paraV.atom_begin_row[iat + 1];
            int numberofRow = end_row - begin_row;
            for (int jat = 0; jat < ucell->nat; jat++)
            {
                int begin_col = paraV.atom_begin_col[jat];
                int end_col = (jat == ucell->nat - 1) ? paraV.get_col_size() : paraV.atom_begin_col[jat + 1];
                int numberofCol = end_col - begin_col;
                hamilt::BaseMatrix<double> tmp_matrix(numberofRow, numberofCol);
                int nnz = 0;
                for (const auto& element: sparse_matrix.getElements())
                {
                    int row = element.first.first;
                    int col = element.first.second;
                    if (row < begin_row || row >= end_row || col < begin_col || col >= end_col)
                    {
                        continue;
                    }
                    tmp_matrix.add_element(row - begin_row, col - begin_col, element.second);
                    nnz++;
                }
                if (nnz != 0)
                {
                    auto tmp_ap = hamilt::AtomPair<double>(iat, jat, RCoord[0], RCoord[1], RCoord[2], &paraV);
                    tmp_ap.convert_add(tmp_matrix, RCoord[0], RCoord[1], RCoord[2]);
                    SR.insert_pair(tmp_ap);
                }
            }
        }
    }
    // output SR
    std::ofstream ofs_out("SR.out");
    double sparse_threshold = 1e-10;
    hamilt::Output_HContainer<double> output_SR(&SR, &paraV, *ucell, ofs_out, sparse_threshold, 8);
    // std::cout << SR.size_R_loop() << std::endl;
    // output_SR.write(-2, -1, 0);
    ofs_out << "STEP: " << 0 << std::endl;
    ofs_out << "Matrix Dimension of S(R): " << paraV.get_col_size() << std::endl;
    ofs_out << "Matrix number of S(R): " << SR.size_R_loop() << std::endl;
    output_SR.write();
    ofs_out.close();
    //
    // read in SR.out again
    std::string fileout = "SR.out";
    ModuleIO::csrFileReader<double> csr_out(fileout);
    EXPECT_EQ(csr.getStep(), csr_out.getStep());
    EXPECT_EQ(csr.getMatrixDimension(), csr_out.getMatrixDimension());
    // EXPECT_EQ(csr.getNumberOfR(), csr_out.getNumberOfR());
    //
    int numberofR_out = csr_out.getNumberOfR();
    // compare csr and csr_out
    for (int i = 0; i < numberofR; i++)
    {
        std::vector<int> RCoord = csr.getRCoordinate(i);
        ModuleIO::SparseMatrix<double> sparse_matrix = csr.getMatrix(i);
        for (int j = 0; j < numberofR_out; j++)
        {
            std::vector<int> RCoord_out = csr_out.getRCoordinate(j);
            if (RCoord[0] == RCoord_out[0] && RCoord[1] == RCoord_out[1] && RCoord[2] == RCoord_out[2])
            {
                ModuleIO::SparseMatrix<double> sparse_matrix_out = csr_out.getMatrix(j);
                // EXPECT_EQ(sparse_matrix.getElements().size(), sparse_matrix_out.getElements().size());
                for (const auto& element: sparse_matrix_out.getElements())
                {
                    int row = element.first.first;
                    int col = element.first.second;
                    double value = element.second;
                    auto it = sparse_matrix.getElements().find(std::make_pair(row, col));
                    EXPECT_NE(it, sparse_matrix.getElements().end());
                    EXPECT_NEAR(value, it->second, 1e-10);
                }
            }
        }
    }
    // remove("SR.out");

    // read HR
    filename = "./support/data-HR-sparse_SPIN0.csr";
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
    ModuleIO::csrFileReader<double> hr_csr(filename);
    std::cout << "csr.getStep " << hr_csr.getStep() << std::endl;
    std::cout << "csr.getMatrixDimension " << hr_csr.getMatrixDimension() << std::endl;
    std::cout << "csr.getNumberOfR " << hr_csr.getNumberOfR() << std::endl;
    // construct HR
    hamilt::HContainer<double> HR(&paraV);
    numberofR = hr_csr.getNumberOfR();
    std::cout << "numberofR " << numberofR << std::endl;
    for (int i = 0; i < numberofR; i++)
    {
        std::vector<int> RCoord = hr_csr.getRCoordinate(i);
        //std::cout << i << " RCoord " << RCoord[0] << " " << RCoord[1] << " " << RCoord[2] << std::endl;
        ModuleIO::SparseMatrix<double> sparse_matrix = hr_csr.getMatrix(i);
        //std::cout << "sparse_matrix.getElements().size() " << sparse_matrix.getElements().size() << std::endl;
        for (int iat = 0; iat < ucell->nat; iat++)
        {
            int begin_row = paraV.atom_begin_row[iat];
            int end_row = (iat == ucell->nat - 1) ? paraV.get_row_size() : paraV.atom_begin_row[iat + 1];
            int numberofRow = end_row - begin_row;
            for (int jat = 0; jat < ucell->nat; jat++)
            {
                int begin_col = paraV.atom_begin_col[jat];
                int end_col = (jat == ucell->nat - 1) ? paraV.get_col_size() : paraV.atom_begin_col[jat + 1];
                int numberofCol = end_col - begin_col;
                hamilt::BaseMatrix<double> tmp_matrix(numberofRow, numberofCol);
                int nnz = 0;
                for (const auto& element: sparse_matrix.getElements())
                {
                    int row = element.first.first;
                    int col = element.first.second;
                    //std::cout << "row " << row << " col " << col << " value " << element.second << std::endl;
                    //std::cout << "begin_row " << begin_row << " end_row " << end_row << " begin_col " << begin_col << " end_col " << end_col << std::endl;
                    if (row < begin_row || row >= end_row || col < begin_col || col >= end_col)
                    {
                        continue;
                    }
                    tmp_matrix.add_element(row - begin_row, col - begin_col, element.second);
                    nnz++;
                }
                if (nnz != 0)
                {
                    //std::cout << "nnz " << nnz << std::endl;
                    auto tmp_ap = hamilt::AtomPair<double>(iat, jat, RCoord[0], RCoord[1], RCoord[2], &paraV);
                    tmp_ap.convert_add(tmp_matrix, RCoord[0], RCoord[1], RCoord[2]);
                    HR.insert_pair(tmp_ap);
                }
            }
        }
    }

    // output HR
    std::cout << HR.size_R_loop() << std::endl;

    ModuleBase::ComplexMatrix SGamma(nlocal, nlocal);
    ModuleBase::ComplexMatrix HGamma(nlocal, nlocal);

    // std::complex<double> SGamma[nlocal * nlocal] = {(0.0, 0.0)};
    // std::complex<double> HGamma[nlocal * nlocal] = {(0.0, 0.0)};
    const std::complex<double> kphase = std::complex<double>(1.0, 0.0);
    for (int i = 0; i < HR.size_R_loop(); ++i)
    {
        // get R index
        int rx, ry, rz;
        HR.loop_R(i, rx, ry, rz);
        // only deal with current_R for hR
        HR.fix_R(rx, ry, rz);

        // loop_atom_pairs
        for (int j = 0; j < HR.size_atom_pairs(); ++j)
        {
            // Hk += HR * e^ikR
            HR.get_atom_pair(j).add_to_matrix(HGamma.c, nlocal, kphase, 0);
        }
    }
    /*
    std::cout << "HGamma " << std::endl;
    for (int i = 0; i < nlocal; i++)
    {
        for (int j = i; j < nlocal; j++)
        {
            if (abs(HGamma(i, j)) > 1e-10)
                std::cout << "i = " << i << " j = " << j << HGamma(i, j) << " " << HGamma(j, i) << std::endl;
        }
        std::cout << std::endl;
    }
    */

    for (int i = 0; i < SR.size_R_loop(); ++i)
    {
        // get R index
        int rx, ry, rz;
        SR.loop_R(i, rx, ry, rz);
        // only deal with current_R for hR
        SR.fix_R(rx, ry, rz);

        // loop_atom_pairs
        for (int j = 0; j < SR.size_atom_pairs(); ++j)
        {
            // Hk += HR * e^ikR
            SR.get_atom_pair(j).add_to_matrix(SGamma.c, nlocal, kphase, 0);
        }
    }
    /*
     std::cout << "SGamma " << std::endl;
     for (int i = 0; i < nlocal; i++)
     {
         for (int j = i; j < nlocal; j++)
         {
            if (abs(SGamma(i,j)) > 1e-10)
             std::cout << "i = " << i << " j = " << j << SGamma(i,j) << " " << SGamma(j,i) << std::endl;
         }
         std::cout << std::endl;
     }
     */
}

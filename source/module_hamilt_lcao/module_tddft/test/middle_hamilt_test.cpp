#include "module_hamilt_lcao/module_tddft/middle_hamilt.h"

#include <gtest/gtest.h>
#include <module_base/scalapack_connector.h>
#include <mpi.h>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "tddft_test.h"

/************************************************
 *  unit test of functions in middle_hamilt.h
 ***********************************************/

/**
 * - Tested Function
 *   - half_Hmatrix
 *     - compute H(t+dt/2).
 */

#define doublethreshold 1e-8
Parallel_2D::Parallel_2D(){}
Parallel_2D::~Parallel_2D(){}
Parallel_Orbitals::Parallel_Orbitals()
{
}
Parallel_Orbitals::~Parallel_Orbitals()
{
}

TEST(MiddleHamiltTest, testMiddleHamilt)
{
    std::complex<double>* Htmp;
    std::complex<double>* Hlaststep;
    int nband = 3;
    int nlocal = 4;
    bool print_matrix = false;
    Parallel_Orbitals* pv;
    pv = new Parallel_Orbitals();
    pv->nloc = nlocal * nlocal;
    pv->ncol = nlocal;
    pv->nrow = nlocal;

    // Initialize input matrices
    int info;
    int mb = 1, nb = 1, lda = nlocal, ldc = nlocal;
    int irsrc = 0, icsrc = 0, lld = numroc_(&nlocal, &mb, &myprow, &irsrc, &nprow);
    descinit_(pv->desc, &nlocal, &nlocal, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);

    // Initialize data
    Htmp = new std::complex<double>[nlocal * nlocal];
    Hlaststep = new std::complex<double>[nlocal * nlocal];

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (i == j)
            {
                Htmp[i * nlocal + j] = std::complex<double>(1.0, 0.0);
                Hlaststep[i * nlocal + j] = std::complex<double>(1.0 + 0.2 * (i * nlocal + j), 0.0);
            }
            else
            {
                Hlaststep[i * nlocal + j] = std::complex<double>(0.2 * (i * nlocal + j), 0.0);
            }
        }
    }

    // Call the function
    module_tddft::half_Hmatrix(pv, nband, nlocal, Htmp, Hlaststep, print_matrix);

    // Check the results
    EXPECT_NEAR(Htmp[0].real(), 1.0, doublethreshold);
    EXPECT_NEAR(Htmp[0].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[1].real(), 0.1, doublethreshold);
    EXPECT_NEAR(Htmp[1].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[2].real(), 0.2, doublethreshold);
    EXPECT_NEAR(Htmp[2].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[3].real(), 0.3, doublethreshold);
    EXPECT_NEAR(Htmp[3].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[4].real(), 0.4, doublethreshold);
    EXPECT_NEAR(Htmp[4].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[5].real(), 1.5, doublethreshold);
    EXPECT_NEAR(Htmp[5].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[6].real(), 0.6, doublethreshold);
    EXPECT_NEAR(Htmp[6].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[7].real(), 0.7, doublethreshold);
    EXPECT_NEAR(Htmp[7].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[8].real(), 0.8, doublethreshold);
    EXPECT_NEAR(Htmp[8].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[9].real(), 0.9, doublethreshold);
    EXPECT_NEAR(Htmp[9].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[10].real(), 2.0, doublethreshold);
    EXPECT_NEAR(Htmp[10].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[11].real(), 1.1, doublethreshold);
    EXPECT_NEAR(Htmp[11].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[12].real(), 1.2, doublethreshold);
    EXPECT_NEAR(Htmp[12].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[13].real(), 1.3, doublethreshold);
    EXPECT_NEAR(Htmp[13].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[14].real(), 1.4, doublethreshold);
    EXPECT_NEAR(Htmp[14].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(Htmp[15].real(), 2.5, doublethreshold);
    EXPECT_NEAR(Htmp[15].imag(), 0.0, doublethreshold);

    delete[] Htmp;
    delete[] Hlaststep;
}

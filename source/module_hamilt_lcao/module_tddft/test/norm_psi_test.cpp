#include "module_hamilt_lcao/module_tddft/norm_psi.h"

#include <gtest/gtest.h>
#include <module_base/scalapack_connector.h>
#include <mpi.h>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "tddft_test.h"

/************************************************
 *  unit test of functions in norm_psi.h
 ***********************************************/

/**
 * - Tested Function
 *   - norm_psi
 *     - normalize the wave function.
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

TEST(NormPsiTest, testNormPsi)
{
    std::complex<double>* psi_k;
    std::complex<double>* Stmp;
    int nband = 3;
    int nlocal = 4;
    bool print_matrix = false;
    Parallel_Orbitals* pv;
    pv = new Parallel_Orbitals();
    pv->nloc = nlocal * nlocal;
    pv->nloc_wfc = nlocal * nband;
    pv->ncol = nlocal;
    pv->nrow = nlocal;
    pv->ncol_bands = nband;
    pv->dim0 = 1;
    pv->dim1 = 1;
    pv->nb = 1;

    int dim[2];
    int period[2] = {1, 1};
    int reorder = 0;
    dim[0] = nprow;
    dim[1] = npcol;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &pv->comm_2D);

    // Initialize input matrices
    int info;
    int mb = 1, nb = 1, lda = nband, ldc = nlocal;
    int irsrc = 0, icsrc = 0, lld = numroc_(&nlocal, &mb, &myprow, &irsrc, &nprow),
        lld1 = numroc_(&nband, &mb, &myprow, &irsrc, &nprow);
    descinit_(pv->desc, &nlocal, &nlocal, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);
    descinit_(pv->desc_wfc, &nlocal, &nband, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);
    descinit_(pv->desc_Eij, &nband, &nband, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);

    // Initialize data
    Stmp = new std::complex<double>[nlocal * nlocal];
    psi_k = new std::complex<double>[nlocal * nband];

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (i == j)
            {
                Stmp[i * nlocal + j] = std::complex<double>(1.0, 0.0);
            }
        }
    }
    Stmp[1] = 0.5;
    Stmp[4] = 0.5;

    psi_k[0] = 1.0;
    psi_k[1] = 1.0;
    psi_k[2] = 0.0;
    psi_k[3] = 0.0;
    psi_k[4] = 2.0;
    psi_k[5] = 1.0;
    psi_k[6] = 1.0;
    psi_k[7] = 0.0;
    psi_k[8] = 3.0;
    psi_k[9] = 0.0;
    psi_k[10] = 0.0;
    psi_k[11] = 1.0;

    // Call the function
    module_tddft::norm_psi(pv, nband, nlocal, Stmp, psi_k, print_matrix);

    // Check the results
    EXPECT_NEAR(psi_k[0].real(), 0.577350269189626, doublethreshold);
    EXPECT_NEAR(psi_k[0].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[1].real(), 0.577350269189626, doublethreshold);
    EXPECT_NEAR(psi_k[1].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[2].real(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[2].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[3].real(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[3].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[4].real(), 0.707106781186547, doublethreshold);
    EXPECT_NEAR(psi_k[4].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[5].real(), 0.353553390593274, doublethreshold);
    EXPECT_NEAR(psi_k[5].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[6].real(), 0.353553390593274, doublethreshold);
    EXPECT_NEAR(psi_k[6].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[7].real(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[7].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[8].real(), 0.948683298050514, doublethreshold);
    EXPECT_NEAR(psi_k[8].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[9].real(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[9].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[10].real(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[10].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(psi_k[11].real(), 0.316227766016838, doublethreshold);
    EXPECT_NEAR(psi_k[11].imag(), 0.0, doublethreshold);

    delete[] psi_k;
    delete[] Stmp;
}

#include <gtest/gtest.h>
#include <module_base/scalapack_connector.h>
#include <mpi.h>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/module_tddft/upsi.h"
#include "tddft_test.h"

#define doublethreshold 1e-8

/************************************************
 *  unit test of functions in upsi.h
 ***********************************************/

/**
 * - Tested Function
 *   - upsi
 *     - apply U_operator to the wave function of the previous step for new wave function.
 */

TEST(UpsiTest, testUpsi3)
{
    std::complex<double>* U_operator;
    std::complex<double>* psi_k_laststep;
    std::complex<double>* psi_k;
    int nband = 3;
    int nlocal = 4;
    bool print_matrix = false;
    Parallel_Orbitals* pv;
    pv = new Parallel_Orbitals();
    pv->nloc = nlocal * nband;
    pv->ncol = nlocal;
    pv->ncol_bands = nband;

    // Initialize input matrices
    int info;
    int mb = 1, nb = 1, lda = nband, ldc = nlocal;
    int irsrc = 0, icsrc = 0, lld = numroc_(&nlocal, &mb, &myprow, &irsrc, &nprow);
    descinit_(pv->desc, &nlocal, &nlocal, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);
    descinit_(pv->desc_wfc, &nlocal, &nband, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);

    // Initialize data
    U_operator = new std::complex<double>[nlocal * nlocal];
    psi_k_laststep = new std::complex<double>[nlocal * nband];
    psi_k = new std::complex<double>[nlocal * nband];

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (i == j)
            {
                U_operator[i * nlocal + j] = std::complex<double>(1.0, 0.0);
            }
            else
            {
                U_operator[i * nlocal + j] = std::complex<double>(0.0, 0.0);
            }
        }
    }

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nband; ++j)
        {
            psi_k_laststep[i * nband + j] = std::complex<double>(i * nband + j, 0.0);
            psi_k[i * nband + j] = std::complex<double>(0.0, 0.0);
        }
    }

    // Call the function
    module_tddft::upsi(pv, nband, nlocal, U_operator, psi_k_laststep, psi_k, print_matrix);

    // Check the results
    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nband; ++j)
        {
            EXPECT_NEAR(psi_k[i * nband + j].real(), i * nband + j, doublethreshold);
            EXPECT_NEAR(psi_k[i].imag(), 0.0, doublethreshold);
        }
    }
    delete[] U_operator;
    delete[] psi_k;
    delete[] psi_k_laststep;
}

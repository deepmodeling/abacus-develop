#include <gtest/gtest.h>
#include <module_base/scalapack_connector.h>
#include <mpi.h>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/module_tddft/propagator.h"
#include "module_io/input.h"
#include "tddft_test.h"

/************************************************
 *  unit test of functions in propagator.h
 ***********************************************/

/**
 * - Tested Function
 *   - Propagator::compute_propagator_etrs
 *     - compute propagator of method ETRS.
 */

#define doublethreshold 1e-8

TEST(PropagatorTest, testPropagatorETRS)
{
    std::complex<double>* U_operator;
    std::complex<double>* Stmp;
    std::complex<double>* Htmp;
    std::complex<double>* Hlaststep;
    int nlocal = 4;
    bool print_matrix = false;
    Parallel_Orbitals* pv;
    pv = new Parallel_Orbitals();
    pv->nloc = nlocal * nlocal;
    pv->ncol = nlocal;
    pv->nrow = nlocal;
    pv->dim0 = 1;
    pv->dim1 = 1;
    pv->nb = 1;

    int dim[2];
    int period[2] = {1, 1};
    int reorder = 0;
    dim[0] = nprow;
    dim[1] = npcol;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &pv->comm_2D);

    INPUT.mdp.md_dt = 4;

    // Initialize input matrices
    int info;
    int mb = 1, nb = 1, lda = nlocal, ldc = nlocal;
    int irsrc = 0, icsrc = 0, lld = numroc_(&nlocal, &mb, &myprow, &irsrc, &nprow);
    descinit_(pv->desc, &nlocal, &nlocal, &mb, &nb, &irsrc, &icsrc, &ictxt, &lld, &info);

    // Initialize data
    U_operator = new std::complex<double>[nlocal * nlocal];
    Stmp = new std::complex<double>[nlocal * nlocal];
    Htmp = new std::complex<double>[nlocal * nlocal];
    Hlaststep = new std::complex<double>[nlocal * nlocal];

    for (int i = 0; i < nlocal * nlocal; ++i)
    {
        U_operator[i] = std::complex<double>(0.0, 0.0);
    }
    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (i == j)
            {
                Htmp[i * nlocal + j] = std::complex<double>(1.0, 0.0);
                Hlaststep[i * nlocal + j] = std::complex<double>(1.0, 0.0);
                Stmp[i * nlocal + j] = std::complex<double>(1.0, 0.0);
            }
        }
    }
    Stmp[1] = 0.5;
    Stmp[4] = 0.5;
    Hlaststep[0] = 1.1;

    // Call the function
    int propagator = 2;
    module_tddft::Propagator prop(propagator, pv);
    prop.compute_propagator(nlocal, Stmp, Htmp, Hlaststep, U_operator, print_matrix);

    // Check the results
    EXPECT_NEAR(U_operator[0].real(), -0.0105865569272976, doublethreshold);
    EXPECT_NEAR(U_operator[0].imag(), -0.228336412132297, doublethreshold);
    EXPECT_NEAR(U_operator[1].real(), 0.195527023319616, doublethreshold);
    EXPECT_NEAR(U_operator[1].imag(), -0.728701844231062, doublethreshold);
    EXPECT_NEAR(U_operator[2].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[2].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[3].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[3].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[4].real(), 0.247662246608749, doublethreshold);
    EXPECT_NEAR(U_operator[4].imag(), -0.694263679317177, doublethreshold);
    EXPECT_NEAR(U_operator[5].real(), -0.0219180003048316, doublethreshold);
    EXPECT_NEAR(U_operator[5].imag(), -0.300090839811004, doublethreshold);
    EXPECT_NEAR(U_operator[6].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[6].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[7].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[7].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[8].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[8].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[9].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[9].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[10].real(), -0.401041666666667, doublethreshold);
    EXPECT_NEAR(U_operator[10].imag(), -0.902777777777778, doublethreshold);
    EXPECT_NEAR(U_operator[11].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[11].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[12].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[12].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[13].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[13].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[14].real(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[14].imag(), 0.0, doublethreshold);
    EXPECT_NEAR(U_operator[15].real(), -0.401041666666667, doublethreshold);
    EXPECT_NEAR(U_operator[15].imag(), -0.902777777777778, doublethreshold);

    delete[] U_operator;
    delete[] Htmp;
    delete[] Stmp;
    delete[] Hlaststep;
}
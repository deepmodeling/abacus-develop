#include "module_hamilt_lcao/module_tddft/bandenergy.h"

#include <gtest/gtest.h>
#include <module_base/scalapack_connector.h>
#include <mpi.h>

#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "tddft_test.h"

/************************************************
 *  unit test of functions in bandenergy.h
 ***********************************************/

/**
 * - Tested Function
 *   - compute_ekb
 *     - compute band energy ekb <psi_i|H|psi_i>.
 */

#define doublethreshold 1e-8
double module_tddft::Evolve_elec::td_print_eij = -1;

Parallel_Orbitals::Parallel_Orbitals()
{
}
Parallel_Orbitals::~Parallel_Orbitals()
{
}

TEST(BandEnergyTest, testBandEnergy)
{
    std::complex<double>* psi_k;
    std::complex<double>* Htmp;
    double* ekb;
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
    Htmp = new std::complex<double>[nlocal * nlocal];
    psi_k = new std::complex<double>[nlocal * nband];
    ekb = new double[nband];

    for (int i = 0; i < nlocal; ++i)
    {
        for (int j = 0; j < nlocal; ++j)
        {
            if (i == j)
            {
                Htmp[i * nlocal + j] = std::complex<double>(1.0, 0.0);
            }
        }
    }
    Htmp[1] = 0.5;
    Htmp[4] = 0.5;

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
    module_tddft::compute_ekb(pv, nband, nlocal, Htmp, psi_k, ekb);

    // Check the results
    EXPECT_NEAR(ekb[0], 3.0, doublethreshold);
    EXPECT_NEAR(ekb[1], 8.0, doublethreshold);
    EXPECT_NEAR(ekb[2], 10.0, doublethreshold);

    delete[] psi_k;
    delete[] Htmp;
    delete[] ekb;
}

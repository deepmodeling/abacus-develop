#include "source_lcao/module_rt/band_energy.h"

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_io/module_parameter/parameter.h"
#include "tddft_test.h"

#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <mpi.h>
#include <source_base/module_external/scalapack_connector.h>

class TestParameters
{
  public:
    static void set_td_print_eij(const double threshold)
    {
        PARAM.input.td_print_eij = threshold;
    }
};

/************************************************
 *  unit test of functions in band_energy.h
 ***********************************************/

/**
 * - Tested Function
 *   - compute_ekb
 *     - compute band energy ekb <psi_i|H|psi_i>.
 */

#define doublethreshold 1e-8

TEST(BandEnergyTest, testBandEnergy)
{
    std::complex<double>* psi_k;
    std::complex<double>* Htmp;
    double* ekb;
    int nband = 3;
    int nlocal = 4;
    Parallel_Orbitals* pv;
    pv = new Parallel_Orbitals();
    pv->nloc = nlocal * nlocal;
    pv->nloc_wfc = nlocal * nband;
    pv->ncol = nlocal;
    pv->nrow = nlocal;
    pv->ncol_bands = nband;
    pv->nrow_bands = nband;
    pv->dim0 = 1;
    pv->dim1 = 1;
    pv->nb = 1;
    pv->blacs_ctxt = 0;
    pv->set_coord(0, 0);

    int dim[2];
    dim[0] = nprow;
    dim[1] = npcol;

    // Initialize input matrices
    int info;
    int mb = 1, nb = 1, lda = nband, ldc = nlocal;
    int irsrc = 0, icsrc = 0, lld = numroc_(&nlocal, &mb, &myprow, &irsrc, &nprow), lld1 = numroc_(&nband, &mb, &myprow, &irsrc, &nprow);
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

    // Call the function with local ofstream
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string output_path = "band_energy_test_output_" + std::to_string(rank) + ".txt";
    std::ofstream ofs(output_path.c_str());
    TestParameters::set_td_print_eij(0.0);
    module_rt::compute_ekb(pv, nband, nlocal, Htmp, psi_k, ekb, ofs);
    ofs.close();
    TestParameters::set_td_print_eij(-1.0);

    // Check the results
    EXPECT_NEAR(ekb[0], 3.0, doublethreshold);
    EXPECT_NEAR(ekb[1], 8.0, doublethreshold);
    EXPECT_NEAR(ekb[2], 10.0, doublethreshold);

    std::ifstream output(output_path.c_str());
    const std::string output_text((std::istreambuf_iterator<char>(output)), std::istreambuf_iterator<char>());
    EXPECT_NE(output_text.find("i =  1"), std::string::npos);
    EXPECT_NE(output_text.find("j =  1"), std::string::npos);
    EXPECT_EQ(output_text.find("i =  0"), std::string::npos);
    EXPECT_EQ(output_text.find("j =  0"), std::string::npos);
    output.close();
    std::remove(output_path.c_str());

    delete[] psi_k;
    delete[] Htmp;
    delete[] ekb;
}

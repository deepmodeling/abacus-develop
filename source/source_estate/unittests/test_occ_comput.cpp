#include "gtest/gtest.h"
#include "source_estate/occ_comput.h"

#include <complex>
#include <vector>

/***************************************************************
 * unit test of elecstate::occ_from_proj
 *
 * The function accumulates, for one k-point, the per-projector 2x2
 * occupation blocks
 *   rho^{ss'}_{iprj} = sum_i w_i conj(proj^s_{i,iprj}) proj^{s'}_{i,iprj}
 * stored as occ_block[iprj*4 + {0,1,2,3}] = {up-up, up-dn, dn-up, dn-dn}.
 ****************************************************************/

class OccComputTest : public ::testing::Test
{
  protected:
    // two atoms with nh = {2, 3} -> tot_nproj = 5
    const int nat = 2;
    const int nh[2] = {2, 3};
    const int nkb = 5;
    const int nbands = 2;

    // proj layout: (nbands*npol) x nkb, spinor components offset by nkb.
    // Values are simple distinct numbers so wrong indexing is caught.
    std::vector<std::complex<double>> make_proj(const int npol)
    {
        std::vector<std::complex<double>> proj(nbands * npol * nkb);
        for (size_t i = 0; i < proj.size(); i++)
        {
            proj[i] = std::complex<double>(0.1 * (i + 1), 0.01 * (i + 1));
        }
        return proj;
    }

    // reference implementation copied from the original inline loops in
    // OnsiteProjector::cal_occupations (onsite_proj_overlap.cpp), used as
    // an oracle for the extracted free function.
    void reference(const std::complex<double>* proj,
                   const double* wg_ik,
                   const int npol,
                   const int nspin,
                   const int isk,
                   std::vector<std::complex<double>>& occs)
    {
        for (int ib = 0; ib < nbands; ib++)
        {
            const double weight = wg_ik[ib];
            int begin_iprj = 0;
            for (int iat = 0; iat < nat; iat++)
            {
                const int nprj = nh[iat];
                for (int iprj = 0; iprj < nprj; iprj++)
                {
                    const int occ_index = (begin_iprj + iprj) * 4;
                    if (npol == 1)
                    {
                        const int index = ib * nkb + begin_iprj + iprj;
                        const double occ = weight * (std::conj(proj[index]) * proj[index]).real();
                        if (nspin == 2 && isk == 1)
                        {
                            occs[occ_index + 3] += occ;
                        }
                        else if (nspin == 1)
                        {
                            occs[occ_index] += 0.5 * occ;
                            occs[occ_index + 3] += 0.5 * occ;
                        }
                        else
                        {
                            occs[occ_index] += occ;
                        }
                    }
                    else
                    {
                        const int index = ib * 2 * nkb + begin_iprj + iprj;
                        occs[occ_index] += weight * std::conj(proj[index]) * proj[index];
                        occs[occ_index + 1] += weight * std::conj(proj[index]) * proj[index + nkb];
                        occs[occ_index + 2] += weight * std::conj(proj[index + nkb]) * proj[index];
                        occs[occ_index + 3] += weight * std::conj(proj[index + nkb]) * proj[index + nkb];
                    }
                }
                begin_iprj += nprj;
            }
        }
    }

    void run_and_compare(const int npol, const int nspin, const int isk)
    {
        const double wg[2] = {1.5, 0.5};
        std::vector<std::complex<double>> proj = make_proj(npol);
        std::vector<std::complex<double>> got(nkb * 4, std::complex<double>(0.0, 0.0));
        std::vector<std::complex<double>> want(nkb * 4, std::complex<double>(0.0, 0.0));

        elecstate::occ_from_proj(proj.data(), wg, nbands, npol, nkb,
                                 nspin, isk, nh, nat, got.data());
        reference(proj.data(), wg, npol, nspin, isk, want);

        for (int i = 0; i < nkb * 4; i++)
        {
            EXPECT_NEAR(got[i].real(), want[i].real(), 1e-12) << "element " << i;
            EXPECT_NEAR(got[i].imag(), want[i].imag(), 1e-12) << "element " << i;
        }
    }
};

TEST_F(OccComputTest, Nspin1SplitsEvenly)
{
    run_and_compare(1, 1, 0);

    // additionally verify the even split and zero magnetization explicitly
    const double wg[2] = {1.0, 1.0};
    std::vector<std::complex<double>> proj = make_proj(1);
    std::vector<std::complex<double>> occs(nkb * 4, std::complex<double>(0.0, 0.0));
    elecstate::occ_from_proj(proj.data(), wg, nbands, 1, nkb, 1, 0, nh, nat, occs.data());
    for (int iprj = 0; iprj < nkb; iprj++)
    {
        EXPECT_NEAR(occs[iprj * 4].real(), occs[iprj * 4 + 3].real(), 1e-12);
        EXPECT_NEAR(occs[iprj * 4 + 1].real(), 0.0, 1e-12);
        EXPECT_NEAR(occs[iprj * 4 + 2].real(), 0.0, 1e-12);
    }
}

TEST_F(OccComputTest, Nspin2SpinUpGoesToBlock0)
{
    run_and_compare(1, 2, 0);

    const double wg[2] = {1.0, 1.0};
    std::vector<std::complex<double>> proj = make_proj(1);
    std::vector<std::complex<double>> occs(nkb * 4, std::complex<double>(0.0, 0.0));
    elecstate::occ_from_proj(proj.data(), wg, nbands, 1, nkb, 2, 0, nh, nat, occs.data());
    for (int iprj = 0; iprj < nkb; iprj++)
    {
        EXPECT_GT(occs[iprj * 4].real(), 0.0);
        EXPECT_NEAR(occs[iprj * 4 + 3].real(), 0.0, 1e-12);
    }
}

TEST_F(OccComputTest, Nspin2SpinDownGoesToBlock3)
{
    run_and_compare(1, 2, 1);

    const double wg[2] = {1.0, 1.0};
    std::vector<std::complex<double>> proj = make_proj(1);
    std::vector<std::complex<double>> occs(nkb * 4, std::complex<double>(0.0, 0.0));
    elecstate::occ_from_proj(proj.data(), wg, nbands, 1, nkb, 2, 1, nh, nat, occs.data());
    for (int iprj = 0; iprj < nkb; iprj++)
    {
        EXPECT_NEAR(occs[iprj * 4].real(), 0.0, 1e-12);
        EXPECT_GT(occs[iprj * 4 + 3].real(), 0.0);
    }
}

TEST_F(OccComputTest, Nspin4FillsAllBlocks)
{
    run_and_compare(2, 4, 0);

    // off-diagonal blocks must be nonzero for generic spinor coefficients
    const double wg[2] = {1.0, 1.0};
    std::vector<std::complex<double>> proj = make_proj(2);
    std::vector<std::complex<double>> occs(nkb * 4, std::complex<double>(0.0, 0.0));
    elecstate::occ_from_proj(proj.data(), wg, nbands, 2, nkb, 4, 0, nh, nat, occs.data());
    for (int iprj = 0; iprj < nkb; iprj++)
    {
        EXPECT_GT(occs[iprj * 4].real(), 0.0);
        EXPECT_GT(occs[iprj * 4 + 3].real(), 0.0);
        // hermiticity: rho^{up,dn} = conj(rho^{dn,up})
        EXPECT_NEAR(occs[iprj * 4 + 1].real(), occs[iprj * 4 + 2].real(), 1e-12);
        EXPECT_NEAR(occs[iprj * 4 + 1].imag(), -occs[iprj * 4 + 2].imag(), 1e-12);
    }
}

TEST_F(OccComputTest, AccumulatesAcrossCalls)
{
    // two calls with weight 1 must equal one call with weight 2 per band
    const double wg1[2] = {1.0, 1.0};
    const double wg2[2] = {2.0, 2.0};
    std::vector<std::complex<double>> proj = make_proj(1);
    std::vector<std::complex<double>> twice(nkb * 4, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> once(nkb * 4, std::complex<double>(0.0, 0.0));

    elecstate::occ_from_proj(proj.data(), wg1, nbands, 1, nkb, 2, 0, nh, nat, twice.data());
    elecstate::occ_from_proj(proj.data(), wg1, nbands, 1, nkb, 2, 0, nh, nat, twice.data());
    elecstate::occ_from_proj(proj.data(), wg2, nbands, 1, nkb, 2, 0, nh, nat, once.data());

    for (int i = 0; i < nkb * 4; i++)
    {
        EXPECT_NEAR(twice[i].real(), once[i].real(), 1e-12);
    }
}

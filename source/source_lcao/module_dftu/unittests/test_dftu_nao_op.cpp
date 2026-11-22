#include "gtest/gtest.h"
#include <algorithm>
#include <chrono>
#include <complex>
#include <string>
#include <vector>

// mock of DFTU
#include "../dftu_full_dc.h"
#include "../dftu_nao_op.h"
#include "source_base/constants.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_estate/module_dm/density_matrix.h"

Plus_U_Base dftu;

//---------------------------------------
// Unit test of DFTU_onsite operator (dftu_nao_op.cpp).
// Tests constructHR for d2d and d2cd variants.
//---------------------------------------

int test_size = 5;
int test_nw = 5;

class DFTUTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
#ifdef __MPI
        MPI_Comm_size(MPI_COMM_WORLD, &dsize);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

        ucell.ntype = 1;
        ucell.nat = test_size;
        ucell.atoms = new Atom[ucell.ntype];
        ucell.iat2it.resize(ucell.nat);
        ucell.iat2ia.resize(ucell.nat);
        ucell.atoms[0].tau.resize(ucell.nat);
        ucell.lat0 = 1.0;
        ucell.itia2iat.create(ucell.ntype, ucell.nat);
        for (int iat = 0; iat < ucell.nat; iat++)
        {
            ucell.iat2it[iat] = 0;
            ucell.iat2ia[iat] = iat;
            ucell.atoms[0].tau[iat] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
            ucell.itia2iat(0, iat) = iat;
        }
        ucell.atoms[0].na = test_size;
        ucell.atoms[0].nw = test_nw;
        ucell.atoms[0].iw2l.resize(test_nw);
        ucell.atoms[0].iw2m.resize(test_nw);
        ucell.atoms[0].iw2n.resize(test_nw);
        for (int iw = 0; iw < test_nw; ++iw)
        {
            ucell.atoms[0].iw2l[iw] = 2;
            ucell.atoms[0].iw2m[iw] = 0;
            ucell.atoms[0].iw2n[iw] = 0;
        }
        ucell.set_iat2iwt(1);
        init_parav();
        HR = new hamilt::HContainer<double>(ucell, paraV);

        dftu.occmat().data().resize(test_size);
        for (int iat = 0; iat < test_size; iat++)
        {
            dftu.occmat().data()[iat].resize(3);
            for (int l = 0; l < 3; l++)
            {
                dftu.occmat().data()[iat][l].resize(2);
                dftu.occmat().data()[iat][l][0].create(2 * l + 1, 2 * l + 1);
                dftu.occmat().data()[iat][l][1].create(2 * l + 1, 2 * l + 1);
            }
        }
        dftu.u_current = {U_test};
        dftu.l_channel = {orbital_c_test};
    }

    void TearDown() override
    {
        delete HR;
        delete paraV;
        delete[] ucell.atoms;
    }

    double occ_mat_c(int iat, int spin, int icc) const
    {
        return dftu.occmat().data()[iat][2][spin].c[icc];
    }

    template <typename TK>
    void check_density_built_spin_factors()
    {
        const std::vector<std::string> modes = {"fll", "amf", "mixed", "mixed", "mixed"};
        const std::vector<double> weights = {1.0, 0.0, 0.0, 0.35, 1.0};
        const std::vector<double> cutoffs(1, 1.0);
        const int shell_size = 5;
        const double density_terms = test_nw * test_nw * test_size * test_size;
        Grid_Driver grid(0, 0);
        for (int nspin : {1, 2})
        {
            SCOPED_TRACE(nspin);
            std::vector<ModuleBase::Vector3<double>> kvec(nspin, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
            hamilt::HS_Matrix_K<TK> hsk(paraV, true);
            for (std::size_t mode = 0; mode < modes.size(); ++mode)
            {
                SCOPED_TRACE(modes[mode]);
                SCOPED_TRACE(weights[mode]);
                module_dm::DensityMatrix<TK, double> dm(paraV, nspin, nspin);
                dm.init_dmr(*HR);
                hamilt::DFTU_onsite<hamilt::OperatorLCAO<TK, double>> op(
                    &hsk, kvec, HR, ucell, &grid, &intor_, cutoffs, &dftu, nspin,
                    onsite_radius_test, &dm, modes[mode], weights[mode]);
                auto& matrices = dm.get_dmk_vec();
                ASSERT_EQ(matrices.size(), static_cast<std::size_t>(nspin));
                // Existing fixture overlaps are all one. DMK -> DMR must yield
                // a shell with every entry p, not a preloaded occupation.
                for (int round = 0; round < 2; ++round)
                {
                    SCOPED_TRACE(round);
                    std::vector<double> occupations(nspin, 0.08);
                    for (int spin = 0; spin < nspin; ++spin)
                    {
                        if (round == 1)
                        {
                            occupations[spin] = 0.06 - 0.02 * spin;
                        }
                        const double degeneracy = (nspin == 1) ? 2.0 : 1.0;
                        const TK density = degeneracy * occupations[spin] / density_terms;
                        std::fill(matrices[spin].begin(), matrices[spin].end(), density);
                    }
                    dm.cal_dmr(-1);
                    ASSERT_TRUE(dm.is_dmr_ready());
                    dftu.set_occmat_stale();
                    double expected_energy = 0.0;
                    for (int spin = 0; spin < nspin; ++spin)
                    {
                        SCOPED_TRACE(spin);
                        HR->set_zero();
                        op.contributeHR();
                        const double p = occupations[spin];
                        const double alpha = weights[mode];
                        const double q = shell_size * shell_size - (1.0 - alpha) * shell_size;
                        const double degeneracy = (nspin == 1) ? 2.0 : 1.0;
                        expected_energy += degeneracy * test_size * 0.5 * U_test * q * p * p;
                        const double expected_hr = test_size * U_test * (0.5 * shell_size * alpha - q * p);
                        for (int atom = 0; atom < test_size; ++atom)
                        {
                            for (int entry = 0; entry < shell_size * shell_size; ++entry)
                            {
                                EXPECT_NEAR(occ_mat_c(atom, spin, entry), p, 1e-11);
                            }
                        }
                        ASSERT_GT(HR->get_nnr(), 0);
                        for (int entry = 0; entry < HR->get_nnr(); ++entry)
                        {
                            EXPECT_NEAR(HR->get_wrapper()[entry], expected_hr, 1e-10);
                        }
                        EXPECT_NEAR(dftu.get_energy(), expected_energy, 1e-11);
                    }
                }
            }
        }
        dftu.set_occmat_stale();
    }

#ifdef __MPI
    void init_parav()
    {
        int nb = 10;
        int global_row = test_size * test_nw;
        int global_col = test_size * test_nw;
        std::ofstream ofs_running;
        paraV = new Parallel_Orbitals();
        paraV->para_2d().init(global_row, global_col, nb, MPI_COMM_WORLD, false);
        paraV->set_atomic_trace(ucell.get_iat2iwt(), test_size, global_row);
    }
#else
    void init_parav()
    {
    }
#endif

    UnitCell ucell;
    hamilt::HContainer<double>* HR;
    Parallel_Orbitals* paraV;
    TwoCenterIntegrator intor_;

    int dsize;
    int my_rank = 0;
    double U_test = 1.0;
    int orbital_c_test = 2;
    double onsite_radius_test = 1.0;
};

TEST_F(DFTUTest, constructHRd2d)
{
    const int nspin = 1;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    const double factor = 1.0 / test_nw / test_nw / test_size / test_size;
    module_dm::DensityMatrix<double, double> dm(paraV, 1);
    dm.init_dmr(*HR);
    for (int i = 0; i < paraV->para_2d().nrloc(); i++)
    {
        for (int j = 0; j < paraV->para_2d().ncloc(); j++)
        {
            dm.set_dmk(1, 0, i, j, factor);
        }
    }
    dm.cal_dmr(-1);
    // reset HR
    for (int i = 0; i < HR->get_nnr(); i++)
    {
        HR->get_wrapper()[i] = 0.0;
    }
    hamilt::DFTU_onsite<hamilt::OperatorLCAO<double, double>>
        op(&hsk, kvec_d_in, HR, ucell, &gd, &intor_, {1.0}, &dftu, nspin, onsite_radius_test, &dm);
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 0, icc), 0.5, 1e-10);
        }
    }
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0 * test_size, 1e-10);
        }
    }
    op.contributeHk(0);
    double* hk = hsk.get_hk();
    for (int i = 0; i < paraV->para_2d().nrloc() * paraV->para_2d().ncloc(); ++i)
    {
        EXPECT_NEAR(hk[i], -10.0 * test_size, 1e-10);
    }
}

TEST_F(DFTUTest, constructHRd2cd)
{
    const int nspin = 2;
    std::vector<ModuleBase::Vector3<double>> kvec_d_in(2, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<std::complex<double>> hsk(paraV, true);
    hsk.set_zero_hk();
    Grid_Driver gd(0, 0);
    const double factor = 0.5 / test_nw / test_nw / test_size / test_size;
    std::vector<ModuleBase::Vector3<double>> kvec_d_dm(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    module_dm::DensityMatrix<std::complex<double>, double> dm(paraV, 2, kvec_d_dm, 1);
    dm.init_dmr(*HR);
    for (int is = 1; is <= 2; ++is)
    {
        for (int i = 0; i < paraV->para_2d().nrloc(); i++)
        {
            for (int j = 0; j < paraV->para_2d().ncloc(); j++)
            {
                dm.set_dmk(is, 0, i, j, std::complex<double>(factor, 0.0));
            }
        }
    }
    dm.cal_dmr(-1);
    // reset HR
    for (int i = 0; i < HR->get_nnr(); i++)
    {
        HR->get_wrapper()[i] = 0.0;
    }
    hamilt::DFTU_onsite<hamilt::OperatorLCAO<std::complex<double>, double>>
        op(&hsk, kvec_d_in, HR, ucell, &gd, &intor_, {1.0}, &dftu, nspin, onsite_radius_test, &dm);
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 0, icc), 0.5, 1e-10);
        }
    }
    for (int iap = 0; iap < HR->size_atom_pairs(); ++iap)
    {
        hamilt::AtomPair<double>& tmp = HR->get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        std::vector<int> indexes1 = paraV->get_indexes_row(iat1);
        std::vector<int> indexes2 = paraV->get_indexes_col(iat2);
        int nwt = indexes1.size() * indexes2.size();
        for (int i = 0; i < nwt; ++i)
        {
            EXPECT_NEAR(tmp.get_pointer(0)[i], -10.0 * test_size, 1e-10);
        }
    }
    op.contributeHk(0);
    std::complex<double>* hk = hsk.get_hk();
    for (int i = 0; i < paraV->para_2d().nrloc() * paraV->para_2d().ncloc(); ++i)
    {
        EXPECT_NEAR(hk[i].real(), -10.0 * test_size, 1e-10);
        EXPECT_NEAR(hk[i].imag(), 0.0, 1e-10);
    }
    op.contributeHR();
    for (int iat = 0; iat < test_size; iat++)
    {
        for (int icc = 0; icc < 25; icc++)
        {
            EXPECT_NEAR(occ_mat_c(iat, 1, icc), 0.5, 1e-10);
        }
    }
}

TEST_F(DFTUTest, FullDcDispatchEndpointsAndLegacyReset)
{
    std::vector<ModuleBase::Vector3<double>> kvec(1, ModuleBase::Vector3<double>(0.0, 0.0, 0.0));
    hamilt::HS_Matrix_K<double> hsk(paraV, true);
    Grid_Driver grid(0, 0);
    module_dm::DensityMatrix<double, double> dm(paraV, 1);
    dm.init_dmr(*HR);
    // Exercise the real production operator using a nonuniform symmetric shell.
    // Preloaded in memory only; no restart file is used by this unit fixture.
    for (int atom = 0; atom < test_size; ++atom)
    {
        auto& n = dftu.occmat().data()[atom][2][0];
        for (int i = 0; i < 5; ++i)
        {
            for (int j = 0; j < 5; ++j)
            {
                n(i, j) = (i == j) ? 0.15 + 0.1 * i : 0.01 * (i + j + 1);
            }
        }
    }
    std::vector<std::vector<double>> matrices;
    std::vector<double> energies;
    const std::vector<std::string> modes = {"fll", "mixed", "amf", "mixed", "mixed", "fll"};
    const std::vector<double> weights = {1.0, 1.0, 0.0, 0.0, 0.5, 1.0};
    for (std::size_t pass = 0; pass < modes.size(); ++pass)
    {
        std::fill(HR->get_wrapper(), HR->get_wrapper() + HR->get_nnr(), 0.0);
        dftu.set_occmat_ready();
        if (pass == 0)
        {
            hamilt::DFTU_onsite<hamilt::OperatorLCAO<double, double>> op(
                &hsk, kvec, HR, ucell, &grid, &intor_, {1.0}, &dftu, 1, onsite_radius_test, &dm);
            op.contributeHR();
        }
        else
        {
            hamilt::DFTU_onsite<hamilt::OperatorLCAO<double, double>> op(
                &hsk, kvec, HR, ucell, &grid, &intor_, {1.0}, &dftu, 1, onsite_radius_test, &dm,
                modes[pass], weights[pass]);
            op.contributeHR();
        }
        matrices.emplace_back(HR->get_wrapper(), HR->get_wrapper() + HR->get_nnr());
        energies.push_back(dftu.get_energy());
    }
    ASSERT_FALSE(matrices[0].empty());
    double change = 0.0;
    for (std::size_t i = 0; i < matrices[0].size(); ++i)
    {
        EXPECT_NEAR(matrices[0][i], matrices[1][i], 1e-12);
        EXPECT_NEAR(matrices[2][i], matrices[3][i], 1e-12);
        EXPECT_NEAR(matrices[4][i], 0.5 * (matrices[0][i] + matrices[2][i]), 1e-12);
        EXPECT_DOUBLE_EQ(matrices[0][i], matrices[5][i]);
        change += std::abs(matrices[0][i] - matrices[2][i]);
    }
    EXPECT_GT(change, 1e-8);
    EXPECT_NEAR(energies[0], energies[1], 1e-12);
    EXPECT_NEAR(energies[2], energies[3], 1e-12);
    EXPECT_NEAR(energies[4], 0.5 * (energies[0] + energies[2]), 1e-12);
    EXPECT_DOUBLE_EQ(energies[0], energies[5]);
    dftu.set_occmat_stale();
}

TEST_F(DFTUTest, FullDcSpinFactorsFromRealDensity)
{
    check_density_built_spin_factors<double>();
}

TEST_F(DFTUTest, FullDcSpinFactorsFromComplexDensity)
{
    check_density_built_spin_factors<std::complex<double>>();
}

TEST_F(DFTUTest, FullDcConjugateKDensityAndHamiltonian)
{
    using Complex = std::complex<double>;
    const int shell = 5;
    const int shell_entries = shell * shell;
    const int basis = test_size * test_nw;
    const auto& layout = paraV->para_2d();
    const int rows = layout.nrloc();
    const int cols = layout.ncloc();
    const int local_entries = rows * cols;
    const double angle = ModuleBase::TWO_PI * 0.23;
    const Complex phase = std::polar(1.0, angle);
    const std::vector<ModuleBase::Vector3<double>> kvec = {{0.23, 0.0, 0.0}, {-0.23, 0.0, 0.0}};

    // Supply overlaps and neighbors explicitly: the shared grid/integrator mock
    // only yields R=0. The contractions, DMK->DMR and H(R)->H(k) below are real
    // production functions, with atoms on different block-cyclic MPI owners.
    AdjacentAtomInfo adjs;
    adjs.adj_num = 1;
    adjs.ntype = {0, 0};
    adjs.natom = {0, 3};
    adjs.box = {{0, 0, 0}, {1, 0, 0}};
    DFTU_LCAO::NlmTot overlaps(test_size);
    overlaps[0].resize(2);
    hamilt::HContainer<double> hr(paraV);
    std::vector<Complex> coefficients(basis, 0.0);
    const std::vector<Complex> empty_projector(shell, 0.0);
    std::vector<std::vector<Complex>> projectors(basis, empty_projector);
    std::vector<Complex> projected(shell, 0.0);
    std::vector<Complex> gamma_projected(shell, 0.0);
    for (int ad = 0; ad < 2; ++ad)
    {
        const int atom = adjs.natom[ad];
        const Complex factor = (ad == 0) ? Complex(1.0, 0.0) : phase;
        for (int iw = 0; iw < test_nw; ++iw)
        {
            const int g = atom * test_nw + iw;
            const double real = 0.15 + 0.02 * g;
            const double imag = 0.08 * (iw + 1) * (ad == 0 ? 1.0 : -1.0);
            coefficients[g] = Complex(real, imag);
            auto& overlap = overlaps[0][ad][iw];
            overlap.resize(shell);
            for (int m = 0; m < shell; ++m)
            {
                overlap[m] = 0.06 * (1 + iw + m) + 0.02 * ad * (m + 1) + (iw == m ? 0.1 : 0.0);
                projectors[g][m] = overlap[m] * factor;
                projected[m] += coefficients[g] * projectors[g][m];
                gamma_projected[m] += coefficients[g] * overlap[m];
            }
        }
        for (int other = 0; other < 2; ++other)
        {
            const int atom_j = adjs.natom[other];
            if (paraV->get_nrow_atom(atom) == 0 || paraV->get_ncol_atom(atom_j) == 0)
            {
                continue;
            }
            const int rx = other - ad;
            hamilt::AtomPair<double> pair(atom, atom_j, rx, 0, 0, paraV);
            hr.insert_pair(pair);
        }
    }
    hr.allocate(nullptr, true);
    module_dm::DensityMatrix<Complex, double> dm(paraV, 1, kvec, 2, 1);
    dm.init_dmr(hr);
    auto& dmk = dm.get_dmk_vec();
    for (int j = 0; j < cols; ++j)
    {
        const int gj = layout.icglb(j);
        for (int i = 0; i < rows; ++i)
        {
            const int gi = layout.irglb(i);
            const int index = j * rows + i;
            // Stored DMK is conjugated physical D and already includes the
            // k weight (1/2) and nonmagnetic spin degeneracy (2).
            dmk[0][index] = std::conj(coefficients[gi]) * coefficients[gj];
            dmk[1][index] = std::conj(dmk[0][index]);
        }
    }
    dm.cal_dmr(-1);
    const auto& dmr = *dm.get_dmr_ptr(1);
    std::vector<double> occupation(shell_entries, 0.0);
    DFTU_LCAO::compute_occ_from_dmr(ucell, dftu, 0, 2, 0, 1, adjs, *paraV, overlaps, dmr, occupation);

    std::vector<double> expected_occ(shell_entries, 0.0);
    double phase_sensitivity = 0.0;
    double trace = 0.0;
    double square_trace = 0.0;
    for (int m = 0; m < shell; ++m)
    {
        for (int n = 0; n < shell; ++n)
        {
            const int mn = m * shell + n;
            const Complex product = std::conj(projected[m]) * projected[n];
            expected_occ[mn] = product.real();
            const Complex gamma_product = std::conj(gamma_projected[m]) * gamma_projected[n];
            const double difference = expected_occ[mn] - gamma_product.real();
            phase_sensitivity += std::abs(difference);
            EXPECT_NEAR(occupation[mn], expected_occ[mn], 1e-11);
            square_trace += expected_occ[mn] * expected_occ[mn];
        }
        trace += expected_occ[m * shell + m];
    }
    EXPECT_GT(phase_sensitivity, 1e-3);

    for (double alpha : {0.0, 0.35, 1.0})
    {
        SCOPED_TRACE(alpha);
        std::vector<double> potential;
        double energy_correction = 0.0;
        const double energy = DFTU_LCAO::evaluate_full_dc(occupation, shell, U_test, alpha, potential, energy_correction);
        const double shift = alpha * 0.5 + (1.0 - alpha) * trace / shell;
        const double interaction = 0.5 * U_test
            * (alpha * trace + (1.0 - alpha) * trace * trace / shell - square_trace);
        std::vector<double> expected_v(shell_entries, 0.0);
        double shell_trace = 0.0;
        for (int m = 0; m < shell; ++m)
        {
            for (int n = 0; n < shell; ++n)
            {
                const int mn = m * shell + n;
                expected_v[mn] = U_test * ((m == n ? shift : 0.0) - expected_occ[mn]);
                EXPECT_NEAR(potential[mn], expected_v[mn], 1e-11);
                shell_trace += expected_v[mn] * expected_occ[n * shell + m];
            }
        }
        const double correction = interaction - shell_trace;
        EXPECT_NEAR(energy, interaction, 1e-11);
        EXPECT_NEAR(energy_correction, correction, 1e-11);
        // Reset/rebuild, rather than adding twice to an already built H(R).
        for (int rebuild = 0; rebuild < 2; ++rebuild)
        {
            SCOPED_TRACE(rebuild);
            hr.set_zero();
            DFTU_LCAO::accumulate_hr_for_iat0(ucell, &hr, overlaps, 0, adjs, *paraV, potential);
            double band_trace = 0.0;
            for (int ik = 0; ik < 2; ++ik)
            {
                std::vector<Complex> hk(local_entries, 0.0);
                hamilt::folding_HR(hr, hk.data(), kvec[ik], rows, 1);
                for (int j = 0; j < cols; ++j)
                {
                    const int gj = layout.icglb(j);
                    for (int i = 0; i < rows; ++i)
                    {
                        const int gi = layout.irglb(i);
                        const int index = j * rows + i;
                        Complex expected_h = 0.0;
                        for (int m = 0; m < shell; ++m)
                        {
                            for (int n = 0; n < shell; ++n)
                            {
                                expected_h += std::conj(projectors[gi][m])
                                    * expected_v[m * shell + n] * projectors[gj][n];
                            }
                        }
                        if (ik == 1)
                        {
                            expected_h = std::conj(expected_h);
                        }
                        const Complex residual = hk[index] - expected_h;
                        const double error = std::abs(residual);
                        EXPECT_LT(error, 1e-11);
                        const Complex term = dmk[ik][index] * hk[index];
                        band_trace += term.real();
                    }
                }
            }
            Parallel_Reduce::reduce_all(band_trace);
            const double spin_summed_trace = 2.0 * shell_trace;
            EXPECT_NEAR(band_trace, spin_summed_trace, 1e-10);
        }
    }
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}

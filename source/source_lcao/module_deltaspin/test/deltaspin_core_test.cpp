#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

/***********************************************************************
 * Unit tests for DeltaSpin core algorithms.
 *
 * These tests target the most complex and bug-prone logic:
 * 1. pauli_to_moment: spinor -> magnetic moment conversion
 * 2. calculate_delta_hcc: Pauli matrix expansion for npol=1/2
 * 3. accumulate_Mi_from_becp: multi-k-point magnetic moment accumulation
 * 4. Adaptive threshold calculation
 * 5. Gradient decay check
 * 6. Direction-only lambda projection
 ***********************************************************************/

struct Vec3 { double x, y, z; };
struct Vec3i { int x, y, z; };

// =====================================================================
// 1. pauli_to_moment: spinor -> magnetic moment
//
// Mx = w * (occ[1] + occ[2]).real()
// My = w * (occ[1] - occ[2]).imag()
// Mz = w * (occ[0] - occ[3]).real()
// =====================================================================

static Vec3 pauli_to_moment(const std::complex<double> occ[4], double weight)
{
    return {
        weight * (occ[1] + occ[2]).real(),
        weight * (occ[1] - occ[2]).imag(),
        weight * (occ[0] - occ[3]).real()
    };
}

class PauliToMomentTest : public ::testing::Test
{
  protected:
    std::complex<double> occ[4];
    void SetUp() override { for (int i = 0; i < 4; i++) occ[i] = {0.0, 0.0}; }
};

TEST_F(PauliToMomentTest, PureSpinUp)
{
    occ[0] = {1.0, 0.0}; // |a|^2 = 1
    auto M = pauli_to_moment(occ, 1.0);
    EXPECT_NEAR(M.x, 0.0, 1e-15);
    EXPECT_NEAR(M.y, 0.0, 1e-15);
    EXPECT_NEAR(M.z, 1.0, 1e-15);
}

TEST_F(PauliToMomentTest, PureSpinDown)
{
    occ[3] = {1.0, 0.0}; // |b|^2 = 1
    auto M = pauli_to_moment(occ, 1.0);
    EXPECT_NEAR(M.z, -1.0, 1e-15);
}

TEST_F(PauliToMomentTest, EqualSuperposition_SpinX)
{
    // |psi> = 1/sqrt(2) * (|up> + |down>)
    double a = 1.0 / std::sqrt(2.0);
    occ[0] = {a * a, 0.0};
    occ[1] = {a * a, 0.0}; // a*b
    occ[2] = {a * a, 0.0}; // b*a
    occ[3] = {a * a, 0.0};
    auto M = pauli_to_moment(occ, 1.0);
    // Mx = (0.5 + 0.5) = 1.0, My = 0, Mz = 0
    EXPECT_NEAR(M.x, 1.0, 1e-15);
    EXPECT_NEAR(M.y, 0.0, 1e-15);
    EXPECT_NEAR(M.z, 0.0, 1e-15);
}

TEST_F(PauliToMomentTest, GeneralCase_AllComponents)
{
    occ[0] = {0.6, 0.0};
    occ[1] = {0.1, 0.2};
    occ[2] = {0.1, -0.2}; // conj of occ[1]
    occ[3] = {0.4, 0.0};
    auto M = pauli_to_moment(occ, 1.0);
    // Mx = (0.1+0.2i + 0.1-0.2i).real = 0.2
    // My = (0.1+0.2i - (0.1-0.2i)).imag = (0+0.4i).imag = 0.4
    // Mz = (0.6 - 0.4) = 0.2
    EXPECT_NEAR(M.x, 0.2, 1e-15);
    EXPECT_NEAR(M.y, 0.4, 1e-15);
    EXPECT_NEAR(M.z, 0.2, 1e-15);
}

// =====================================================================
// 2. calculate_delta_hcc: Pauli matrix expansion
//
// npol=2: H += becp^H * lambda * becp
//   lambda in Pauli basis: |lambda_z    lambda_x+i*lambda_y|
//                          |lambda_x-i*lambda_y   -lambda_z |
//
// npol=1: H += becp^H * lambda_z * sign * becp
// =====================================================================

class DeltaHCCTest : public ::testing::Test
{
  protected:
    static void compute_delta_hcc_npol2(
        std::vector<std::complex<double>>& ps,
        const std::vector<std::complex<double>>& becp,
        const Vec3& delta_lambda, int nbands, int nkb, int npol)
    {
        const std::complex<double> c0(delta_lambda.z, 0.0);
        const std::complex<double> c1(delta_lambda.x, delta_lambda.y);
        const std::complex<double> c2(delta_lambda.x, -delta_lambda.y);
        const std::complex<double> c3(-delta_lambda.z, 0.0);
        for (int ib = 0; ib < nbands * npol; ib += npol)
        {
            for (int ip = 0; ip < nkb; ip++)
            {
                const int becpind = ib * nkb + ip;
                const std::complex<double> b1 = becp[becpind];
                const std::complex<double> b2 = becp[becpind + nkb];
                ps[becpind] += c0 * b1 + c2 * b2;
                ps[becpind + nkb] += c1 * b1 + c3 * b2;
            }
        }
    }
    static void compute_delta_hcc_npol1(
        std::vector<std::complex<double>>& ps,
        const std::vector<std::complex<double>>& becp,
        double lambda_z, int sign, int nbands, int nkb)
    {
        double coeff = lambda_z * sign;
        for (int ib = 0; ib < nbands; ib++)
            for (int ip = 0; ip < nkb; ip++)
                ps[ib * nkb + ip] += coeff * becp[ib * nkb + ip];
    }
};

TEST_F(DeltaHCCTest, Npol2_ZeroLambda)
{
    const int nkb = 2, nbands = 1, npol = 2;
    std::vector<std::complex<double>> becp(nbands * npol * nkb, {1.0, 0.0});
    std::vector<std::complex<double>> ps(nbands * npol * nkb, {0.0, 0.0});
    compute_delta_hcc_npol2(ps, becp, {0.0, 0.0, 0.0}, nbands, nkb, npol);
    for (auto& val : ps) {
        EXPECT_NEAR(val.real(), 0.0, 1e-15);
        EXPECT_NEAR(val.imag(), 0.0, 1e-15);
    }
}

TEST_F(DeltaHCCTest, Npol2_PureZ)
{
    const int nkb = 1, nbands = 1, npol = 2;
    std::vector<std::complex<double>> becp(nbands * npol * nkb);
    becp[0] = {1.0, 0.0}; becp[1] = {0.0, 0.0};
    std::vector<std::complex<double>> ps(nbands * npol * nkb, {0.0, 0.0});
    compute_delta_hcc_npol2(ps, becp, {0.0, 0.0, 2.0}, nbands, nkb, npol);
    // c0 = (2,0), c3 = (-2,0); ps_up = 2, ps_dn = 0
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 0.0, 1e-15);
}

TEST_F(DeltaHCCTest, Npol2_PureX)
{
    const int nkb = 1, nbands = 1, npol = 2;
    std::vector<std::complex<double>> becp(nbands * npol * nkb);
    becp[0] = {1.0, 0.0}; becp[1] = {1.0, 0.0};
    std::vector<std::complex<double>> ps(nbands * npol * nkb, {0.0, 0.0});
    compute_delta_hcc_npol2(ps, becp, {3.0, 0.0, 0.0}, nbands, nkb, npol);
    // c1 = c2 = (3,0); ps_up = 3, ps_dn = 3
    EXPECT_NEAR(ps[0].real(), 3.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 3.0, 1e-15);
}

TEST_F(DeltaHCCTest, Npol1_SpinUpPositive)
{
    const int nkb = 2, nbands = 1;
    std::vector<std::complex<double>> becp = {{1.0, 0.5}, {0.0, -0.5}};
    std::vector<std::complex<double>> ps(nkb, {0.0, 0.0});
    compute_delta_hcc_npol1(ps, becp, 2.0, 1, nbands, nkb);
    // ps[0] = 2.0 * (1.0, 0.5) = (2.0, 1.0)
    // ps[1] = 2.0 * (0.0, -0.5) = (0.0, -1.0)
    EXPECT_NEAR(ps[0].real(), 2.0, 1e-15);
    EXPECT_NEAR(ps[0].imag(), 1.0, 1e-15);
    EXPECT_NEAR(ps[1].real(), 0.0, 1e-15);
    EXPECT_NEAR(ps[1].imag(), -1.0, 1e-15);
}

// =====================================================================
// 3. accumulate_Mi_from_becp: multi-k-point accumulation
//
// npol=1: Mi.z += sign[ik] * weight[ik] * sum(|becp|^2)
// npol=2: Mi.xyz from Pauli decomposition of becp
// =====================================================================

class AccumulateMiTest : public ::testing::Test
{
  protected:
    static double accumulate_mi_npol1(
        const std::vector<std::complex<double>>& becp,
        int nbands, int nkb, int nk,
        const double* weights, const int* isk, int nproj)
    {
        double Mi_z = 0.0;
        for (int ik = 0; ik < nk; ik++) {
            int sign = (isk[ik] == 0) ? 1 : -1;
            int offset = ik * nbands * nkb;
            for (int ib = 0; ib < nbands; ib++) {
                double occ = 0.0;
                for (int ip = 0; ip < nproj; ip++) {
                    int idx = offset + ib * nkb + ip;
                    occ += (std::conj(becp[idx]) * becp[idx]).real();
                }
                Mi_z += sign * weights[ik] * occ;
            }
        }
        return Mi_z;
    }
    static Vec3 accumulate_mi_npol2(
        const std::vector<std::complex<double>>& becp,
        int nbands, int nkb, int nk,
        const double* weights, int nproj)
    {
        Vec3 Mi = {0.0, 0.0, 0.0};
        for (int ik = 0; ik < nk; ik++) {
            int offset = ik * nbands * nkb * 2;
            for (int ib = 0; ib < nbands; ib++) {
                std::complex<double> occ[4] = {{0,0},{0,0},{0,0},{0,0}};
                for (int ip = 0; ip < nproj; ip++) {
                    int up_idx = offset + ib * nkb * 2 + ip;
                    int dn_idx = up_idx + nkb;
                    occ[0] += std::conj(becp[up_idx]) * becp[up_idx];
                    occ[1] += std::conj(becp[up_idx]) * becp[dn_idx];
                    occ[2] += std::conj(becp[dn_idx]) * becp[up_idx];
                    occ[3] += std::conj(becp[dn_idx]) * becp[dn_idx];
                }
                Mi.x += weights[ik] * (occ[1] + occ[2]).real();
                Mi.y += weights[ik] * (occ[1] - occ[2]).imag();
                Mi.z += weights[ik] * (occ[0] - occ[3]).real();
            }
        }
        return Mi;
    }
};

TEST_F(AccumulateMiTest, Npol1_SingleK_SpinUp)
{
    const int nbands = 1, nkb = 2, nk = 1, nproj = 2;
    double weights[] = {1.0}; int isk[] = {0};
    std::vector<std::complex<double>> becp = {{0.8, 0.0}, {0.0, 0.6}};
    double Mi_z = accumulate_mi_npol1(becp, nbands, nkb, nk, weights, isk, nproj);
    // occ = 0.64 + 0.36 = 1.0
    EXPECT_NEAR(Mi_z, 1.0, 1e-15);
}

TEST_F(AccumulateMiTest, Npol1_MultiK_MixedSpin)
{
    const int nbands = 1, nkb = 1, nk = 2, nproj = 1;
    double weights[] = {1.0, 0.5}; int isk[] = {0, 1};
    std::vector<std::complex<double>> becp = {{1.0, 0.0}, {1.0, 0.0}};
    double Mi_z = accumulate_mi_npol1(becp, nbands, nkb, nk, weights, isk, nproj);
    // k=0: +1 * 1.0 * 1.0 = 1.0; k=1: -1 * 0.5 * 1.0 = -0.5
    EXPECT_NEAR(Mi_z, 0.5, 1e-15);
}

TEST_F(AccumulateMiTest, Npol2_PureZMag)
{
    const int nbands = 1, nkb = 1, nk = 1, nproj = 1;
    double weights[] = {1.0};
    std::vector<std::complex<double>> becp(nbands * nkb * 2, {0.0, 0.0});
    becp[0] = {0.7, 0.0}; becp[1] = {0.0, 0.0};
    Vec3 Mi = accumulate_mi_npol2(becp, nbands, nkb, nk, weights, nproj);
    // occ[0] = 0.49, occ[3] = 0; Mz = 0.49
    EXPECT_NEAR(Mi.x, 0.0, 1e-15);
    EXPECT_NEAR(Mi.y, 0.0, 1e-15);
    EXPECT_NEAR(Mi.z, 0.49, 1e-15);
}

// =====================================================================
// 4. Adaptive threshold calculation
//
// current_sc_thr = max(initial_rms * sc_drop_thr, sc_thr)
// =====================================================================

class AdaptiveThresholdTest : public ::testing::Test
{
  protected:
    static double calc_adaptive(double initial_rms, double sc_thr, double sc_drop_thr)
    { return std::max(initial_rms * sc_drop_thr, sc_thr); }
};

TEST_F(AdaptiveThresholdTest, DefaultAndTightThreshold)
{
    // Default case: adaptive threshold dominates
    double thr1 = calc_adaptive(0.5, 1e-6, 1e-3);
    EXPECT_NEAR(thr1, 5e-4, 1e-15);
    // Tight threshold: sc_thr dominates
    double thr2 = calc_adaptive(0.1, 1e-4, 1e-3);
    EXPECT_NEAR(thr2, 1e-4, 1e-15);
}

TEST_F(AdaptiveThresholdTest, ZeroScDropThr)
{
    double thr = calc_adaptive(0.5, 1e-6, 0.0);
    EXPECT_NEAR(thr, 1e-6, 1e-15);
}

// =====================================================================
// 5. Gradient decay check
//
// dM/dlambda diagonal: (Mi_new - Mi_old) / (lambda_new - lambda_old)
// If max gradient < decay_grad, stop lambda optimization
// =====================================================================

class GradientDecayTest : public ::testing::Test
{
  protected:
    static bool check_gradient_decay(
        const std::vector<Vec3>& new_spin, const std::vector<Vec3>& old_spin,
        const std::vector<Vec3>& new_lambda, const std::vector<Vec3>& old_lambda,
        const std::vector<Vec3i>& constrain, double decay_grad, int nat)
    {
        for (int iat = 0; iat < nat; iat++) {
            for (int ic = 0; ic < 3; ic++) {
                if (!constrain[iat].x && ic == 0) continue;
                if (!constrain[iat].y && ic == 1) continue;
                if (!constrain[iat].z && ic == 2) continue;
                double dM = 0.0, dL = 0.0;
                if (ic == 0) { dM = new_spin[iat].x - old_spin[iat].x; dL = new_lambda[iat].x - old_lambda[iat].x; }
                if (ic == 1) { dM = new_spin[iat].y - old_spin[iat].y; dL = new_lambda[iat].y - old_lambda[iat].y; }
                if (ic == 2) { dM = new_spin[iat].z - old_spin[iat].z; dL = new_lambda[iat].z - old_lambda[iat].z; }
                if (std::abs(dL) < 1e-30) continue;
                if (std::abs(dM / dL) > decay_grad) return false;
            }
        }
        return true;
    }
};

TEST_F(GradientDecayTest, GradientBelowThreshold)
{
    std::vector<Vec3> new_spin = {{0.0, 0.0, 0.1}};
    std::vector<Vec3> old_spin = {{0.0, 0.0, 0.2}};
    std::vector<Vec3> new_lambda = {{0.0, 0.0, 1.0}};
    std::vector<Vec3> old_lambda = {{0.0, 0.0, 2.0}};
    std::vector<Vec3i> constrain = {{0, 0, 1}};
    // gradient = |0.1 - 0.2| / |1.0 - 2.0| = 0.1
    bool decayed = check_gradient_decay(new_spin, old_spin, new_lambda, old_lambda, constrain, 0.9, 1);
    EXPECT_TRUE(decayed);
}

TEST_F(GradientDecayTest, GradientAboveThreshold)
{
    std::vector<Vec3> new_spin = {{0.0, 0.0, 10.0}};
    std::vector<Vec3> old_spin = {{0.0, 0.0, 0.2}};
    std::vector<Vec3> new_lambda = {{0.0, 0.0, 1.0}};
    std::vector<Vec3> old_lambda = {{0.0, 0.0, 2.0}};
    std::vector<Vec3i> constrain = {{0, 0, 1}};
    // gradient = |10.0 - 0.2| / |1.0 - 2.0| = 9.8
    bool decayed = check_gradient_decay(new_spin, old_spin, new_lambda, old_lambda, constrain, 0.9, 1);
    EXPECT_FALSE(decayed);
}

// =====================================================================
// 6. Direction-only lambda projection
//
// lambda_perp = lambda - (lambda . target_hat) * target_hat
// =====================================================================

static Vec3 project_lambda_perpendicular(const Vec3& lambda, const Vec3& target)
{
    double mag = std::sqrt(target.x * target.x + target.y * target.y + target.z * target.z);
    if (mag < 1e-15) return lambda;
    double tx = target.x / mag, ty = target.y / mag, tz = target.z / mag;
    double dot = lambda.x * tx + lambda.y * ty + lambda.z * tz;
    return { lambda.x - dot * tx, lambda.y - dot * ty, lambda.z - dot * tz };
}

class DirectionOnlyTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(DirectionOnlyTest, ParallelLambda_ZeroResult)
{
    Vec3 lambda = {1.0, 2.0, 3.0};
    Vec3 target = {2.0, 4.0, 6.0};
    auto perp = project_lambda_perpendicular(lambda, target);
    EXPECT_NEAR(perp.x, 0.0, 1e-14);
    EXPECT_NEAR(perp.y, 0.0, 1e-14);
    EXPECT_NEAR(perp.z, 0.0, 1e-14);
}

TEST_F(DirectionOnlyTest, GeneralCase_Perpendicular)
{
    Vec3 lambda = {3.0, 2.0, 1.0};
    Vec3 target = {1.0, 1.0, 1.0};
    auto perp = project_lambda_perpendicular(lambda, target);
    // Verify result is perpendicular to target
    double dot = perp.x * target.x + perp.y * target.y + perp.z * target.z;
    EXPECT_NEAR(dot, 0.0, 1e-14);
}

// =====================================================================
// 7. Constraint energy: E_scon = -sum(lambda . Mi)
// =====================================================================

class ConstraintEnergyTest : public ::testing::Test
{
  protected:
    static double calc_escon(
        const std::vector<Vec3>& lambda, const std::vector<Vec3>& Mi,
        const std::vector<Vec3i>& constrain, int nat)
    {
        double escon = 0.0;
        for (int iat = 0; iat < nat; iat++) {
            if (constrain[iat].x) escon -= lambda[iat].x * Mi[iat].x;
            if (constrain[iat].y) escon -= lambda[iat].y * Mi[iat].y;
            if (constrain[iat].z) escon -= lambda[iat].z * Mi[iat].z;
        }
        return escon;
    }
};

TEST_F(ConstraintEnergyTest, SingleAtom_ZConstrained)
{
    std::vector<Vec3> lambda = {{0.0, 0.0, 1.5}};
    std::vector<Vec3> Mi = {{0.0, 0.0, 2.0}};
    std::vector<Vec3i> constrain = {{0, 0, 1}};
    double escon = calc_escon(lambda, Mi, constrain, 1);
    EXPECT_DOUBLE_EQ(escon, -3.0);
}

TEST_F(ConstraintEnergyTest, MultiAtom_MixedConstraints)
{
    std::vector<Vec3> lambda = {{0.0, 0.0, 1.0}, {0.5, 0.0, 0.0}, {0.0, 0.0, 2.0}};
    std::vector<Vec3> Mi = {{0.0, 0.0, 3.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, -1.0}};
    std::vector<Vec3i> constrain = {{0, 0, 1}, {1, 0, 0}, {0, 0, 1}};
    // atom0: -(1.0 * 3.0) = -3.0
    // atom1: -(0.5 * 1.0) = -0.5
    // atom2: -(2.0 * (-1.0)) = 2.0
    double escon = calc_escon(lambda, Mi, constrain, 3);
    EXPECT_DOUBLE_EQ(escon, -1.5);
}

// =====================================================================
// 8. RMS error calculation
// =====================================================================

static double compute_rms(
    const std::vector<Vec3>& Mi, const std::vector<Vec3>& target,
    const std::vector<Vec3i>& constrain, int nat)
{
    double sum_sq = 0.0; int n_constrained = 0;
    for (int iat = 0; iat < nat; iat++) {
        if (constrain[iat].x) { sum_sq += std::pow(Mi[iat].x - target[iat].x, 2); n_constrained++; }
        if (constrain[iat].y) { sum_sq += std::pow(Mi[iat].y - target[iat].y, 2); n_constrained++; }
        if (constrain[iat].z) { sum_sq += std::pow(Mi[iat].z - target[iat].z, 2); n_constrained++; }
    }
    if (n_constrained == 0) return 0.0;
    return std::sqrt(sum_sq / n_constrained);
}

class RmsErrorTest : public ::testing::Test
{
  protected:
    void SetUp() override {}
};

TEST_F(RmsErrorTest, SingleComponent)
{
    std::vector<Vec3> Mi = {{0.0, 0.0, 1.5}};
    std::vector<Vec3> target = {{0.0, 0.0, 2.0}};
    std::vector<Vec3i> constrain = {{0, 0, 1}};
    double rms = compute_rms(Mi, target, constrain, 1);
    EXPECT_NEAR(rms, 0.5, 1e-15);
}

TEST_F(RmsErrorTest, MultiAtom)
{
    std::vector<Vec3> Mi = {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}};
    std::vector<Vec3> target = {{2.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
    std::vector<Vec3i> constrain = {{1, 0, 0}, {0, 1, 0}};
    // atom0: error_x = 1.0, atom1: error_y = 1.0
    // RMS = sqrt((1+1)/2) = 1.0
    double rms = compute_rms(Mi, target, constrain, 2);
    EXPECT_NEAR(rms, 1.0, 1e-15);
}

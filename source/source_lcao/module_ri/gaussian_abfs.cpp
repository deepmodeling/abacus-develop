//=======================
// AUTHOR : jiyy
// DATE :   2024-02-27
//=======================

#ifndef GAUSSIAN_ABFS_CPP
#define GAUSSIAN_ABFS_CPP

#include "gaussian_abfs.h"

#include <algorithm>
#include <cmath>
#include <numeric>
// #include <chrono>

#include "lri_cv_tools.h"
#include "source_base/global_variable.h"
#include "source_base/math_ylmreal.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
//#include "source_pw/hamilt_pwdft/global.h"

#include <RI/global/Global_Func-1.h>

namespace
{
constexpr int GH_N32 = 32;
constexpr double GH_X32[GH_N32] = {
    -7.125813909830727, -6.409498149269660, -5.812225949515914, -5.275550986515880,
    -4.777164503502596, -4.305547953351198, -3.853755485471445, -3.417167492818570,
    -2.992490825002374, -2.577249537732317, -2.169499183606112, -1.767654109463202,
    -1.370376410952871, -0.976500463589683, -0.584978765435932, -0.194840741569399,
     0.194840741569399,  0.584978765435932,  0.976500463589683,  1.370376410952871,
     1.767654109463202,  2.169499183606112,  2.577249537732317,  2.992490825002374,
     3.417167492818570,  3.853755485471445,  4.305547953351198,  4.777164503502596,
     5.275550986515880,  5.812225949515914,  6.409498149269660,  7.125813909830727};
constexpr double GH_W32[GH_N32] = {
    7.310676427384163e-23, 9.231736536518292e-19, 1.197344017092849e-15, 4.215010211326448e-13,
    5.933291463396639e-11, 4.098832164770897e-09, 1.574167792545595e-07, 3.650585129562376e-06,
    5.416584061819986e-05, 5.362683655279965e-04, 3.654890326654428e-03, 1.755342883157974e-02,
    6.045813095591261e-02, 1.511269427958280e-01, 2.774581423025293e-01, 3.752383525928023e-01,
    3.752383525928023e-01, 2.774581423025293e-01, 1.511269427958280e-01, 6.045813095591261e-02,
    1.755342883157974e-02, 3.654890326654428e-03, 5.362683655279965e-04, 5.416584061819986e-05,
    3.650585129562376e-06, 1.574167792545595e-07, 4.098832164770897e-09, 5.933291463396639e-11,
    4.215010211326448e-13, 1.197344017092849e-15, 9.231736536518292e-19, 7.310676427384163e-23};
}

void Gaussian_Abfs::init(const UnitCell& ucell,
                         const int& Lmax,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                         const ModuleBase::Matrix3& G,
                         const double& lambda)
{
    ModuleBase::TITLE("Gaussian_Abfs", "init");
    ModuleBase::timer::start("Gaussian_Abfs", "init");

    this->kvec_c = kvec_c;
    const int nks0 = kvec_c.size();

    this->lambda = lambda;
    this->tpiba = ucell.tpiba;
    this->lat0 = ucell.lat0;
    this->omega = ucell.omega;
    const ModuleBase::Vector3<double> a1_bohr = ucell.a1 * this->lat0;
    const ModuleBase::Vector3<double> a2_bohr = ucell.a2 * this->lat0;
    this->area_parallel_bohr2 = (a1_bohr ^ a2_bohr).norm();
    this->lz_bohr = this->omega / this->area_parallel_bohr2;
    const double eta = 35;
    std::vector<ModuleBase::Vector3<double>> Gvec;
    Gvec.resize(3);
    Gvec[0].x = G.e11;
    Gvec[0].y = G.e12;
    Gvec[0].z = G.e13;

    Gvec[1].x = G.e21;
    Gvec[1].y = G.e22;
    Gvec[1].z = G.e23;

    Gvec[2].x = G.e31;
    Gvec[2].y = G.e32;
    Gvec[2].z = G.e33;
    this->gvecs_direct = Gvec;

    this->n_cells.resize(nks0);
    this->n_supercells_2d.resize(nks0);
    this->qGvecs.resize(nks0);
    this->check_gamma.resize(nks0);
    this->ylm.resize(nks0);
    const int total_lm = (Lmax + 1) * (Lmax + 1);

#pragma omp parallel for schedule(dynamic)
    for (size_t ik = 0; ik != nks0; ++ik)
    {
        ModuleBase::Vector3<double> qvec = this->kvec_c[ik];
        const double Gmax = std::sqrt(eta * this->lambda) + qvec.norm() * this->tpiba;
        std::vector<int> n_supercells = get_n_supercells(this->lat0, G, Gmax);
        int total_cells = std::accumulate(n_supercells.begin(), n_supercells.end(), 1, [](int a, int b) {
            return a * (2 * b + 1);
        });

        std::vector<ModuleBase::Vector3<double>> qGvec_ik(total_cells);
        std::vector<bool> check_gamma_ik(total_cells);
        for (int idx = 0; idx < total_cells; ++idx)
        {
            int G0 = (idx / ((2 * n_supercells[1] + 1) * (2 * n_supercells[2] + 1))) - n_supercells[0];
            int G1 = ((idx / (2 * n_supercells[2] + 1)) % (2 * n_supercells[1] + 1)) - n_supercells[1];
            int G2 = (idx % (2 * n_supercells[2] + 1)) - n_supercells[2];
            ModuleBase::Vector3<double> qGvec
                = -(qvec + Gvec[0] * static_cast<double>(G0) + Gvec[1] * static_cast<double>(G1)
                    + Gvec[2] * static_cast<double>(G2));
            qGvec_ik[idx] = qGvec;
            if (G0 == 0 && G1 == 0 && G2 == 0)
                check_gamma_ik[idx] = true;
            else
                check_gamma_ik[idx] = false;
        }
        ModuleBase::matrix ylm_ik(total_lm, total_cells);
        ModuleBase::YlmReal::Ylm_Real(total_lm, total_cells, qGvec_ik.data(), ylm_ik);

#pragma omp critical(Gaussian_Abfs_init)
        {
            this->n_cells[ik] = total_cells;
            this->n_supercells_2d[ik] = {n_supercells[0], n_supercells[1]};
            this->qGvecs[ik] = qGvec_ik;
            this->check_gamma[ik] = check_gamma_ik;
            this->ylm[ik] = ylm_ik;
        }
    }

    ModuleBase::timer::end("Gaussian_Abfs", "init");
}

auto Gaussian_Abfs::get_Vq(const int& lp_max,
                           const int& lq_max, // Maximum L for which to calculate interaction.
                           const size_t& ik,
                           const double& chi, // Singularity corrected value at q=0.
                           const ModuleBase::Vector3<double>& tau,
                           const ModuleBase::realArray& gaunt) -> RI::Tensor<std::complex<double>>
{
    ModuleBase::TITLE("Gaussian_Abfs", "get_Vq");
    ModuleBase::timer::start("Gaussian_Abfs", "get_Vq");

    const T_func_DPcal_lattice_sum<std::complex<double>> func_DPcal_lattice_sum
        = std::bind(&Gaussian_Abfs::get_lattice_sum,
                    this,
                    this->tpiba,
                    ik,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3,
                    std::placeholders::_4,
                    tau);
    auto res = this->DPcal_Vq_dVq<RI::Tensor<std::complex<double>>>(this->omega,
                                                                    lp_max,
                                                                    lq_max,
                                                                    ik,
                                                                    chi,
                                                                    tau,
                                                                    gaunt,
                                                                    func_DPcal_lattice_sum);

    ModuleBase::timer::end("Gaussian_Abfs", "get_Vq");
    return res;
}

auto Gaussian_Abfs::get_Vq_2d(const int& lp_max,
                              const int& lq_max,
                              const size_t& ik,
                              const double& chi2d,
                              const ModuleBase::Vector3<double>& tau,
                              const ModuleBase::realArray& gaunt) -> RI::Tensor<std::complex<double>>
{
    ModuleBase::TITLE("Gaussian_Abfs", "get_Vq_2d");
    ModuleBase::timer::tick("Gaussian_Abfs", "get_Vq_2d");

    auto res = this->DPcal_Vq_2d(lp_max, lq_max, ik, chi2d, tau, gaunt);

    ModuleBase::timer::tick("Gaussian_Abfs", "get_Vq_2d");
    return res;
}

auto Gaussian_Abfs::get_dVq(const int& lp_max,
                            const int& lq_max, // Maximum L for which to calculate interaction.
                            const size_t& ik,
                            const double& chi, // Singularity corrected value at q=0.
                            const ModuleBase::Vector3<double>& tau,
                            const ModuleBase::realArray& gaunt) -> std::array<RI::Tensor<std::complex<double>>, 3>
{
    ModuleBase::TITLE("Gaussian_Abfs", "get_dVq");
    ModuleBase::timer::start("Gaussian_Abfs", "get_dVq");

    const T_func_DPcal_lattice_sum<std::array<std::complex<double>, 3>> func_DPcal_d_lattice_sum
        = std::bind(&Gaussian_Abfs::get_d_lattice_sum,
                    this,
                    this->tpiba,
                    ik,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3,
                    std::placeholders::_4,
                    tau);
    auto res = this->DPcal_Vq_dVq<std::array<RI::Tensor<std::complex<double>>, 3>>(this->omega,
                                                                                   lp_max,
                                                                                   lq_max,
                                                                                   ik,
                                                                                   chi,
                                                                                   tau,
                                                                                   gaunt,
                                                                                   func_DPcal_d_lattice_sum);

    ModuleBase::timer::end("Gaussian_Abfs", "get_dVq");
    return res;
}

template <typename Tout, typename Tin>
auto Gaussian_Abfs::DPcal_Vq_dVq(const double& omega,
                                 const int& lp_max,
                                 const int& lq_max, // Maximum L for which to calculate interaction.
                                 const size_t& ik,
                                 const double& chi, // Singularity corrected value at q=0.
                                 const ModuleBase::Vector3<double>& tau,
                                 const ModuleBase::realArray& gaunt,
                                 const T_func_DPcal_lattice_sum<Tin>& func_DPcal_lattice_sum) -> Tout
{
    const int Lmax = lp_max + lq_max;
    const int n_LM = (Lmax + 1) * (Lmax + 1);
    const size_t vq_ndim0 = (lp_max + 1) * (lp_max + 1);
    const size_t vq_ndim1 = (lq_max + 1) * (lq_max + 1);
    Tout Vq_dVq;
    LRI_CV_Tools::init_elem(Vq_dVq, vq_ndim0, vq_ndim1);
    /*
     n_add_ksq * 2 = lp_max + lq_max - abs(lp_max - lq_max)
        if lp_max < lq_max
            n_add_ksq * 2 = lp_max + lq_max - (lq_max - lp_max)
                          = lp_max * 2
        if lp_max > lq_max
            n_add_ksq * 2 = lp_max + lq_max - (lp_max - lq_max)
                          = lq_max * 2
        thus,
            n_add_ksq = min(lp_max, lq_max)
    */
    const int n_add_ksq = std::min(lp_max, lq_max);
    std::vector<std::vector<Tin>> lattice_sum;
    lattice_sum.resize(n_add_ksq + 1);

    const double exponent = 1.0 / this->lambda;
    ModuleBase::Vector3<double> qvec = this->kvec_c[ik];

    for (int i_add_ksq = 0; i_add_ksq != n_add_ksq + 1; ++i_add_ksq) // integrate lp, lq, L to one index i_add_ksq, i.e.
                                                                     // (lp+lq-L)/2
    {
        const double power = -2.0 + 2 * i_add_ksq;
        const int this_Lmax = Lmax - 2 * i_add_ksq;                         // calculate Lmax at current lp+lq
        const bool exclude_Gamma = (qvec.norm() < 1e-10 && i_add_ksq == 0); // only Gamma point and lq+lp-2>0 need to be
                                                                            // corrected
        lattice_sum[i_add_ksq] = func_DPcal_lattice_sum(power, exponent, exclude_Gamma, this_Lmax);
    }

    /* The exponent term comes in from Taylor expanding the
        Gaussian at zero to first order in k^2, which cancels the k^-2 from the
        Coulomb interaction.  While terms of this order are in principle
        neglected, we make one exception here.  Without this, the final result
        would (slightly) depend on the Ewald lambda.*/
    if (qvec.norm() < 1e-10)
    {
        std::complex<double> val = chi - exponent;
        std::complex<double> frac = 1.0 / std::sqrt(ModuleBase::FOUR_PI);
        LRI_CV_Tools::add_elem(lattice_sum[0][0], val, frac);
    }

    for (int lp = 0; lp != lp_max + 1; ++lp)
    {
        double norm_1 = double_factorial(2 * lp - 1) * std::sqrt(ModuleBase::PI * 0.5);
        for (int lq = 0; lq != lq_max + 1; ++lq)
        {
            double norm_2 = double_factorial(2 * lq - 1) * std::sqrt(ModuleBase::PI * 0.5);
            std::complex<double> phase = std::pow(ModuleBase::IMAG_UNIT, lp - lq);
            std::complex<double> cfac
                = ModuleBase::FOUR_PI * phase * std::pow(ModuleBase::TWO_PI, 3) / (norm_1 * norm_2) / omega;
            for (int L = std::abs(lp - lq); L <= lp + lq; L += 2) // if lp+lq-L == odd, then Gaunt_Coefficients = 0
            {
                const int i_add_ksq = (lp + lq - L) / 2;
                for (int mp = 0; mp != 2 * lp + 1; ++mp)
                {
                    const int lmp = lp * lp + mp;
                    for (int mq = 0; mq != 2 * lq + 1; ++mq)
                    {
                        const int lmq = lq * lq + mq;
                        for (int m = 0; m != 2 * L + 1; ++m)
                        {
                            const int lm = L * L + m;
                            double triple_Y = gaunt(lmp, lmq, lm);
                            std::complex<double> fac = triple_Y * cfac;
                            LRI_CV_Tools::add_elem(Vq_dVq, lmp, lmq, lattice_sum[i_add_ksq][lm], fac);
                        }
                    }
                }
            }
        }
    }

    return Vq_dVq;
}

auto Gaussian_Abfs::DPcal_Vq_2d(const int& lp_max,
                                const int& lq_max,
                                const size_t& ik,
                                const double& chi2d,
                                const ModuleBase::Vector3<double>& tau,
                                const ModuleBase::realArray& gaunt) -> RI::Tensor<std::complex<double>>
{
    const int Lmax = lp_max + lq_max;
    const size_t vq_ndim0 = (lp_max + 1) * (lp_max + 1);
    const size_t vq_ndim1 = (lq_max + 1) * (lq_max + 1);
    RI::Tensor<std::complex<double>> Vq;
    LRI_CV_Tools::init_elem(Vq, vq_ndim0, vq_ndim1);

    const int n_add_ksq = std::min(lp_max, lq_max);
    std::vector<std::vector<std::complex<double>>> lattice_sum(n_add_ksq + 1);
    const double exponent = 1.0 / this->lambda;
    const ModuleBase::Vector3<double> qvec = this->kvec_c[ik];

    for (int i_add_ksq = 0; i_add_ksq != n_add_ksq + 1; ++i_add_ksq)
    {
        const double power = -2.0 + 2 * i_add_ksq;
        const int this_Lmax = Lmax - 2 * i_add_ksq;
        const bool exclude_Gamma = (qvec.norm() < 1e-10 && i_add_ksq == 0);
        lattice_sum[i_add_ksq]
            = this->get_lattice_sum_2d(this->tpiba, ik, power, exponent, exclude_Gamma, this_Lmax, tau);
    }

    if (qvec.norm() < 1e-10)
    {
        const double tau_z_cart = tau.z * this->lat0;
        const double regular = this->C_reg_2d(tau_z_cart, exponent);
        const std::complex<double> val = ModuleBase::PI * chi2d + regular;
        const std::complex<double> frac = 1.0 / std::sqrt(ModuleBase::FOUR_PI);
        LRI_CV_Tools::add_elem(lattice_sum[0][0], val, frac * this->lz_bohr / ModuleBase::TWO_PI);
    }

    for (int lp = 0; lp != lp_max + 1; ++lp)
    {
        const double norm_1 = double_factorial(2 * lp - 1) * std::sqrt(ModuleBase::PI * 0.5);
        for (int lq = 0; lq != lq_max + 1; ++lq)
        {
            const double norm_2 = double_factorial(2 * lq - 1) * std::sqrt(ModuleBase::PI * 0.5);
            const std::complex<double> phase = std::pow(ModuleBase::IMAG_UNIT, lp - lq);
            const std::complex<double> cfac
                = ModuleBase::FOUR_PI * phase * std::pow(ModuleBase::TWO_PI, 3) / (norm_1 * norm_2) / this->omega;
            for (int L = std::abs(lp - lq); L <= lp + lq; L += 2)
            {
                const int i_add_ksq = (lp + lq - L) / 2;
                for (int mp = 0; mp != 2 * lp + 1; ++mp)
                {
                    const int lmp = lp * lp + mp;
                    for (int mq = 0; mq != 2 * lq + 1; ++mq)
                    {
                        const int lmq = lq * lq + mq;
                        for (int m = 0; m != 2 * L + 1; ++m)
                        {
                            const int lm = L * L + m;
                            const double triple_Y = gaunt(lmp, lmq, lm);
                            const std::complex<double> fac = triple_Y * cfac;
                            LRI_CV_Tools::add_elem(Vq, lmp, lmq, lattice_sum[i_add_ksq][lm], fac);
                        }
                    }
                }
            }
        }
    }

    return Vq;
}

Numerical_Orbital_Lm Gaussian_Abfs::Gauss(const Numerical_Orbital_Lm& orb, const double& lambda)
{
    Numerical_Orbital_Lm gaussian;
    const int angular_momentum_l = orb.getL();
    const double eta = 35;
    const double rcut = std::sqrt(eta / lambda);
    const double dr = orb.get_rab().back();
    int Nr = std::ceil(rcut / dr);
    if (Nr % 2 == 0)
        Nr += 1;

    std::vector<double> rab(Nr);
    for (size_t ir = 0; ir < Nr; ++ir)
        rab[ir] = dr;
    std::vector<double> r_radial(Nr);
    for (size_t ir = 0; ir < Nr; ++ir)
        r_radial[ir] = ir * dr;

    const double frac = std::pow(lambda, angular_momentum_l + 1.5) / double_factorial(2 * angular_momentum_l - 1)
                        / std::sqrt(ModuleBase::PI * 0.5);

    std::vector<double> psi(Nr);

    for (size_t ir = 0; ir != Nr; ++ir)
        psi[ir]
            = frac * std::pow(r_radial[ir], angular_momentum_l) * std::exp(-lambda * r_radial[ir] * r_radial[ir] * 0.5);

    gaussian.set_orbital_info(orb.getLabel(),
                              orb.getType(),
                              angular_momentum_l,
                              orb.getChi(),
                              Nr,
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial),
                              Numerical_Orbital_Lm::Psi_Type::Psi,
                              ModuleBase::GlobalFunc::VECTOR_TO_PTR(psi),
                              orb.getNk(),
                              orb.getDk(),
                              orb.getDruniform(),
                              false,
                              true,
                              PARAM.inp.cal_force);

    return gaussian;
}

double Gaussian_Abfs::double_factorial(const int& n)
{
    double result = 1.0;
    for (int i = n; i > 0; i -= 2)
    {
        if (i == 1)
            result *= 1.0;
        else
            result *= static_cast<double>(i);
    }
    return result;
}

auto Gaussian_Abfs::get_lattice_sum(const double& tpiba,
                                    const size_t& ik,
                                    const double& power, // Will be 0. for straight GTOs and -2. for Coulomb interaction
                                    const double& exponent,
                                    const bool& exclude_Gamma, // The R==0. can be excluded by this flag.
                                    const int& lmax,           // Maximum angular momentum the sum is needed for.
                                    const ModuleBase::Vector3<double>& tau) -> std::vector<std::complex<double>>
{
    const T_func_DPcal_phase<std::complex<double>> func_DPcal_phase
        = [&tau](const ModuleBase::Vector3<double>& vec) -> std::complex<double> {
        return std::exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (vec * tau));
    };

    return this
        ->DPcal_lattice_sum<std::complex<double>>(tpiba, ik, power, exponent, exclude_Gamma, lmax, func_DPcal_phase);
}

auto Gaussian_Abfs::get_lattice_sum_2d(const double& tpiba,
                                       const size_t& ik,
                                       const double& power,
                                       const double& exponent,
                                       const bool& exclude_Gamma,
                                       const int& lmax,
                                       const ModuleBase::Vector3<double>& tau) -> std::vector<std::complex<double>>
{
    if (power < 0.0 && !exclude_Gamma && this->kvec_c[ik].norm() < 1e-10)
    {
        ModuleBase::WARNING_QUIT("Gaussian_Abfs::lattice_sum_2d", "Gamma point for power<0.0 cannot be evaluated!");
    }

    const int total_lm = (lmax + 1) * (lmax + 1);
    std::vector<std::complex<double>> result(total_lm, std::complex<double>{});
    const int n_super0 = this->n_supercells_2d[ik][0];
    const int n_super1 = this->n_supercells_2d[ik][1];
    const double conversion = this->lz_bohr / ModuleBase::TWO_PI;
    const double tau_z_bohr = tau.z * this->lat0;

    for (int i0 = -n_super0; i0 <= n_super0; ++i0)
    {
        for (int i1 = -n_super1; i1 <= n_super1; ++i1)
        {
            const bool excluded_gamma_term = exclude_Gamma && i0 == 0 && i1 == 0;

            const ModuleBase::Vector3<double> qg_direct
                = -(this->kvec_c[ik] + this->gvecs_direct[0] * static_cast<double>(i0)
                    + this->gvecs_direct[1] * static_cast<double>(i1));
            const ModuleBase::Vector3<double> qg_cart = qg_direct * tpiba;
            const double q_parallel_sq = qg_cart.x * qg_cart.x + qg_cart.y * qg_cart.y;
            const std::complex<double> phase
                = std::exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (qg_direct * tau));
            const double gaussian_parallel = std::exp(-exponent * q_parallel_sq);

            for (int L = 0; L != lmax + 1; ++L)
            {
                // Only L=0 is singular and L=1 has a direction-dependent q->0
                // limit. The L>=2 Coulomb channels have finite Gamma limits.
                if (excluded_gamma_term && (std::abs(power + 2.0) >= 1e-12 || L < 2))
                {
                    continue;
                }
                for (int m = 0; m != 2 * L + 1; ++m)
                {
                    const int lm = L * L + m;
                    const std::complex<double> k_lm = this->K_LM_2d(
                        L,
                        m,
                        power,
                        ModuleBase::Vector3<double>(qg_cart.x, qg_cart.y, 0.0),
                        tau_z_bohr,
                        exponent);
                    result[lm] += conversion * gaussian_parallel * phase * k_lm;
                }
            }
        }
    }

    return result;
}

double Gaussian_Abfs::I2_2d(const double& q_parallel_abs, const double& tau_z, const double& beta) const
{
    if (q_parallel_abs < 1e-12)
    {
        return 0.0;
    }
    if (std::abs(tau_z) < 1e-14)
    {
        const double x = std::sqrt(beta) * q_parallel_abs;
        return ModuleBase::PI / q_parallel_abs * std::exp(x * x) * std::erfc(x);
    }

    const double a = std::abs(tau_z);
    const double q = q_parallel_abs;
    const double sqrt_beta = std::sqrt(beta);
    const double term_plus = std::exp(q * a) * std::erfc(sqrt_beta * q + a / (2.0 * sqrt_beta));
    const double term_minus = std::exp(-q * a) * std::erfc(sqrt_beta * q - a / (2.0 * sqrt_beta));
    return ModuleBase::PI / (2.0 * q) * std::exp(beta * q * q) * (term_plus + term_minus);
}

double Gaussian_Abfs::C_reg_2d(const double& tau_z, const double& beta) const
{
    const double a = std::abs(tau_z);
    const double sqrt_beta = std::sqrt(beta);
    return ModuleBase::PI * a * (std::erfc(a / (2.0 * sqrt_beta)) - 1.0)
           - 2.0 * std::sqrt(ModuleBase::PI * beta) * std::exp(-a * a / (4.0 * beta));
}

std::complex<double> Gaussian_Abfs::K_LM_2d(const int& L,
                                            const int& M,
                                            const double& power,
                                            const ModuleBase::Vector3<double>& q_parallel_cart,
                                            const double& tau_z,
                                            const double& beta) const
{
    const double q_parallel_abs
        = std::sqrt(q_parallel_cart.x * q_parallel_cart.x + q_parallel_cart.y * q_parallel_cart.y);
    if (L == 0 && M == 0 && std::abs(power + 2.0) < 1e-12)
    {
        return this->I2_2d(q_parallel_abs, tau_z, beta) / std::sqrt(ModuleBase::FOUR_PI);
    }

    const int total_lm = (L + 1) * (L + 1);
    const int lm = L * L + M;
    const double sqrt_beta = std::sqrt(beta);
    const double inv_sqrt_beta = 1.0 / sqrt_beta;
    std::complex<double> sum = 0.0;

    for (int i = 0; i != GH_N32; ++i)
    {
        const double kz = GH_X32[i] * inv_sqrt_beta;
        const double norm_sq = q_parallel_abs * q_parallel_abs + kz * kz;
        if (norm_sq < 1e-24 && power < 0.0)
        {
            continue;
        }
        ModuleBase::Vector3<double> q3(q_parallel_cart.x, q_parallel_cart.y, kz);
        ModuleBase::matrix ylm_here(total_lm, 1);
        ModuleBase::YlmReal::Ylm_Real(total_lm, 1, &q3, ylm_here);
        const double radial = std::pow(norm_sq, 0.5 * (power + L));
        const std::complex<double> phase = std::exp(ModuleBase::IMAG_UNIT * kz * tau_z);
        sum += GH_W32[i] * radial * ylm_here(lm, 0) * phase;
    }

    return inv_sqrt_beta * sum;
}

auto Gaussian_Abfs::get_d_lattice_sum(
    const double& tpiba,
    const size_t& ik,
    const double& power, // Will be 0. for straight GTOs and -2. for Coulomb interaction
    const double& exponent,
    const bool& exclude_Gamma, // The R==0. can be excluded by this flag.
    const int& lmax,           // Maximum angular momentum the sum is needed for.
    const ModuleBase::Vector3<double>& tau) -> std::vector<std::array<std::complex<double>, 3>>
{
    const T_func_DPcal_phase<std::array<std::complex<double>, 3>> func_DPcal_d_phase
        = [&tau, &tpiba](const ModuleBase::Vector3<double>& vec) -> std::array<std::complex<double>, 3> {
        using namespace RI::Array_Operator;
        std::complex<double> phase = std::exp(ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (vec * tau));
        std::array<std::complex<double>, 3> ip_vec = {phase * vec.x, phase * vec.y, phase * vec.z};
        std::array<std::complex<double>, 3> d_phase = tpiba * ModuleBase::IMAG_UNIT * ip_vec;

        return d_phase;
    };

    return this->DPcal_lattice_sum<std::array<std::complex<double>, 3>>(tpiba,
                                                                        ik,
                                                                        power,
                                                                        exponent,
                                                                        exclude_Gamma,
                                                                        lmax,
                                                                        func_DPcal_d_phase);
}

template <typename Tresult>
auto Gaussian_Abfs::DPcal_lattice_sum(
    const double& tpiba,
    const size_t& ik,
    const double& power, // Will be 0. for straight GTOs and -2. for Coulomb interaction
    const double& exponent,
    const bool& exclude_Gamma, // The R==0. can be excluded by this flag.
    const int& lmax,           // Maximum angular momentum the sum is needed for.
    const T_func_DPcal_phase<Tresult>& func_DPcal_phase) -> std::vector<Tresult>
{
    if (power < 0.0 && !exclude_Gamma && this->kvec_c[ik].norm() < 1e-10)
        ModuleBase::WARNING_QUIT("Gaussian_Abfs::lattice_sum", "Gamma point for power<0.0 cannot be evaluated!");

    using namespace RI::Array_Operator;

    const int total_lm = (lmax + 1) * (lmax + 1);
    std::vector<Tresult> result(total_lm, Tresult{});
    const int total_cells = this->n_cells[ik];

#pragma omp declare reduction(vec_plus : std::vector<Tresult> : std::transform(omp_out.begin(),                        \
                                                                                   omp_out.end(),                      \
                                                                                   omp_in.begin(),                     \
                                                                                   omp_out.begin(),                    \
                                                                                   LRI_CV_Tools::plus<Tresult>()))     \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
    auto accumulate_cell = [&](const int idx, std::vector<Tresult>& out) {
        if (exclude_Gamma && this->check_gamma[ik][idx])
            return;

        ModuleBase::Vector3<double> vec = this->qGvecs[ik][idx];
        const double vec_sq = vec.norm2() * tpiba * tpiba;
        const double vec_abs = std::sqrt(vec_sq);

        const double val_s = std::exp(-exponent * vec_sq) * std::pow(vec_abs, power);

        Tresult phase = func_DPcal_phase(vec);
        for (int L = 0; L != lmax + 1; ++L)
        {
            const double val_l = val_s * std::pow(vec_abs, L);
            for (int m = 0; m != 2 * L + 1; ++m)
            {
                const int lm = L * L + m;
                const double val_lm = val_l * this->ylm[ik](lm, idx);
                out[lm] = out[lm] + RI::Global_Func::convert<std::complex<double>>(val_lm) * phase;
            }
        }
    };

//     // auto start0 = std::chrono::system_clock::now();
#pragma omp parallel for reduction(vec_plus : result)
    for (int idx = 0; idx < total_cells; ++idx)
    {
        accumulate_cell(idx, result);
    }

    // auto end0 = std::chrono::system_clock::now();
    // auto duration0 =
    // std::chrono::duration_cast<std::chrono::microseconds>(end0 - start0);
    // std::cout << "lattice Time: "
    //           << double(duration0.count()) *
    //           std::chrono::microseconds::period::num
    //                  / std::chrono::microseconds::period::den
    //           << " s" << std::endl;

    return result;
}

std::vector<int> Gaussian_Abfs::get_n_supercells(const double& lat0, const ModuleBase::Matrix3& G, const double& Gmax)
{
    std::vector<int> n_supercells(3);
    ModuleBase::Matrix3 GI = G.Inverse();
    ModuleBase::Matrix3 latvec = GI.Transpose();
    latvec *= lat0;
    std::vector<ModuleBase::Vector3<double>> lat;
    lat.resize(3);
    lat[0].x = latvec.e11;
    lat[0].y = latvec.e12;
    lat[0].z = latvec.e13;
    lat[1].x = latvec.e21;
    lat[1].y = latvec.e22;
    lat[1].z = latvec.e23;
    lat[2].x = latvec.e31;
    lat[2].y = latvec.e32;
    lat[2].z = latvec.e33;

    n_supercells[0] = static_cast<int>(std::floor(lat[0].norm() * Gmax / ModuleBase::TWO_PI + 1e-5));
    n_supercells[1] = static_cast<int>(std::floor(lat[1].norm() * Gmax / ModuleBase::TWO_PI + 1e-5));
    n_supercells[2] = static_cast<int>(std::floor(lat[2].norm() * Gmax / ModuleBase::TWO_PI + 1e-5));

    return n_supercells;
}

#endif

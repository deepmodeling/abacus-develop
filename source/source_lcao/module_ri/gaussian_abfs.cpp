//=======================
// AUTHOR : jiyy
// DATE :   2024-02-27
//=======================

#ifndef GAUSSIAN_ABFS_CPP
#define GAUSSIAN_ABFS_CPP

#include "gaussian_abfs.h"

#include <algorithm>
// #include <chrono>

#include "LRI_CV_Tools.h"
#include "source_base/global_variable.h"
#include "source_base/math_ylmreal.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
//#include "source_pw/hamilt_pwdft/global.h"

#include <RI/global/Global_Func-1.h>

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

    this->n_cells.resize(nks0);
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
//     // auto start0 = std::chrono::system_clock::now();
#pragma omp parallel for reduction(vec_plus : result)
    for (int idx = 0; idx < total_cells; ++idx)
    {
        if (exclude_Gamma && this->check_gamma[ik][idx])
            continue;

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
                result[lm] = result[lm] + RI::Global_Func::convert<std::complex<double>>(val_lm) * phase;
            }
        }
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
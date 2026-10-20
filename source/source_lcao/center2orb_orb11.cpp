//=========================================================
// AUTHOR : Peize Lin
// DATE : 2016-01-24
//=========================================================

#include "center2orb_orb11.h"

#include "center2orb.h"
#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "source_base/mathzone_add1.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/timer.h"

#ifdef __CUDA
#include "source_base/module_device/memory_op.h"

#include <ATen/core/tensor.h>
#include <ATen/kernels/blas.h>
#endif
#include "source_base/math_polyint.h"
#include "source_base/sph_bessel_recursive.h"
#include "source_base/ylm.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <tuple>

Center2_Orb::Orb11::Orb11(const Numerical_Orbital_Lm& nA_in,
                          const Numerical_Orbital_Lm& nB_in,
                          const ModuleBase::Sph_Bessel_Recursive::D2* psb,
                          const ORB_gaunt_table& MGT_in)
    : nA(nA_in), nB(nB_in), psb_(psb), MGT(MGT_in)
{
}

void Center2_Orb::Orb11::init_radial_table()
{
    const int rmesh = Center2_Orb::get_rmesh(this->nA.getRcut(), this->nB.getRcut(), dr_);

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    for (int LAB = std::abs(LA - LB); LAB <= LA + LB; ++LAB)
    {
        if ((LAB - std::abs(LA - LB)) % 2 == 1) // if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
        {
            continue;
        }

        this->Table_r[LAB].resize(rmesh, 0);
        this->Table_dr[LAB].resize(rmesh, 0);

        Center2_Orb::cal_ST_Phi12_R(1,
                                    LAB,
                                    this->nA,
                                    this->nB,
                                    rmesh,
                                    this->Table_r[LAB],
                                    this->Table_dr[LAB],
                                    psb_);
    }
    return;
}

void Center2_Orb::Orb11::init_radial_table(const std::set<size_t>& radials)
{
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();

    const size_t rmesh = Center2_Orb::get_rmesh(this->nA.getRcut(), this->nB.getRcut(), dr_);

    std::set<size_t> radials_used;
    for (const size_t& ir: radials)
    {
        if (ir < rmesh)
        {
            radials_used.insert(ir);
        }
    }

    for (int LAB = std::abs(LA - LB); LAB <= LA + LB; ++LAB)
    {
        if ((LAB - std::abs(LA - LB)) % 2 == 1) // if LA+LB-LAB == odd, then Gaunt_Coefficients = 0
        { 
            continue;
        }

        this->Table_r[LAB].resize(rmesh, 0);
        this->Table_dr[LAB].resize(rmesh, 0);

        Center2_Orb::cal_ST_Phi12_R(1,
                                    LAB,
                                    this->nA,
                                    this->nB,
                                    radials_used,
                                    this->Table_r[LAB],
                                    this->Table_dr[LAB],
                                    psb_);
    }
}

namespace
{
double radial_product(const Numerical_Orbital_Lm& left, const Numerical_Orbital_Lm& right, const int ik)
{
    if (!left.get_psif().empty() && !right.get_psi_k2().empty())
    {
        return left.getPsif(ik) * right.getPsi_k2(ik);
    }
    if (!left.get_psi_k().empty() && !right.get_psi_k().empty())
    {
        return left.getPsi_k(ik) * right.getPsi_k(ik);
    }
    if (!left.get_psi_k2().empty() && !right.get_psif().empty())
    {
        return left.getPsi_k2(ik) * right.getPsif(ik);
    }
    return 0.0;
}
} // namespace

void Center2_Orb::Orb11::init_radial_table(const std::vector<Orb11*>& orbitals, const bool use_gpu)
{
    ModuleBase::timer::start("Center2_Orb", "init_radial_table");
#ifndef __CUDA
    if (use_gpu)
    {
        throw std::invalid_argument("GPU radial tables require a CUDA build");
    }
#else
    if (use_gpu)
    {
        container::kernels::createGpuBlasHandle();
    }
#endif
    // Group by quadrature grid and angular momentum. All orbitals in a batch
    // belong to one Matrix_Orbs object and share its immutable Bessel table.
    using Grid = std::tuple<int, int, int, double>;
    std::map<Grid, std::vector<Orb11*>> groups;
    for (Orb11* orbital: orbitals)
    {
        if (orbital->psb_ != orbitals.front()->psb_)
        {
            throw std::invalid_argument("A radial batch must share its Bessel table");
        }
        const int rmesh = Center2_Orb::get_rmesh(orbital->nA.getRcut(), orbital->nB.getRcut(), orbital->dr_);
        const int kmesh = orbital->nA.getNk();
        const double dk = orbital->nA.getDk();
        const int first_l = std::abs(orbital->nA.getL() - orbital->nB.getL());
        const int last_l = orbital->nA.getL() + orbital->nB.getL();
        for (int l = first_l; l <= last_l; l += 2)
        {
            orbital->Table_r[l].resize(rmesh);
            orbital->Table_dr[l].resize(rmesh);
            const Grid grid(l, kmesh, rmesh, dk);
            groups[grid].push_back(orbital);
        }
    }
    for (const auto& group: groups)
    {
        const int l = std::get<0>(group.first);
        const int kmesh = std::get<1>(group.first);
        const int rmesh = std::get<2>(group.first);
        const double dk = std::get<3>(group.first);
        const auto& tasks = group.second;
        const auto& bessel = tasks.front()->psb_->get_jlx();
        const double* kpoint = tasks.front()->nA.getKpoint();
        // Bound all three device buffers together: at most 160 MiB, released
        // before returning to the SCF/LibRI contraction stages.
        constexpr std::size_t input_limit = 16 * 1024 * 1024 / sizeof(double);
        constexpr std::size_t table_limit = 128 * 1024 * 1024 / sizeof(double);
        constexpr std::size_t output_limit = 16 * 1024 * 1024 / sizeof(double);
        if (kmesh < 3 || kmesh % 2 == 0 || static_cast<std::size_t>(kmesh) > input_limit)
        {
            throw std::invalid_argument("Unsupported radial quadrature grid for batched tables");
        }
        const std::size_t radial_limit = std::min(table_limit / (2 * kmesh), output_limit / 2);
        const int radial_block = std::min<std::size_t>(rmesh, radial_limit);
        const std::size_t batch_limit = std::min(input_limit / kmesh, output_limit / (2 * radial_block));
        const int batch_size = std::min<std::size_t>(tasks.size(), std::min<std::size_t>(128, batch_limit));
        const std::size_t input_size = static_cast<std::size_t>(batch_size) * kmesh;
        const std::size_t table_size = static_cast<std::size_t>(2 * radial_block) * kmesh;
        const std::size_t output_size = static_cast<std::size_t>(batch_size) * 2 * radial_block;
        std::vector<double> input(input_size);
        std::vector<double> table(table_size);
        std::vector<double> output(output_size);
        double* input_ptr = input.data();
        double* table_ptr = table.data();
        double* output_ptr = output.data();
#ifdef __CUDA
        std::unique_ptr<container::Tensor> workspace;
        if (use_gpu)
        {
            const int64_t workspace_size = input_size + table_size + output_size;
            const container::TensorShape shape({workspace_size});
            workspace.reset(
                new container::Tensor(container::DataType::DT_DOUBLE, container::DeviceType::GpuDevice, shape));
            input_ptr = workspace->data<double>();
            table_ptr = input_ptr + input_size;
            output_ptr = table_ptr + table_size;
        }
#endif
        for (int radial_start = 0; radial_start < rmesh; radial_start += radial_block)
        {
            const int rows = std::min(radial_block, rmesh - radial_start);
            const int columns = 2 * rows;
            const double derivative_factor = -(l + 1.0) / (2 * l + 1.0);
            const double lower_factor = l / (l + 1.0);
            for (int ir = 0; ir < rows; ++ir)
            {
                const int radial = radial_start + ir;
                for (int ik = 0; ik < kmesh; ++ik)
                {
                    const double lower = l == 0 ? 0.0 : bessel[l - 1][radial][ik];
                    table[(2 * ir) * kmesh + ik] = bessel[l][radial][ik];
                    table[(2 * ir + 1) * kmesh + ik]
                        = derivative_factor * kpoint[ik] * (bessel[l + 1][radial][ik] - lower_factor * lower);
                }
            }
#ifdef __CUDA
            if (use_gpu)
            {
                const std::size_t count = static_cast<std::size_t>(columns) * kmesh;
                base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(
                    table_ptr,
                    table.data(),
                    count);
            }
#endif
            for (std::size_t start = 0; start < tasks.size(); start += batch_size)
            {
                const int count = std::min<std::size_t>(batch_size, tasks.size() - start);
                for (int it = 0; it < count; ++it)
                {
                    const Orb11& orbital = *tasks[start + it];
                    for (int ik = 0; ik < kmesh; ++ik)
                    {
                        const int weight = (ik == 0 || ik == kmesh - 1) ? 1 : (ik % 2 == 1 ? 4 : 2);
                        const double product = radial_product(orbital.nA, orbital.nB, ik);
                        input[it * kmesh + ik] = product * weight * dk * ModuleBase::FOUR_PI / 3.0;
                    }
                }
#ifdef __CUDA
                if (use_gpu)
                {
                    const std::size_t size = static_cast<std::size_t>(count) * kmesh;
                    base_device::memory::
                        synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_CPU>()(input_ptr,
                                                                                                          input.data(),
                                                                                                          size);
                }
#endif
#ifdef __CUDA
                if (use_gpu)
                {
                    const double one = 1.0;
                    const double zero = 0.0;
                    container::kernels::blas_gemm<double, container::DEVICE_GPU>()('T',
                                                                                   'N',
                                                                                   columns,
                                                                                   count,
                                                                                   kmesh,
                                                                                   &one,
                                                                                   table_ptr,
                                                                                   kmesh,
                                                                                   input_ptr,
                                                                                   kmesh,
                                                                                   &zero,
                                                                                   output_ptr,
                                                                                   columns);
                }
                else
#endif
                {
                    BlasConnector::gemm_cm('T',
                                           'N',
                                           columns,
                                           count,
                                           kmesh,
                                           1.0,
                                           table_ptr,
                                           kmesh,
                                           input_ptr,
                                           kmesh,
                                           0.0,
                                           output_ptr,
                                           columns,
                                           base_device::AbacusDevice_t::CpuDevice);
                }
#ifdef __CUDA
                if (use_gpu)
                {
                    const std::size_t size = static_cast<std::size_t>(count) * columns;
                    base_device::memory::
                        synchronize_memory_op<double, base_device::DEVICE_CPU, base_device::DEVICE_GPU>()(output.data(),
                                                                                                          output_ptr,
                                                                                                          size);
                }
#endif
                for (int it = 0; it < count; ++it)
                {
                    auto& values = tasks[start + it]->Table_r.at(l);
                    auto& derivatives = tasks[start + it]->Table_dr.at(l);
                    for (int ir = 0; ir < rows; ++ir)
                    {
                        values[radial_start + ir] = output[it * columns + 2 * ir];
                        derivatives[radial_start + ir] = output[it * columns + 2 * ir + 1];
                    }
                }
            }
        }
        // Preserve the analytic R=0 convention used by cal_overlap.
        if (l > 0)
        {
            std::vector<double> integrand(kmesh);
            for (Orb11* orbital: tasks)
            {
                for (int ik = 0; ik < kmesh; ++ik)
                {
                    integrand[ik] = radial_product(orbital->nA, orbital->nB, ik) * std::pow(kpoint[ik], l);
                }
                double integral = 0.0;
                ModuleBase::Integral::Simpson_Integral(kmesh, integrand.data(), dk, integral);
                orbital->Table_r.at(l)[0]
                    = ModuleBase::FOUR_PI / ModuleBase::Mathzone_Add1::dualfac(2 * l + 1) * integral;
            }
        }
    }
    ModuleBase::timer::end("Center2_Orb", "init_radial_table");
}

double Center2_Orb::Orb11::cal_overlap(const ModuleBase::Vector3<double>& RA,
                                       const ModuleBase::Vector3<double>& RB,
                                       const int& mA,
                                       const int& mB) const
{
    const double tiny1 = 1e-12; // why 1e-12?
    const double tiny2 = 1e-10; // why 1e-10?

    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > (RcutA + RcutB))
    {
        return 0.0;
    }

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();

    std::vector<double> rly;
    ModuleBase::Ylm::rl_sph_harm(LA + LB, // max LAB
                                 delta_R.x,
                                 delta_R.y,
                                 delta_R.z,
                                 rly);

    double overlap = 0.0;
    const int idx1 = this->MGT.get_lm_index(LA, mA);
    const int idx2 = this->MGT.get_lm_index(LB, mB);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, dr_);
    for (const auto& tb_r: this->Table_r)
    {
        const int LAB = tb_r.first;
        const double distance_power = std::pow(distance, LAB);
        const int phase = (LA - LB - LAB) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double Interp_Tlm
            = (distance > tiny2)
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(tb_r.second.data(), rmesh, dr_, distance)
                  : tb_r.second.at(0);

        for (int mAB = 0; mAB != 2 * LAB + 1; ++mAB)
        // const int mAB = mA + mB;
        {
            const int idx3 = this->MGT.get_lm_index(LAB, mAB);
            const double Gaunt_real_A_B_AB = *(Gaunt_Coefficients_ptr + idx3);
            if (0 == Gaunt_real_A_B_AB)
            {
                continue;
            }

            const double ylm_solid = rly[idx3];
            if (0 == ylm_solid)
            {
                continue;
            }
            const double ylm_real = (distance > tiny2) ? ylm_solid / distance_power : ylm_solid;

            overlap +=
                //			pow(2*PI,1.5)
                i_exp * Gaunt_real_A_B_AB * Interp_Tlm * ylm_real;
        }
    }

    return overlap;
}

ModuleBase::Vector3<double> Center2_Orb::Orb11::cal_grad_overlap( // caoyu add 2021-11-19
    const ModuleBase::Vector3<double>& RA,
    const ModuleBase::Vector3<double>& RB,
    const int& mA,
    const int& mB) const
{
    const double tiny1 = 1e-12; // same as `cal_overlap`
    const double tiny2 = 1e-10; // same as `cal_overlap`

    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = (distance_true >= tiny1) ? distance_true : distance_true + tiny1;
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > (RcutA + RcutB))
    {
        return ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }

    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const int idx1 = this->MGT.get_lm_index(LA, mA);
    const int idx2 = this->MGT.get_lm_index(LB, mB);
    const double* Gaunt_Coefficients_ptr = &(this->MGT.Gaunt_Coefficients(idx1, idx2, 0));

    const int LAB2 = (LA + LB + 1) * (LA + LB + 1);
    std::vector<double> rly(LAB2);
    std::vector<ModuleBase::Vector3<double>> grly;
    std::vector<double> tmp_grly(LAB2 * 3);
    ModuleBase::Ylm::grad_rl_sph_harm(LA + LB, delta_R.x, delta_R.y, delta_R.z, rly.data(), tmp_grly.data());
    grly.reserve(LAB2);
    for (int i=0; i<LAB2; ++i)
    {
        grly.emplace_back(tmp_grly[i*3], tmp_grly[i*3+1], tmp_grly[i*3+2]);
    }

    ModuleBase::Vector3<double> grad_overlap(0.0, 0.0, 0.0);

    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, dr_);
    for (const auto& tb_r: this->Table_r)
    {
        const int LAB = tb_r.first;
        const double distance_power = std::pow(distance, LAB);
        const int phase = (LA - LB - LAB) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double Interp_Tlm
            = (distance > tiny2)
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(tb_r.second.data(), rmesh, dr_, distance)
                  : tb_r.second.at(0);
        const double grad_Interp_Tlm
            = (distance > tiny2)
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(this->Table_dr.at(LAB).data(), rmesh, dr_, distance)
                        - Interp_Tlm * LAB / distance
                  : 0.0;
        for (int mAB = 0; mAB != 2 * LAB + 1; ++mAB)
        // const int mAB = mA + mB;
        {
            const int idx3 = this->MGT.get_lm_index(LAB, mAB);
            const double Gaunt_real_A_B_AB = *(Gaunt_Coefficients_ptr + idx3);
            if (0 == Gaunt_real_A_B_AB)
            {
                continue;
            }

            const double ylm_solid = rly[idx3];
            const double ylm_real = (distance > tiny2) ? ylm_solid / distance_power : ylm_solid;

            const ModuleBase::Vector3<double> gylm_solid = grly[idx3];
            const ModuleBase::Vector3<double> gylm_real = (distance > tiny2) ? gylm_solid / distance_power : gylm_solid;

            grad_overlap += i_exp //			pow(2*PI,1.5)
                            * Gaunt_real_A_B_AB
                            * (Interp_Tlm * gylm_real + grad_Interp_Tlm * ylm_real * delta_R / distance);
        }
    }

    return grad_overlap;
}

template <typename T>
std::vector<T> Center2_Orb::Orb11::contract_angular(const std::vector<T>& components) const
{
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const int count_A = 2 * LA + 1;
    const int count_B = 2 * LB + 1;
    const int count = count_A * count_B;
    std::vector<T> result(count);
    for (int ma = 0; ma < count_A; ++ma)
    {
        const int idx_A = this->MGT.get_lm_index(LA, ma);
        for (int mb = 0; mb < count_B; ++mb)
        {
            const int idx_B = this->MGT.get_lm_index(LB, mb);
            const double* coefficients = &this->MGT.Gaunt_Coefficients(idx_A, idx_B, 0);
            T& value = result[ma * count_B + mb];
            // Keep the scalar path's angular-momentum accumulation order.
            for (const std::pair<const int, std::vector<double>>& table: this->Table_r)
            {
                const int l = table.first;
                for (int m = 0; m <= 2 * l; ++m)
                {
                    const int idx = this->MGT.get_lm_index(l, m);
                    if (coefficients[idx] != 0.0)
                    {
                        value += coefficients[idx] * components[idx];
                    }
                }
            }
        }
    }
    return result;
}

std::vector<double> Center2_Orb::Orb11::cal_overlap(const ModuleBase::Vector3<double>& RA,
                                                 const ModuleBase::Vector3<double>& RB) const
{
    const double tiny1 = 1e-12;
    const double tiny2 = 1e-10;
    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = distance_true >= tiny1 ? distance_true : distance_true + tiny1;
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > RcutA + RcutB)
    {
        const int count = (2 * LA + 1) * (2 * LB + 1);
        return std::vector<double>(count);
    }

    const int lmax = LA + LB;
    std::vector<double> rly;
    ModuleBase::Ylm::rl_sph_harm(lmax, delta_R.x, delta_R.y, delta_R.z, rly);
    const int angular_count = (lmax + 1) * (lmax + 1);
    std::vector<double> components(angular_count);
    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, this->dr_);
    for (const std::pair<const int, std::vector<double>>& table: this->Table_r)
    {
        const int l = table.first;
        const double distance_power = std::pow(distance, l);
        const int phase = (LA - LB - l) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double radial
            = distance > tiny2
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(table.second.data(), rmesh, dr_, distance)
                  : table.second.at(0);
        for (int m = 0; m <= 2 * l; ++m)
        {
            const int idx = this->MGT.get_lm_index(l, m);
            const double angular = distance > tiny2 ? rly[idx] / distance_power : rly[idx];
            components[idx] = i_exp * radial * angular;
        }
    }
    return this->contract_angular(components);
}

std::vector<ModuleBase::Vector3<double>> Center2_Orb::Orb11::cal_grad_overlap(
    const ModuleBase::Vector3<double>& RA,
    const ModuleBase::Vector3<double>& RB) const
{
    const double tiny1 = 1e-12;
    const double tiny2 = 1e-10;
    const ModuleBase::Vector3<double> delta_R = RB - RA;
    const double distance_true = delta_R.norm();
    const double distance = distance_true >= tiny1 ? distance_true : distance_true + tiny1;
    const int LA = this->nA.getL();
    const int LB = this->nB.getL();
    const double RcutA = this->nA.getRcut();
    const double RcutB = this->nB.getRcut();
    if (distance > RcutA + RcutB)
    {
        const int count = (2 * LA + 1) * (2 * LB + 1);
        return std::vector<ModuleBase::Vector3<double>>(count);
    }

    const int lmax = LA + LB;
    const int angular_count = (lmax + 1) * (lmax + 1);
    const int gradient_count = 3 * angular_count;
    std::vector<double> rly(angular_count);
    std::vector<double> grly(gradient_count);
    ModuleBase::Ylm::grad_rl_sph_harm(lmax, delta_R.x, delta_R.y, delta_R.z, rly.data(), grly.data());
    std::vector<ModuleBase::Vector3<double>> components(angular_count);
    const int rmesh = Center2_Orb::get_rmesh(RcutA, RcutB, this->dr_);
    for (const std::pair<const int, std::vector<double>>& table: this->Table_r)
    {
        const int l = table.first;
        const double distance_power = std::pow(distance, l);
        const int phase = (LA - LB - l) / 2;
        const double i_exp = std::pow(-1.0, phase);
        const double radial
            = distance > tiny2
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(table.second.data(), rmesh, dr_, distance)
                  : table.second.at(0);
        const double derivative
            = distance > tiny2
                  ? ModuleBase::PolyInt::Polynomial_Interpolation(this->Table_dr.at(l).data(), rmesh, dr_, distance)
                        - radial * l / distance
                  : 0.0;
        for (int m = 0; m <= 2 * l; ++m)
        {
            const int idx = this->MGT.get_lm_index(l, m);
            const double angular = distance > tiny2 ? rly[idx] / distance_power : rly[idx];
            const int offset = 3 * idx;
            const double gx = grly[offset];
            const double gy = grly[offset + 1];
            const double gz = grly[offset + 2];
            const ModuleBase::Vector3<double> solid_gradient(gx, gy, gz);
            const ModuleBase::Vector3<double> angular_gradient
                = distance > tiny2 ? solid_gradient / distance_power : solid_gradient;
            components[idx] = i_exp * (radial * angular_gradient + derivative * angular * delta_R / distance);
        }
    }
    return this->contract_angular(components);
}

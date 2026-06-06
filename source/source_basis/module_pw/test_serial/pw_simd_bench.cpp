#include <chrono>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "source_base/matrix3.h"
#include "source_base/timer.h"

#include "../pw_basis.h"
#include "../pw_basis_k.h"

namespace
{

using Clock = std::chrono::steady_clock;

template <typename Func>
double measure_seconds(Func&& func)
{
    const auto start = Clock::now();
    func();
    const auto end = Clock::now();
    return std::chrono::duration<double>(end - start).count();
}

void print_metric(const std::string& name, const double value)
{
    std::cout << "METRIC " << name << " " << std::fixed << std::setprecision(9) << value << '\n';
}

void print_timer_metric(const std::string& class_name, const std::string& timer_name)
{
    const auto class_it = ModuleBase::timer::timer_pool.find(class_name);
    if (class_it == ModuleBase::timer::timer_pool.end())
    {
        return;
    }
    const auto timer_it = class_it->second.find(timer_name);
    if (timer_it == class_it->second.end())
    {
        return;
    }
    print_metric("timer." + class_name + "." + timer_name + ".seconds", timer_it->second.cpu_second);
    print_metric("timer." + class_name + "." + timer_name + ".calls", static_cast<double>(timer_it->second.calls));
}

double run_pw_basis_roundtrip(ModulePW::PW_Basis& basis,
                              const int npw,
                              const int nrxx,
                              const int repeat)
{
    std::vector<std::complex<double>> recip_in(npw);
    std::vector<std::complex<double>> real_space(nrxx);
    std::vector<std::complex<double>> recip_out(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        recip_in[ig] = std::complex<double>((ig % 17 - 8) / 13.0, (ig % 19 - 9) / 17.0);
    }

    for (int warmup = 0; warmup < 3; ++warmup)
    {
        basis.recip2real(recip_in.data(), real_space.data());
        basis.real2recip(real_space.data(), recip_out.data());
    }

    return measure_seconds([&]() {
        for (int i = 0; i < repeat; ++i)
        {
            basis.recip2real(recip_in.data(), real_space.data());
            basis.real2recip(real_space.data(), recip_out.data());
        }
    });
}

double run_pw_basis_k_roundtrip(ModulePW::PW_Basis_K& basis,
                                const int npw,
                                const int nrxx,
                                const int ik,
                                const int repeat)
{
    std::vector<std::complex<double>> recip_in(npw);
    std::vector<std::complex<double>> real_space(nrxx);
    std::vector<std::complex<double>> recip_out(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        recip_in[ig] = std::complex<double>((ig % 17 - 8) / 13.0, (ig % 19 - 9) / 17.0);
    }

    for (int warmup = 0; warmup < 3; ++warmup)
    {
        basis.recip2real(recip_in.data(), real_space.data(), ik);
        basis.real2recip(real_space.data(), recip_out.data(), ik);
    }

    return measure_seconds([&]() {
        for (int i = 0; i < repeat; ++i)
        {
            basis.recip2real(recip_in.data(), real_space.data(), ik);
            basis.real2recip(real_space.data(), recip_out.data(), ik);
        }
    });
}

void bench_pw_basis_medium()
{
    ModuleBase::timer::timer_pool.clear();
    ModulePW::PW_Basis basis;
    basis.initgrids(2.0, ModuleBase::Matrix3(1.0, 0.0, 1.0,
                                             0.0, 2.0, 0.0,
                                             0.0, 0.0, 2.0),
                    30.0);
    basis.initparameters(false, 20.0, 2, false);
    basis.setuptransform();

    const int repeat = 4096;
    const double elapsed = run_pw_basis_roundtrip(basis, basis.npw, basis.nrxx, repeat);
    print_metric("PW_Basis.medium.roundtrip.wall", elapsed);
    print_metric("PW_Basis.medium.roundtrip.ms_per_op", elapsed / repeat * 1000.0);
    print_metric("PW_Basis.medium.nrxx", static_cast<double>(basis.nrxx));
    print_metric("PW_Basis.medium.npw", static_cast<double>(basis.npw));

    print_timer_metric("PW_Basis", "real2recip");
    print_timer_metric("PW_Basis", "recip2real");
    print_timer_metric("PW_Basis", "gatherp_copy_serial");
    print_timer_metric("PW_Basis", "gathers_copy_serial");
}

void bench_pw_basis_large()
{
    ModuleBase::timer::timer_pool.clear();
    ModulePW::PW_Basis basis;
    basis.initgrids(2.0, ModuleBase::Matrix3(2.0, 0.0, 0.0,
                                             0.0, 2.0, 0.0,
                                             0.0, 0.0, 2.0),
                    40.0);
    basis.initparameters(false, 25.0, 2, false);
    basis.setuptransform();

    const int repeat = 2048;
    const double elapsed = run_pw_basis_roundtrip(basis, basis.npw, basis.nrxx, repeat);
    print_metric("PW_Basis.large.roundtrip.wall", elapsed);
    print_metric("PW_Basis.large.roundtrip.ms_per_op", elapsed / repeat * 1000.0);
    print_metric("PW_Basis.large.nrxx", static_cast<double>(basis.nrxx));
    print_metric("PW_Basis.large.npw", static_cast<double>(basis.npw));

    print_timer_metric("PW_Basis", "real2recip");
    print_timer_metric("PW_Basis", "recip2real");
    print_timer_metric("PW_Basis", "gatherp_copy_serial");
    print_timer_metric("PW_Basis", "gathers_copy_serial");
}

void bench_pw_basis_k_medium()
{
    ModuleBase::timer::timer_pool.clear();
    ModulePW::PW_Basis_K basis("cpu", "double");
    const ModuleBase::Vector3<double> kvec_d[1] = {{0.0, 0.0, 0.0}};
    basis.initgrids(2.0, ModuleBase::Matrix3(1.0, 0.0, 1.0,
                                             0.0, 2.0, 0.0,
                                             0.0, 0.0, 2.0),
                    30.0);
    basis.initparameters(false, 20.0, 1, kvec_d, 2, false);
    basis.setuptransform();

    const int repeat = 4096;
    const int ik = 0;
    const double elapsed = run_pw_basis_k_roundtrip(basis, basis.npwk[ik], basis.nrxx, ik, repeat);
    print_metric("PW_Basis_K.medium.roundtrip.wall", elapsed);
    print_metric("PW_Basis_K.medium.roundtrip.ms_per_op", elapsed / repeat * 1000.0);
    print_metric("PW_Basis_K.medium.nrxx", static_cast<double>(basis.nrxx));
    print_metric("PW_Basis_K.medium.npw", static_cast<double>(basis.npwk[ik]));

    print_timer_metric("PW_Basis_K", "real2recip");
    print_timer_metric("PW_Basis_K", "recip2real");
    print_timer_metric("PW_Basis_K", "gatherp_copy_serial");
    print_timer_metric("PW_Basis_K", "gathers_copy_serial");
}

} // namespace

int main()
{
    ModuleBase::timer::enable();

    bench_pw_basis_medium();
    bench_pw_basis_large();
    bench_pw_basis_k_medium();

    std::ofstream nullstream("/dev/null");
    ModuleBase::timer::finish(nullstream, false, false);
    return 0;
}

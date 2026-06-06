#include "source_base/matrix3.h"
#include "source_base/timer.h"

#include "../pw_basis.h"
#include "../pw_basis_k.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef __MPI
#include "mpi.h"
#endif

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

void bench_pw_basis()
{
    constexpr int repeat_calls = 2000;
    ModuleBase::timer::timer_pool.clear();

    ModulePW::PW_Basis basis;
    const ModuleBase::Matrix3 latvec(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const double lat0 = 10.0;
    const double wfcecut = 50.0;
    const double rhoecut = 4.0 * wfcecut;
    const int distribution_type = 1;

    basis.initgrids(lat0, latvec, rhoecut);
    basis.initparameters(false, wfcecut, distribution_type, true);

    print_metric("PW_Basis.setuptransform.wall", measure_seconds([&]() { basis.setuptransform(); }));
    print_metric("PW_Basis.collect_local_pw.first.wall", measure_seconds([&]() { basis.collect_local_pw(); }));
    print_metric("PW_Basis.collect_local_pw.repeat.wall",
                 measure_seconds([&]() {
                     for (int i = 0; i < repeat_calls; ++i)
                     {
                         basis.collect_local_pw();
                     }
                 }));
    print_metric("PW_Basis.collect_uniqgg.first.wall", measure_seconds([&]() { basis.collect_uniqgg(); }));
    print_metric("PW_Basis.collect_uniqgg.repeat.wall",
                 measure_seconds([&]() {
                     for (int i = 0; i < repeat_calls; ++i)
                     {
                         basis.collect_uniqgg();
                     }
                 }));

    print_timer_metric("PW_Basis", "setuptransform");
    print_timer_metric("PW_Basis", "collect_local_pw");
    print_timer_metric("PW_Basis", "collect_local_pw_cache_hit");
    print_timer_metric("PW_Basis", "collect_local_pw_cache_build");
    print_timer_metric("PW_Basis", "collect_uniqgg");
    print_timer_metric("PW_Basis", "collect_uniqgg_cache_hit");
    print_timer_metric("PW_Basis", "collect_uniqgg_cache_build");
}

void bench_pw_basis_k()
{
    constexpr int repeat_calls = 2000;
    ModuleBase::timer::timer_pool.clear();

    ModulePW::PW_Basis_K basis("cpu", "double");
    const ModuleBase::Matrix3 latvec(10.0, 0.0, 0.0,
                                     0.0, 10.0, 0.0,
                                     0.0, 0.0, 10.0);
    const double lat0 = 1.8897261254578281;
    const double gridecut = 10.0;
    const bool gamma_only = true;
    const double gk_ecut = 11.0;
    const int nks = 3;
    const ModuleBase::Vector3<double> kvec_d[3] = {{0.0, 0.0, 0.0}, {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}};
    const int distribution_type = 1;
    const bool xprime = true;

    basis.initgrids(lat0, latvec, gridecut);
    basis.initparameters(gamma_only, gk_ecut, nks, kvec_d, distribution_type, xprime);

    print_metric("PW_Basis_K.setuptransform.wall", measure_seconds([&]() { basis.setuptransform(); }));
    print_metric("PW_Basis_K.collect_local_pw.first.wall", measure_seconds([&]() { basis.collect_local_pw(); }));
    print_metric("PW_Basis_K.collect_local_pw.repeat.wall",
                 measure_seconds([&]() {
                     for (int i = 0; i < repeat_calls; ++i)
                     {
                         basis.collect_local_pw();
                     }
                 }));
    print_metric("PW_Basis_K.collect_local_pw.gk2_rebuild.wall",
                 measure_seconds([&]() {
                     for (int i = 0; i < repeat_calls; ++i)
                     {
                         basis.collect_local_pw(1.0, 0.5, 0.2);
                     }
                 }));

    print_timer_metric("PW_Basis_K", "setuptransform");
    print_timer_metric("PW_Basis_K", "collect_local_pw");
    print_timer_metric("PW_Basis_K", "collect_local_pw_cache_hit");
    print_timer_metric("PW_Basis_K", "collect_local_pw_build_gcar");
    print_timer_metric("PW_Basis_K", "collect_local_pw_build_gk2");
}

} // namespace

int main()
{
#ifdef __MPI
    int argc = 0;
    char** argv = nullptr;
    MPI_Init(&argc, &argv);
#endif
    bench_pw_basis();
    bench_pw_basis_k();
#ifdef __MPI
    MPI_Finalize();
#endif
    return 0;
}

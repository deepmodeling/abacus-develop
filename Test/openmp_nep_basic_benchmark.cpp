#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{

struct Vec3
{
    double x;
    double y;
    double z;
};

struct Matrix3
{
    double m[3][3];
};

struct BenchmarkData
{
    int natom = 2000000;
    int ntype = 4;
    double md_dt = 0.5;
    double lat0 = 10.0;
    double fact_f = 0.0433641153087705;
    Matrix3 gt{};

    std::vector<Vec3> pos;
    std::vector<Vec3> vel;
    std::vector<Vec3> base_vel;
    std::vector<Vec3> force;
    std::vector<double> mass;
    std::vector<std::array<int, 3>> ionmbl;
    std::vector<double> force_matrix;

    std::vector<int> atom_type_index;
    std::vector<int> atom_local_index;
    std::vector<int> type_offset;
    std::vector<double> tau_x;
    std::vector<double> tau_y;
    std::vector<double> tau_z;
    std::vector<double> nep_coord;
    std::vector<double> nep_e;
    std::vector<double> nep_f;
    std::vector<double> nep_v;
    std::vector<double> nep_force;
};

struct Measurement
{
    std::string kernel;
    int threads = 1;
    double serial_ms = 0.0;
    double omp_ms = 0.0;
    double max_abs_diff = 0.0;
    double checksum = 0.0;
};

using Clock = std::chrono::steady_clock;
volatile double global_sink = 0.0;

double elapsed_ms(const Clock::time_point& start, const Clock::time_point& end)
{
    return std::chrono::duration<double, std::milli>(end - start).count();
}

int parse_int_arg(int argc, char** argv, const std::string& name, int fallback)
{
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (argv[i] == name)
        {
            return std::atoi(argv[i + 1]);
        }
    }
    return fallback;
}

double value_for_index(const int i, const double scale)
{
    return scale * (1.0 + static_cast<double>((i * 17) % 97) / 97.0);
}

BenchmarkData make_data(const int natom)
{
    BenchmarkData data;
    data.natom = natom;
    data.gt = {{{1.02, 0.01, 0.03}, {0.02, 0.98, 0.04}, {0.01, 0.02, 1.01}}};

    data.pos.resize(natom);
    data.vel.resize(natom);
    data.base_vel.resize(natom);
    data.force.resize(natom);
    data.mass.resize(natom);
    data.ionmbl.resize(natom);
    data.force_matrix.resize(static_cast<std::size_t>(natom) * 3);
    data.atom_type_index.resize(natom);
    data.atom_local_index.resize(natom);
    data.nep_coord.resize(static_cast<std::size_t>(natom) * 3);
    data.nep_e.resize(natom);
    data.nep_f.resize(static_cast<std::size_t>(natom) * 3);
    data.nep_v.resize(static_cast<std::size_t>(natom) * 9);
    data.nep_force.resize(static_cast<std::size_t>(natom) * 3);

    data.type_offset.assign(data.ntype + 1, 0);
    for (int it = 0; it < data.ntype; ++it)
    {
        data.type_offset[it + 1] = data.type_offset[it] + natom / data.ntype + (it < natom % data.ntype ? 1 : 0);
    }
    data.tau_x.resize(natom);
    data.tau_y.resize(natom);
    data.tau_z.resize(natom);

    for (int i = 0; i < natom; ++i)
    {
        data.base_vel[i] = {value_for_index(i, 0.001), value_for_index(i + 11, 0.0012), value_for_index(i + 23, 0.0009)};
        data.vel[i] = data.base_vel[i];
        data.force[i] = {value_for_index(i + 3, 0.004), value_for_index(i + 5, 0.003), value_for_index(i + 7, 0.005)};
        data.mass[i] = 10.0 + static_cast<double>(i % 13);
        data.ionmbl[i] = {1, (i % 11) == 0 ? 0 : 1, (i % 17) == 0 ? 0 : 1};
        data.force_matrix[3 * static_cast<std::size_t>(i)] = data.force[i].x;
        data.force_matrix[3 * static_cast<std::size_t>(i) + 1] = data.force[i].y;
        data.force_matrix[3 * static_cast<std::size_t>(i) + 2] = data.force[i].z;
        data.nep_e[i] = value_for_index(i + 13, 0.02);
        data.nep_f[i] = value_for_index(i + 29, 0.03);
        data.nep_f[static_cast<std::size_t>(i) + natom] = value_for_index(i + 31, 0.02);
        data.nep_f[static_cast<std::size_t>(i) + 2 * static_cast<std::size_t>(natom)] = value_for_index(i + 37, 0.04);
    }

    for (int it = 0; it < data.ntype; ++it)
    {
        for (int local = 0; local < data.type_offset[it + 1] - data.type_offset[it]; ++local)
        {
            const int global = data.type_offset[it] + local;
            data.atom_type_index[global] = it;
            data.atom_local_index[global] = local;
            data.tau_x[global] = value_for_index(global + 41, 0.2);
            data.tau_y[global] = value_for_index(global + 43, 0.3);
            data.tau_z[global] = value_for_index(global + 47, 0.4);
        }
    }

    for (int j = 0; j < 9; ++j)
    {
        for (int i = 0; i < natom; ++i)
        {
            data.nep_v[static_cast<std::size_t>(j) * natom + i] = value_for_index(i + j * 19, 0.0003 * (j + 1));
        }
    }

    return data;
}

template <typename F>
double time_repeated(const int repeat, F&& f)
{
    const auto start = Clock::now();
    for (int r = 0; r < repeat; ++r)
    {
        f();
    }
    return elapsed_ms(start, Clock::now()) / repeat;
}

void update_pos_serial(const BenchmarkData& data, std::vector<Vec3>& pos)
{
    for (int i = 0; i < data.natom; ++i)
    {
        double px = data.ionmbl[i][0] ? data.vel[i].x * data.md_dt / data.lat0 : 0.0;
        double py = data.ionmbl[i][1] ? data.vel[i].y * data.md_dt / data.lat0 : 0.0;
        double pz = data.ionmbl[i][2] ? data.vel[i].z * data.md_dt / data.lat0 : 0.0;
        pos[i].x = px * data.gt.m[0][0] + py * data.gt.m[1][0] + pz * data.gt.m[2][0];
        pos[i].y = px * data.gt.m[0][1] + py * data.gt.m[1][1] + pz * data.gt.m[2][1];
        pos[i].z = px * data.gt.m[0][2] + py * data.gt.m[1][2] + pz * data.gt.m[2][2];
    }
}

void update_pos_omp(const BenchmarkData& data, std::vector<Vec3>& pos)
{
#pragma omp parallel for schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        double px = data.ionmbl[i][0] ? data.vel[i].x * data.md_dt / data.lat0 : 0.0;
        double py = data.ionmbl[i][1] ? data.vel[i].y * data.md_dt / data.lat0 : 0.0;
        double pz = data.ionmbl[i][2] ? data.vel[i].z * data.md_dt / data.lat0 : 0.0;
        pos[i].x = px * data.gt.m[0][0] + py * data.gt.m[1][0] + pz * data.gt.m[2][0];
        pos[i].y = px * data.gt.m[0][1] + py * data.gt.m[1][1] + pz * data.gt.m[2][1];
        pos[i].z = px * data.gt.m[0][2] + py * data.gt.m[1][2] + pz * data.gt.m[2][2];
    }
}

void update_vel_serial(const BenchmarkData& data, std::vector<Vec3>& vel)
{
    for (int i = 0; i < data.natom; ++i)
    {
        if (data.ionmbl[i][0])
        {
            vel[i].x += 0.5 * data.force[i].x * data.md_dt / data.mass[i];
        }
        if (data.ionmbl[i][1])
        {
            vel[i].y += 0.5 * data.force[i].y * data.md_dt / data.mass[i];
        }
        if (data.ionmbl[i][2])
        {
            vel[i].z += 0.5 * data.force[i].z * data.md_dt / data.mass[i];
        }
    }
}

void update_vel_omp(const BenchmarkData& data, std::vector<Vec3>& vel)
{
#pragma omp parallel for schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        if (data.ionmbl[i][0])
        {
            vel[i].x += 0.5 * data.force[i].x * data.md_dt / data.mass[i];
        }
        if (data.ionmbl[i][1])
        {
            vel[i].y += 0.5 * data.force[i].y * data.md_dt / data.mass[i];
        }
        if (data.ionmbl[i][2])
        {
            vel[i].z += 0.5 * data.force[i].z * data.md_dt / data.mass[i];
        }
    }
}

double kinetic_serial(const BenchmarkData& data)
{
    double ke = 0.0;
    for (int i = 0; i < data.natom; ++i)
    {
        ke += 0.5 * data.mass[i]
              * (data.vel[i].x * data.vel[i].x + data.vel[i].y * data.vel[i].y + data.vel[i].z * data.vel[i].z);
    }
    return ke;
}

double kinetic_omp(const BenchmarkData& data)
{
    double ke = 0.0;
#pragma omp parallel for reduction(+:ke) schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        ke += 0.5 * data.mass[i]
              * (data.vel[i].x * data.vel[i].x + data.vel[i].y * data.vel[i].y + data.vel[i].z * data.vel[i].z);
    }
    return ke;
}

std::array<double, 9> temp_vector_serial(const BenchmarkData& data)
{
    std::array<double, 9> t{};
    for (int i = 0; i < data.natom; ++i)
    {
        const double m = data.mass[i];
        const double vx = data.vel[i].x;
        const double vy = data.vel[i].y;
        const double vz = data.vel[i].z;
        t[0] += m * vx * vx;
        t[1] += m * vx * vy;
        t[2] += m * vx * vz;
        t[3] += m * vy * vx;
        t[4] += m * vy * vy;
        t[5] += m * vy * vz;
        t[6] += m * vz * vx;
        t[7] += m * vz * vy;
        t[8] += m * vz * vz;
    }
    return t;
}

std::array<double, 9> temp_vector_omp(const BenchmarkData& data)
{
    double t00 = 0.0;
    double t01 = 0.0;
    double t02 = 0.0;
    double t10 = 0.0;
    double t11 = 0.0;
    double t12 = 0.0;
    double t20 = 0.0;
    double t21 = 0.0;
    double t22 = 0.0;

#pragma omp parallel for reduction(+:t00, t01, t02, t10, t11, t12, t20, t21, t22) schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        const double m = data.mass[i];
        const double vx = data.vel[i].x;
        const double vy = data.vel[i].y;
        const double vz = data.vel[i].z;
        t00 += m * vx * vx;
        t01 += m * vx * vy;
        t02 += m * vx * vz;
        t10 += m * vy * vx;
        t11 += m * vy * vy;
        t12 += m * vy * vz;
        t20 += m * vz * vx;
        t21 += m * vz * vy;
        t22 += m * vz * vz;
    }
    return {t00, t01, t02, t10, t11, t12, t20, t21, t22};
}

void force_copy_serial(const BenchmarkData& data, std::vector<Vec3>& out)
{
    for (int i = 0; i < data.natom; ++i)
    {
        out[i].x = data.force_matrix[3 * static_cast<std::size_t>(i)];
        out[i].y = data.force_matrix[3 * static_cast<std::size_t>(i) + 1];
        out[i].z = data.force_matrix[3 * static_cast<std::size_t>(i) + 2];
    }
}

void force_copy_omp(const BenchmarkData& data, std::vector<Vec3>& out)
{
#pragma omp parallel for schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        out[i].x = data.force_matrix[3 * static_cast<std::size_t>(i)];
        out[i].y = data.force_matrix[3 * static_cast<std::size_t>(i) + 1];
        out[i].z = data.force_matrix[3 * static_cast<std::size_t>(i) + 2];
    }
}

void nep_coord_serial(const BenchmarkData& data, std::vector<double>& coord)
{
    const int nat = data.natom;
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = data.atom_type_index[iat];
        const int ia = data.atom_local_index[iat];
        const int index = data.type_offset[it] + ia;
        coord[iat] = data.tau_x[index] * data.lat0;
        coord[static_cast<std::size_t>(iat) + nat] = data.tau_y[index] * data.lat0;
        coord[static_cast<std::size_t>(iat) + 2 * static_cast<std::size_t>(nat)] = data.tau_z[index] * data.lat0;
    }
}

void nep_coord_omp(const BenchmarkData& data, std::vector<double>& coord)
{
    const int nat = data.natom;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int iat = 0; iat < nat; ++iat)
    {
        const int it = data.atom_type_index[iat];
        const int ia = data.atom_local_index[iat];
        const int index = data.type_offset[it] + ia;
        coord[iat] = data.tau_x[index] * data.lat0;
        coord[static_cast<std::size_t>(iat) + nat] = data.tau_y[index] * data.lat0;
        coord[static_cast<std::size_t>(iat) + 2 * static_cast<std::size_t>(nat)] = data.tau_z[index] * data.lat0;
    }
}

double nep_energy_serial(const BenchmarkData& data)
{
    return std::accumulate(data.nep_e.begin(), data.nep_e.end(), 0.0);
}

double nep_energy_omp(const BenchmarkData& data)
{
    double energy = 0.0;
#pragma omp parallel for reduction(+:energy) schedule(static) if (data.natom >= 256)
    for (int i = 0; i < data.natom; ++i)
    {
        energy += data.nep_e[i];
    }
    return energy;
}

void nep_force_serial(const BenchmarkData& data, std::vector<double>& out)
{
    const int nat = data.natom;
    for (int i = 0; i < nat; ++i)
    {
        out[3 * static_cast<std::size_t>(i)] = data.nep_f[i] * data.fact_f;
        out[3 * static_cast<std::size_t>(i) + 1] = data.nep_f[static_cast<std::size_t>(i) + nat] * data.fact_f;
        out[3 * static_cast<std::size_t>(i) + 2] = data.nep_f[static_cast<std::size_t>(i) + 2 * static_cast<std::size_t>(nat)] * data.fact_f;
    }
}

void nep_force_omp(const BenchmarkData& data, std::vector<double>& out)
{
    const int nat = data.natom;
#pragma omp parallel for schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        out[3 * static_cast<std::size_t>(i)] = data.nep_f[i] * data.fact_f;
        out[3 * static_cast<std::size_t>(i) + 1] = data.nep_f[static_cast<std::size_t>(i) + nat] * data.fact_f;
        out[3 * static_cast<std::size_t>(i) + 2] = data.nep_f[static_cast<std::size_t>(i) + 2 * static_cast<std::size_t>(nat)] * data.fact_f;
    }
}

std::array<double, 9> nep_virial_serial(const BenchmarkData& data)
{
    std::array<double, 9> sum{};
    const int nat = data.natom;
    for (int j = 0; j < 9; ++j)
    {
        for (int i = 0; i < nat; ++i)
        {
            sum[j] += data.nep_v[static_cast<std::size_t>(j) * nat + i];
        }
    }
    return sum;
}

std::array<double, 9> nep_virial_omp(const BenchmarkData& data)
{
    const int nat = data.natom;
    double v0 = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double v3 = 0.0;
    double v4 = 0.0;
    double v5 = 0.0;
    double v6 = 0.0;
    double v7 = 0.0;
    double v8 = 0.0;
#pragma omp parallel for reduction(+:v0, v1, v2, v3, v4, v5, v6, v7, v8) schedule(static) if (nat >= 256)
    for (int i = 0; i < nat; ++i)
    {
        v0 += data.nep_v[i];
        v1 += data.nep_v[static_cast<std::size_t>(nat) + i];
        v2 += data.nep_v[2 * static_cast<std::size_t>(nat) + i];
        v3 += data.nep_v[3 * static_cast<std::size_t>(nat) + i];
        v4 += data.nep_v[4 * static_cast<std::size_t>(nat) + i];
        v5 += data.nep_v[5 * static_cast<std::size_t>(nat) + i];
        v6 += data.nep_v[6 * static_cast<std::size_t>(nat) + i];
        v7 += data.nep_v[7 * static_cast<std::size_t>(nat) + i];
        v8 += data.nep_v[8 * static_cast<std::size_t>(nat) + i];
    }
    return {v0, v1, v2, v3, v4, v5, v6, v7, v8};
}

double max_abs_diff(const std::vector<Vec3>& a, const std::vector<Vec3>& b)
{
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        diff = std::max(diff, std::abs(a[i].x - b[i].x));
        diff = std::max(diff, std::abs(a[i].y - b[i].y));
        diff = std::max(diff, std::abs(a[i].z - b[i].z));
    }
    return diff;
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b)
{
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        diff = std::max(diff, std::abs(a[i] - b[i]));
    }
    return diff;
}

double max_abs_diff(const std::array<double, 9>& a, const std::array<double, 9>& b)
{
    double diff = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        diff = std::max(diff, std::abs(a[i] - b[i]));
    }
    return diff;
}

double checksum(const std::vector<Vec3>& a)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); i += 4096)
    {
        sum += a[i].x + 0.5 * a[i].y + 0.25 * a[i].z;
    }
    return sum;
}

double checksum(const std::vector<double>& a)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); i += 4096)
    {
        sum += a[i];
    }
    return sum;
}

double checksum(const std::array<double, 9>& a)
{
    return std::accumulate(a.begin(), a.end(), 0.0);
}

void print_measurement(const Measurement& m)
{
    const double speedup = m.omp_ms > 0.0 ? m.serial_ms / m.omp_ms : 0.0;
    const double efficiency = m.threads > 0 ? speedup / m.threads : 0.0;
    std::cout << m.kernel << ',' << m.threads << ','
              << std::fixed << std::setprecision(6)
              << m.serial_ms << ',' << m.omp_ms << ',' << speedup << ',' << efficiency << ','
              << std::scientific << m.max_abs_diff << ',' << m.checksum << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    const int threads = parse_int_arg(argc, argv, "--threads", 1);
    const int natom = parse_int_arg(argc, argv, "--natoms", 2000000);
    const int repeat = parse_int_arg(argc, argv, "--repeat", 5);

#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif

    BenchmarkData data = make_data(natom);
    std::vector<Vec3> serial_vec(natom);
    std::vector<Vec3> omp_vec(natom);
    std::vector<double> serial_doubles(static_cast<std::size_t>(natom) * 3);
    std::vector<double> omp_doubles(static_cast<std::size_t>(natom) * 3);

    std::cout << "kernel,threads,serial_ms,omp_ms,speedup,efficiency,max_abs_diff,checksum\n";

    {
        update_pos_serial(data, serial_vec);
        update_pos_omp(data, omp_vec);
        const double serial_ms = time_repeated(repeat, [&] { update_pos_serial(data, serial_vec); });
        const double omp_ms = time_repeated(repeat, [&] { update_pos_omp(data, omp_vec); });
        print_measurement({"md_update_pos", threads, serial_ms, omp_ms, max_abs_diff(serial_vec, omp_vec), checksum(omp_vec)});
    }
    {
        std::vector<Vec3> serial_vel = data.base_vel;
        std::vector<Vec3> omp_vel = data.base_vel;
        update_vel_serial(data, serial_vel);
        update_vel_omp(data, omp_vel);
        const double serial_ms = time_repeated(repeat, [&] { update_vel_serial(data, serial_vel); });
        const double omp_ms = time_repeated(repeat, [&] { update_vel_omp(data, omp_vel); });
        print_measurement({"md_update_vel", threads, serial_ms, omp_ms, max_abs_diff(serial_vel, omp_vel), checksum(omp_vel)});
    }
    {
        const double serial_value = kinetic_serial(data);
        const double omp_value = kinetic_omp(data);
        const double serial_ms = time_repeated(repeat, [&] { global_sink += kinetic_serial(data); });
        const double omp_ms = time_repeated(repeat, [&] { global_sink += kinetic_omp(data); });
        print_measurement({"md_kinetic_energy", threads, serial_ms, omp_ms, std::abs(serial_value - omp_value), omp_value});
    }
    {
        const auto serial_value = temp_vector_serial(data);
        const auto omp_value = temp_vector_omp(data);
        const double serial_ms = time_repeated(repeat, [&] { global_sink += checksum(temp_vector_serial(data)); });
        const double omp_ms = time_repeated(repeat, [&] { global_sink += checksum(temp_vector_omp(data)); });
        print_measurement({"md_temp_vector", threads, serial_ms, omp_ms, max_abs_diff(serial_value, omp_value), checksum(omp_value)});
    }
    {
        force_copy_serial(data, serial_vec);
        force_copy_omp(data, omp_vec);
        const double serial_ms = time_repeated(repeat, [&] { force_copy_serial(data, serial_vec); });
        const double omp_ms = time_repeated(repeat, [&] { force_copy_omp(data, omp_vec); });
        print_measurement({"md_force_copy", threads, serial_ms, omp_ms, max_abs_diff(serial_vec, omp_vec), checksum(omp_vec)});
    }
    {
        nep_coord_serial(data, serial_doubles);
        nep_coord_omp(data, omp_doubles);
        const double serial_ms = time_repeated(repeat, [&] { nep_coord_serial(data, serial_doubles); });
        const double omp_ms = time_repeated(repeat, [&] { nep_coord_omp(data, omp_doubles); });
        print_measurement({"nep_coord_fill", threads, serial_ms, omp_ms, max_abs_diff(serial_doubles, omp_doubles), checksum(omp_doubles)});
    }
    {
        const double serial_value = nep_energy_serial(data);
        const double omp_value = nep_energy_omp(data);
        const double serial_ms = time_repeated(repeat, [&] { global_sink += nep_energy_serial(data); });
        const double omp_ms = time_repeated(repeat, [&] { global_sink += nep_energy_omp(data); });
        print_measurement({"nep_energy_sum", threads, serial_ms, omp_ms, std::abs(serial_value - omp_value), omp_value});
    }
    {
        nep_force_serial(data, serial_doubles);
        nep_force_omp(data, omp_doubles);
        const double serial_ms = time_repeated(repeat, [&] { nep_force_serial(data, serial_doubles); });
        const double omp_ms = time_repeated(repeat, [&] { nep_force_omp(data, omp_doubles); });
        print_measurement({"nep_force_fill", threads, serial_ms, omp_ms, max_abs_diff(serial_doubles, omp_doubles), checksum(omp_doubles)});
    }
    {
        const auto serial_value = nep_virial_serial(data);
        const auto omp_value = nep_virial_omp(data);
        const double serial_ms = time_repeated(repeat, [&] { global_sink += checksum(nep_virial_serial(data)); });
        const double omp_ms = time_repeated(repeat, [&] { global_sink += checksum(nep_virial_omp(data)); });
        print_measurement({"nep_virial_sum", threads, serial_ms, omp_ms, max_abs_diff(serial_value, omp_value), checksum(omp_value)});
    }

    if (global_sink == -1.0)
    {
        std::cerr << "unreachable sink: " << global_sink << '\n';
    }

    return 0;
}

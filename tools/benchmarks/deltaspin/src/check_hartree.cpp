#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/potentials/H_Hartree_pw.h"
#include "module_base/parallel_comm.h"
#include "module_base/constants.h"

#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    POOL_WORLD = MPI_COMM_WORLD;
    GlobalV::NPROC_IN_POOL = ranks;
    GlobalV::RANK_IN_POOL = rank;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    double error = 0, mean_error = 0, analytic_error = 0;
    for (bool oblique : {false, true})
    {
        for (int nx : {12, 19})
        {
            ModulePW::PW_Basis cpu("cpu", "double"), gpu("gpu", "double");
            const ModuleBase::Matrix3 lattice(1, oblique ? 0.1 : 0, 0, 0, 1, oblique ? 0.2 : 0, 0, 0, 1);
            for (auto* basis : {&cpu, &gpu})
            {
                basis->initmpi(ranks, rank, MPI_COMM_WORLD);
                basis->initgrids(10.0, lattice, nx, nx + 2, nx + 4);
                basis->initparameters(false, 4.0);
                basis->setuptransform();
                basis->collect_local_pw();
            }
            UnitCell cell;
            cell.omega = 1000;
            cell.tpiba = 2.0 * std::acos(-1.0) / 10.0;
            cell.tpiba2 = cell.tpiba * cell.tpiba;
            const double kernel = ModuleBase::e2 * ModuleBase::FOUR_PI / cell.tpiba2;
            constexpr double cosine_amplitude = 0.01, sine_amplitude = 0.005;
            const double analytic_energy = cell.omega * kernel
                                           * (cosine_amplitude * cosine_amplitude + sine_amplitude * sine_amplitude) / 4.0;
            for (int nspin : {1, 2, 4})
            {
                std::vector<std::vector<double>> density(nspin, std::vector<double>(cpu.nrxx));
                std::vector<const double*> pointers(nspin);
                for (double constant : {-0.3, 0.3, 0.8})
                {
                    for (int ir = 0; ir < cpu.nrxx; ++ir)
                    {
                        const int ix = ir / (cpu.nplane * cpu.ny);
                        const int iy = ir / cpu.nplane % cpu.ny;
                        const int iz = ir % cpu.nplane + cpu.startz[rank];
                        const int global_index = (ix * cpu.ny + iy) * cpu.nz + iz;
                        const double phase = 2.0 * std::acos(-1.0) * ix / nx;
                        const double value = constant + (oblique ? 0.1 * std::sin(0.17 * global_index)
                            : cosine_amplitude * std::cos(phase) + sine_amplitude * std::sin(phase)
                              + 0.007 * std::cos((nx / 2 - 1) * phase));
                        density[0][ir] = nspin == 2 ? 0.7 * value : value;
                        if (nspin == 2) { density[1][ir] = 0.3 * value; }
                        if (nspin == 4)
                        {
                            for (int spin = 1; spin < 4; ++spin) { density[spin][ir] = 0.2 * std::cos(spin * phase); }
                        }
                    }
                    for (int spin = 0; spin < nspin; ++spin) { pointers[spin] = density[spin].data(); }
                    const auto saved = density;
                    const auto expected = elecstate::H_Hartree_pw::v_hartree(cell, &cpu, nspin, pointers.data());
                    const double expected_energy = elecstate::H_Hartree_pw::hartree_energy;
                    const auto actual = elecstate::H_Hartree_pw::v_hartree(cell, &gpu, nspin, pointers.data());
                    const double actual_energy = elecstate::H_Hartree_pw::hartree_energy;
                    if (density != saved || !std::isfinite(actual_energy) || actual_energy < 0)
                    {
                        MPI_Abort(MPI_COMM_WORLD, 2);
                    }
                    error = std::max(error, std::abs(expected_energy - actual_energy) / (1 + std::abs(expected_energy)));
                    long double sum = 0;
                    for (int spin = 0; spin < nspin; ++spin)
                    {
                        for (int ir = 0; ir < cpu.nrxx; ++ir)
                        {
                            if (!std::isfinite(actual(spin, ir))) { MPI_Abort(MPI_COMM_WORLD, 3); }
                            error = std::max(error, std::abs(actual(spin, ir) - expected(spin, ir))
                                                     / (1 + std::abs(expected(spin, ir))));
                            if (spin == 0) { sum += actual(spin, ir); }
                            if (!oblique)
                            {
                                const int ix = ir / (cpu.nplane * cpu.ny);
                                const double phase = 2.0 * std::acos(-1.0) * ix / nx;
                                const double analytical = nspin == 4 && spin != 0 ? 0.0
                                    : kernel * (cosine_amplitude * std::cos(phase) + sine_amplitude * std::sin(phase));
                                analytic_error = std::max(analytic_error, std::abs(actual(spin, ir) - analytical));
                            }
                        }
                    }
                    double global_sum = double(sum);
                    MPI_Allreduce(MPI_IN_PLACE, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    mean_error = std::max(mean_error, std::abs(global_sum) / cpu.nxyz);
                    if (!oblique) { analytic_error = std::max(analytic_error, std::abs(actual_energy - analytic_energy)); }
                }
            }
        }
    }
    double local[3] = {error, mean_error, analytic_error}, global[3];
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "cpu_gpu_relative_error=" << global[0] << " mean_potential_error=" << global[1]
                  << " analytic_error=" << global[2] << '\n';
    }
    if (std::max({global[0], global[1], global[2]}) > 1e-12) { MPI_Abort(MPI_COMM_WORLD, 4); }
    MPI_Finalize();
}

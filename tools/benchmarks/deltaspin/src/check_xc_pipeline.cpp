#include "module_basis/module_pw/pw_basis.h"
#include "module_base/parallel_comm.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int ranks = 0, rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    POOL_WORLD = MPI_COMM_WORLD;
    GlobalV::NPROC_IN_POOL = ranks;
    GlobalV::RANK_IN_POOL = rank;
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    GlobalV::NSPIN = 4;
    GlobalV::DOMAG = true;
    XC_Functional::set_xc_type("PBE");
    XC_Functional::gga_grad = 2;
    double max_potential = 0.0, max_energy = 0.0, max_stress = 0.0, max_reciprocal = 0.0;
    for (bool gamma : {false, true})
    {
        for (int nx : {12, 19})
        {
            ModulePW::PW_Basis cpu("cpu", "double"), gpu("gpu", "double");
            ModuleBase::Matrix3 lattice(1, 0.1, 0, 0, 1, 0.2, 0, 0, 1);
            for (auto* basis : {&cpu, &gpu})
            {
                basis->initmpi(ranks, rank, MPI_COMM_WORLD);
                basis->initgrids(8.0, lattice, nx, nx + 2, nx + 4);
                basis->initparameters(gamma, 10.0);
                basis->setuptransform();
                basis->collect_local_pw();
            }
            Charge cpu_charge, gpu_charge;
            cpu_charge.set_rhopw(&cpu);
            gpu_charge.set_rhopw(&gpu);
            cpu_charge.allocate(4);
            gpu_charge.allocate(4);
            UnitCell cell;
            cell.tpiba = 2.0 * std::acos(-1.0) / 8.0;
            cell.omega = 512.0;
            std::mt19937 generator(53 + rank);
            std::uniform_real_distribution<double> uniform(-1.0, 1.0);
            for (int mode = 0; mode < 6; ++mode)
            {
                ModuleBase::matrix cpu_v(4, cpu.nrxx), gpu_v(4, gpu.nrxx);
                for (int ir = 0; ir < cpu.nrxx; ++ir)
                {
                    const double density = mode == 2 ? 1e-6 * uniform(generator)
                                           : mode == 4 ? 0.0
                                           : mode == 5 ? (ir % 2 ? 1.01e-10 : -0.99e-10)
                                           : 0.5 + 0.2 * uniform(generator);
                    cpu_charge.rho[0][ir] = gpu_charge.rho[0][ir] = density;
                    for (int spin = 1; spin < 4; ++spin)
                    {
                        const double amplitude = mode == 0 ? 0.08 : mode == 1 ? 0.0 : mode == 2 ? 1e-6
                                                 : mode == 5 ? 2e-10 : 1e-12;
                        cpu_charge.rho[spin][ir] = gpu_charge.rho[spin][ir] = amplitude * uniform(generator);
                    }
                    cpu_charge.rho_core[ir] = gpu_charge.rho_core[ir] = 0.1 * std::abs(density);
                    for (int spin = 0; spin < 4; ++spin)
                    {
                        cpu_v(spin, ir) = gpu_v(spin, ir) = 0.02 * uniform(generator);
                    }
                }
                double cpu_e = 1.2, gpu_e = 1.2, cpu_vt = 0.8, gpu_vt = 0.8;
                std::vector<double> cpu_stress, gpu_stress;
                XC_Functional::gradcorr(cpu_e, cpu_vt, cpu_v, &cpu_charge, &cpu, &cell, cpu_stress, false);
                XC_Functional::gradcorr(gpu_e, gpu_vt, gpu_v, &gpu_charge, &gpu, &cell, gpu_stress, false);
                max_energy = std::max({max_energy, std::abs(cpu_e - gpu_e) / (1.0 + std::abs(cpu_e)),
                                       std::abs(cpu_vt - gpu_vt) / (1.0 + std::abs(cpu_vt))});
                for (int ir = 0; ir < cpu.nrxx; ++ir)
                {
                    for (int spin = 0; spin < 4; ++spin)
                    {
                        if (!std::isfinite(gpu_v(spin, ir))) { MPI_Abort(MPI_COMM_WORLD, 2); }
                        max_potential = std::max(max_potential, std::abs(cpu_v(spin, ir) - gpu_v(spin, ir))
                                                               / (1.0 + std::abs(cpu_v(spin, ir))));
                    }
                }
                for (int ig = 0; ig < cpu.npw; ++ig)
                {
                    max_reciprocal = std::max({max_reciprocal, std::abs(cpu_charge.rhog[0][ig] - gpu_charge.rhog[0][ig]),
                                              std::abs(cpu_charge.rhog_core[ig] - gpu_charge.rhog_core[ig])});
                }
                XC_Functional::gradcorr(cpu_e, cpu_vt, cpu_v, &cpu_charge, &cpu, &cell, cpu_stress, true);
                XC_Functional::gradcorr(gpu_e, gpu_vt, gpu_v, &gpu_charge, &gpu, &cell, gpu_stress, true);
                for (int i = 0; i < 9; ++i)
                {
                    max_stress = std::max(max_stress, std::abs(cpu_stress[i] - gpu_stress[i]) / (1.0 + std::abs(cpu_stress[i])));
                }
                auto cpu_total = XC_Functional::v_xc(cpu.nrxx, &cpu_charge, &cell);
                auto gpu_total = XC_Functional::v_xc(gpu.nrxx, &gpu_charge, &cell);
                if (!std::isfinite(std::get<0>(gpu_total)) || !std::isfinite(std::get<1>(gpu_total)))
                {
                    MPI_Abort(MPI_COMM_WORLD, 3);
                }
                max_energy = std::max({max_energy,
                    std::abs(std::get<0>(cpu_total) - std::get<0>(gpu_total)) / (1.0 + std::abs(std::get<0>(cpu_total))),
                    std::abs(std::get<1>(cpu_total) - std::get<1>(gpu_total)) / (1.0 + std::abs(std::get<1>(cpu_total)))});
                for (int spin = 0; spin < 4; ++spin)
                {
                    for (int ir = 0; ir < cpu.nrxx; ++ir)
                    {
                        const double expected = std::get<2>(cpu_total)(spin, ir);
                        const double actual = std::get<2>(gpu_total)(spin, ir);
                        if (!std::isfinite(actual)) { MPI_Abort(MPI_COMM_WORLD, 4); }
                        max_potential = std::max(max_potential, std::abs(actual - expected) / (1.0 + std::abs(expected)));
                    }
                }
            }
        }
    }
    double local[4] = {max_potential, max_energy, max_stress, max_reciprocal}, errors[4] = {};
    MPI_Allreduce(local, errors, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
    {
        std::cout << "potential_relative_error=" << errors[0] << " energy_relative_error=" << errors[1]
                  << " stress_relative_error=" << errors[2] << " reciprocal_error=" << errors[3] << '\n';
    }
    fftw_cleanup_threads();
    MPI_Finalize();
    return *std::max_element(errors, errors + 4) < 1e-10 ? 0 : 3;
}

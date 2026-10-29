#include "module_elecstate/module_charge/charge_mixing.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

#include <mpi.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

void run(int spin, bool domag_z, bool tau, int length, int repeats)
{
    GlobalV::NSPIN = spin;
    GlobalV::DOMAG_Z = domag_z;
    GlobalV::DOMAG = spin == 4 && !domag_z;
    GlobalV::SCF_THR_TYPE = 2;
    GlobalV::MIXING_ANGLE = 0;
    GlobalV::double_grid = false;
    GlobalV::use_paw = false;
    XC_Functional::set_xc_type(tau ? "SCAN" : "PBE");
    ModulePW::PW_Basis basis("cpu", "double");
    basis.nrxx = basis.nxyz = length;
    basis.npw = 1;
    Charge charge;
    charge.set_rhopw(&basis);
    charge.allocate(spin);
    Charge_Mixing mixing;
    mixing.set_rhopw(&basis, &basis);
    mixing.set_mixing("plain", 0.25, 4, 0.0, tau, 0.25, 0.0, 0.0, 0.0, false);
    mixing.init_mixing();
    std::vector<double> input(spin * length), output(spin * length), kin_input, kin_output;
    if (tau) { kin_input.resize(spin * length); kin_output.resize(spin * length); }
    std::vector<double> timings;
    double error = 0;
    for (int repeat = 0; repeat < repeats; ++repeat)
    {
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < spin * length; ++i)
        {
            input[i] = 1.1 + 0.003 * (i % 41) + 0.01 * repeat;
            output[i] = 0.9 + 0.005 * (i % 29) - 0.02 * repeat;
            charge.rho_save[0][i] = input[i];
            charge.rho[0][i] = output[i];
            if (tau)
            {
                kin_input[i] = 0.7 * input[i];
                kin_output[i] = 0.8 * output[i];
                charge.kin_r_save[0][i] = kin_input[i];
                charge.kin_r[0][i] = kin_output[i];
            }
        }
        const auto start = std::chrono::steady_clock::now();
        mixing.mix_rho(&charge);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (repeat >= 2) { timings.push_back(seconds); }
        for (int i = 0; i < spin * length; ++i)
        {
            const int channel = i / length;
            const bool saved = channel == 0 || channel == 3 || !domag_z;
            if (charge.rho_save[0][i] != (saved ? output[i] : input[i])) { std::exit(2); }
            const double expected = input[i] + 0.25 * (output[i] - input[i]);
            error = std::max(error, std::abs(charge.rho[0][i] - expected));
            if (!std::isfinite(charge.rho[0][i])) { std::exit(3); }
            if (tau)
            {
                if (charge.kin_r_save[0][i] != kin_output[i]) { std::exit(4); }
                const double expected_kin = kin_input[i] + 0.25 * (kin_output[i] - kin_input[i]);
                error = std::max(error, std::abs(charge.kin_r[0][i] - expected_kin));
                if (!std::isfinite(charge.kin_r[0][i])) { std::exit(5); }
            }
        }
    }
    if (error > 1e-12) { std::cerr << "Mixing error " << error << '\n'; std::exit(6); }
    std::sort(timings.begin(), timings.end());
    std::cout << "spin=" << spin << " domag_z=" << domag_z << " tau=" << tau << " length=" << length
              << " max_error=" << error << " median_seconds=" << timings[timings.size() / 2] << std::endl;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    if (argc > 1) { run(4, false, false, 8000000, 7); }
    else
    {
        for (int spin : {1, 2, 4})
        {
            for (bool tau : {false, true})
            {
                if (spin != 4 || !tau) { run(spin, false, tau, 101, 3); }
            }
        }
        run(4, true, false, 101, 3);
        run(4, false, false, 1048593, 3);
    }
    MPI_Finalize();
}

#include "lr_state_io.h"

#include "lr_util_print.h"
#include "source_base/tool_quit.h"

namespace LR
{
template <typename T>
void write_excitation_states(
    const std::string &directory,
    const std::string &label,
    const std::string &suffix,
    int rank,
    const double *energies,
    const T *amplitudes,
    int local_size,
    int state_count)
{
    if (rank == 0)
    {
        const std::string filename = directory + "Excitation_Energy_" + label + suffix + ".dat";
        const int written = LR_Util::write_value(
            filename,
            17,
            energies,
            state_count);
        if (written != state_count)
        {
            ModuleBase::WARNING_QUIT("ESolver_LR", "Incomplete LR excitation-energy output");
        }
    }
    const int output_rank = rank + 1;
    const std::string filename =
        directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(output_rank) + ".dat";
    const int written = LR_Util::write_value(
        filename,
        17,
        amplitudes,
        state_count,
        local_size);
    if (written != state_count * local_size)
    {
        ModuleBase::WARNING_QUIT("ESolver_LR", "Incomplete LR excitation-amplitude output");
    }
}
template void write_excitation_states<double>(
    const std::string &,
    const std::string &,
    const std::string &,
    int,
    const double *,
    const double *,
    int,
    int);
template void write_excitation_states<std::complex<double>>(
    const std::string &,
    const std::string &,
    const std::string &,
    int,
    const double *,
    const std::complex<double> *,
    int,
    int);
} // namespace LR

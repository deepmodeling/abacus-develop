#include "source_io/module_output/td_conserv_io.h"

#include "source_base/global_function.h"

#include <iomanip>

namespace ModuleIO
{
void print_conserv(const pw::TDConserv& values,
                   double npsi0,
                   double nrho0,
                   int step,
                   bool full_gram,
                   const std::string& orth,
                   std::ofstream& log)
{
    if (!log.good())
    {
        return;
    }

    const char* stage = step == 0 ? "initial state" : "propagation";
    const int evolve_step = step + 1;
    const double delta_wave = values.npsi - npsi0;
    const double delta_rho = values.nrho - nrho0;
    const std::ios::fmtflags saved_flags = log.flags();
    const std::streamsize saved_prec = log.precision();
    log << " PW RT-TDDFT conservation: evolution step " << evolve_step << ", " << stage << ", " << orth << '\n';
    log << std::right << std::fixed << std::setprecision(10);
    ModuleBase::GlobalFunc::OUT(log, "Npsi", values.npsi);
    ModuleBase::GlobalFunc::OUT(log, "Nrho", values.nrho);

    log << std::scientific << std::setprecision(6);
    ModuleBase::GlobalFunc::OUT(log, "dNpsi", delta_wave);
    ModuleBase::GlobalFunc::OUT(log, "dNrho", delta_rho);
    if (full_gram)
    {
        if (step == 0)
        {
            ModuleBase::GlobalFunc::OUT(log, "initial_orth_error", values.orth_after);
        }
        else
        {
            ModuleBase::GlobalFunc::OUT(log, "orth_before_step_max", values.orth_before);
            ModuleBase::GlobalFunc::OUT(log, "orth_after_step_max", values.orth_after);
        }
    }
    log.flags(saved_flags);
    log.precision(saved_prec);
}
} // namespace ModuleIO

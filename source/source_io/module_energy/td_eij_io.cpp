#include "source_io/module_energy/td_eij_io.h"

#include <cmath>
#include <iomanip>

namespace ModuleIO
{
void print_eij(const std::complex<double>* eij,
               const int nband,
               const int leading_dimension,
               const EijOutputInfo& output_info,
               std::ostream& ofs)
{
    const std::ios::fmtflags original_flags = ofs.flags();
    const std::streamsize original_precision = ofs.precision();
    ofs << std::right << std::noshowpos;
    ofs << " ----------------------------------------------------------\n";
    ofs << "  Eij (Ry): estep = " << output_info.electronic_step << ", SCF iter = " << output_info.scf_iteration
        << ", k-point = " << output_info.kpoint;
    if (output_info.nspin == 2)
    {
        const char* spin_label = output_info.spin == 1 ? " (up)" : " (down)";
        ofs << ", spin = " << output_info.spin << spin_label;
    }
    else if (output_info.nspin == 4)
    {
        ofs << " (spinor)";
    }
    ofs << '\n' << std::scientific << std::setprecision(8);
    bool printed = false;
    for (int i = 0; i < nband; ++i)
    {
        for (int j = 0; j < nband; ++j)
        {
            const std::complex<double>& value = eij[static_cast<size_t>(j) * leading_dimension + i];
            if (std::abs(value.real()) < output_info.threshold && std::abs(value.imag()) < output_info.threshold)
            {
                continue;
            }
            const char* imaginary_sign = std::signbit(value.imag()) ? " - " : " + ";
            ofs << "  i = " << std::setw(2) << i + 1 << ", j = " << std::setw(2) << j + 1 << ", Eij = " << std::setw(16) << value.real()
                << imaginary_sign << std::abs(value.imag()) << " i\n";
            printed = true;
        }
    }
    if (!printed)
    {
        ofs << "  No elements meet the threshold.\n";
    }
    ofs << " ----------------------------------------------------------" << std::endl;
    ofs.flags(original_flags);
    ofs.precision(original_precision);
}

} // namespace ModuleIO

#include "rpa_lri_band_io.h"

#include "source_base/constants.h"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>

namespace RpaLriDetail
{
void write_band_output(std::ostream& output, const BandOutputData& data)
{
    if (data.nk <= 0 || data.nspin <= 0 || data.nbands <= 0 || data.nbasis <= 0)
    {
        throw std::runtime_error("LibRPA band output has invalid dimensions.");
    }
    const unsigned long long expected = static_cast<unsigned long long>(data.nk)
                                        * static_cast<unsigned long long>(data.nspin)
                                        * static_cast<unsigned long long>(data.nbands);
    if (expected != static_cast<unsigned long long>(data.samples.size()))
    {
        throw std::runtime_error("LibRPA band output has an inconsistent sample count.");
    }
    if (!std::isfinite(data.fermi_ry))
    {
        throw std::runtime_error("LibRPA band output has a non-finite Fermi energy.");
    }
    for (const BandSample& sample: data.samples)
    {
        if (!std::isfinite(sample.occupation) || !std::isfinite(sample.energy_ry))
        {
            throw std::runtime_error("LibRPA band output has a non-finite sample.");
        }
    }
    if (!output.good())
    {
        throw std::runtime_error("LibRPA band output stream is not writable.");
    }

    output << std::fixed << std::setprecision(15);
    output << data.nk << '\n';
    output << data.nspin << '\n';
    output << data.nbands << '\n';
    output << data.nbasis << '\n';
    output << data.fermi_ry / 2.0 << '\n';

    std::size_t sample_index = 0;
    for (int ik = 0; ik < data.nk; ++ik)
    {
        for (int is = 0; is < data.nspin; ++is)
        {
            output << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << '\n';
            for (int ib = 0; ib < data.nbands; ++ib, ++sample_index)
            {
                const BandSample& sample = data.samples[sample_index];
                output << std::setw(5) << ib + 1 << "   " << std::setw(8) << sample.occupation << std::setw(25)
                       << sample.energy_ry / 2.0 << std::setw(25) << sample.energy_ry * ModuleBase::Ry_to_eV << '\n';
            }
        }
    }
    if (!output.good())
    {
        throw std::runtime_error("Failed to write LibRPA band output.");
    }
}
} // namespace RpaLriDetail

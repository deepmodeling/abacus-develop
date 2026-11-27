#ifndef RPA_LRI_BAND_IO_H
#define RPA_LRI_BAND_IO_H

#include <iosfwd>
#include <vector>

namespace RpaLriDetail
{
struct BandSample
{
    double occupation = 0.0;
    double energy_ry = 0.0;
};

struct BandOutputData
{
    int nk = 0;
    int nspin = 0;
    int nbands = 0;
    int nbasis = 0;
    double fermi_ry = 0.0;
    std::vector<BandSample> samples;
};

void write_band_output(std::ostream& output, const BandOutputData& data);
}

#endif

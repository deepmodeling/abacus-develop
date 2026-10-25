#include "../rpa_lri_band_io.h"

#include "source_base/constants.h"

#include <cmath>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

namespace
{
RpaLriDetail::BandOutputData make_band_data()
{
    RpaLriDetail::BandOutputData data;
    data.nk = 2;
    data.nspin = 2;
    data.nbands = 2;
    data.nbasis = 5;
    data.fermi_ry = 3.0;
    for (int ik = 0; ik < data.nk; ++ik)
    {
        for (int is = 0; is < data.nspin; ++is)
        {
            for (int ib = 0; ib < data.nbands; ++ib)
            {
                RpaLriDetail::BandSample sample;
                sample.occupation = static_cast<double>((ik + 1) * 10 + (is + 1) * 2 + ib);
                sample.energy_ry = static_cast<double>(100 * ik + 10 * is + ib + 1);
                data.samples.push_back(sample);
            }
        }
    }
    return data;
}
} // namespace

TEST(RpaLriBandIo, WritesHeaderAndSpinBandRecords)
{
    std::ostringstream output;
    const RpaLriDetail::BandOutputData data = make_band_data();
    RpaLriDetail::write_band_output(output, data);

    std::istringstream input(output.str());
    int value = 0;
    input >> value;
    EXPECT_EQ(value, 2);
    input >> value;
    EXPECT_EQ(value, 2);
    input >> value;
    EXPECT_EQ(value, 2);
    input >> value;
    EXPECT_EQ(value, 5);
    double fermi = 0.0;
    input >> fermi;
    EXPECT_DOUBLE_EQ(fermi, 1.5);

    for (int ik = 0; ik < data.nk; ++ik)
    {
        for (int is = 0; is < data.nspin; ++is)
        {
            int kpoint = 0;
            int spin = 0;
            input >> kpoint >> spin;
            EXPECT_EQ(kpoint, ik + 1);
            EXPECT_EQ(spin, is + 1);
            for (int ib = 0; ib < data.nbands; ++ib)
            {
                int band = 0;
                double occupation = 0.0;
                double energy_ha = 0.0;
                double energy_ev = 0.0;
                input >> band >> occupation >> energy_ha >> energy_ev;
                const RpaLriDetail::BandSample& sample = data.samples[(ik * data.nspin + is) * data.nbands + ib];
                EXPECT_EQ(band, ib + 1);
                EXPECT_DOUBLE_EQ(occupation, sample.occupation);
                EXPECT_DOUBLE_EQ(energy_ha, sample.energy_ry / 2.0);
                EXPECT_DOUBLE_EQ(energy_ev, sample.energy_ry * ModuleBase::Ry_to_eV);
            }
        }
    }
    std::string trailing;
    EXPECT_FALSE(input >> trailing);
}

TEST(RpaLriBandIo, RejectsInvalidDimensionsAndSamples)
{
    RpaLriDetail::BandOutputData data = make_band_data();
    data.nbands = 0;
    EXPECT_THROW(RpaLriDetail::write_band_output(std::cout, data), std::runtime_error);

    data = make_band_data();
    data.samples.pop_back();
    EXPECT_THROW(RpaLriDetail::write_band_output(std::cout, data), std::runtime_error);

    data = make_band_data();
    data.samples[0].energy_ry = std::nan("");
    EXPECT_THROW(RpaLriDetail::write_band_output(std::cout, data), std::runtime_error);
}

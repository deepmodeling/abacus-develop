#include "libxc_r_suffix.h"

#include <stdexcept>

namespace ModuleIO
{

std::vector<std::string> libxc_r_channel_suffixes(const std::string& data_name, const int nspin)
{
    if (nspin != 1 && nspin != 2)
    {
        throw std::invalid_argument("libxc_r_channel_suffixes: nspin must be 1 or 2, got "
                                    + std::to_string(nspin));
    }
    const bool polarized = (nspin == 2);

    if (data_name == "rho" || data_name == "vrho")
    {
        return polarized ? std::vector<std::string>{"s1", "s2"} : std::vector<std::string>{"s1"};
    }
    if (data_name == "exc" || data_name == "amag")
    {
        return {""};
    }
    if (data_name == "sigma" || data_name == "vsigma" || data_name == "v2rho2")
    {
        return polarized ? std::vector<std::string>{"00", "01", "11"} : std::vector<std::string>{"00"};
    }
    if (data_name == "v2rhosigma")
    {
        return polarized ? std::vector<std::string>{"0_00", "0_01", "0_11", "1_00", "1_01", "1_11"}
                         : std::vector<std::string>{"0_00"};
    }
    if (data_name == "v2sigma2")
    {
        return polarized ? std::vector<std::string>{"00_00", "00_01", "00_11", "01_01", "01_11", "11_11"}
                         : std::vector<std::string>{"00_00"};
    }
    if (data_name == "v3rho3")
    {
        return polarized ? std::vector<std::string>{"0_0_0", "0_0_1", "0_1_1", "1_1_1"}
                         : std::vector<std::string>{"0_0_0"};
    }
    if (data_name == "v4rho4")
    {
        return polarized ? std::vector<std::string>{"0_0_0_0", "0_0_0_1", "0_0_1_1", "0_1_1_1", "1_1_1_1"}
                         : std::vector<std::string>{"0_0_0_0"};
    }
    if (data_name == "v3rho2sigma")
    {
        return polarized ? std::vector<std::string>{"00_00", "00_01", "00_11",
                                                    "01_00", "01_01", "01_11",
                                                    "11_00", "11_01", "11_11"}
                         : std::vector<std::string>{"00_00"};
    }
    if (data_name == "v3rhosigma2")
    {
        return polarized ? std::vector<std::string>{"0_00_00", "0_00_01", "0_00_11",
                                                    "0_01_01", "0_01_11", "0_11_11",
                                                    "1_00_00", "1_00_01", "1_00_11",
                                                    "1_01_01", "1_01_11", "1_11_11"}
                         : std::vector<std::string>{"0_00_00"};
    }
    if (data_name == "v3sigma3")
    {
        return polarized ? std::vector<std::string>{"00_00_00", "00_00_01", "00_00_11",
                                                    "00_01_01", "00_01_11", "00_11_11",
                                                    "01_01_01", "01_01_11", "01_11_11",
                                                    "11_11_11"}
                         : std::vector<std::string>{"00_00_00"};
    }
    if (data_name == "v4rho3sigma")
    {
        return polarized ? std::vector<std::string>{"000_00", "000_01", "000_11",
                                                    "001_00", "001_01", "001_11",
                                                    "011_00", "011_01", "011_11",
                                                    "111_00", "111_01", "111_11"}
                         : std::vector<std::string>{"000_00"};
    }
    if (data_name == "v4rho2sigma2")
    {
        return polarized ? std::vector<std::string>{"00_00_00", "00_00_01", "00_00_11",
                                                    "00_01_01", "00_01_11", "00_11_11",
                                                    "01_00_00", "01_00_01", "01_00_11",
                                                    "01_01_01", "01_01_11", "01_11_11",
                                                    "11_00_00", "11_00_01", "11_00_11",
                                                    "11_01_01", "11_01_11", "11_11_11"}
                         : std::vector<std::string>{"00_00_00"};
    }
    if (data_name == "v4rhosigma3")
    {
        return polarized ? std::vector<std::string>{"0_00_00_00", "0_00_00_01", "0_00_00_11",
                                                    "0_00_01_01", "0_00_01_11", "0_00_11_11",
                                                    "0_01_01_01", "0_01_01_11", "0_01_11_11",
                                                    "0_11_11_11", "1_00_00_00", "1_00_00_01",
                                                    "1_00_00_11", "1_00_01_01", "1_00_01_11",
                                                    "1_00_11_11", "1_01_01_01", "1_01_01_11",
                                                    "1_01_11_11", "1_11_11_11"}
                         : std::vector<std::string>{"0_00_00_00"};
    }
    if (data_name == "v4sigma4")
    {
        return polarized ? std::vector<std::string>{"00_00_00_00", "00_00_00_01", "00_00_00_11",
                                                    "00_00_01_01", "00_00_01_11", "00_00_11_11",
                                                    "00_01_01_01", "00_01_01_11", "00_01_11_11",
                                                    "00_11_11_11", "01_01_01_01", "01_01_01_11",
                                                    "01_01_11_11", "01_11_11_11", "11_11_11_11"}
                         : std::vector<std::string>{"00_00_00_00"};
    }

    throw std::invalid_argument("libxc_r_channel_suffixes: unknown quantity " + data_name);
}

} // namespace ModuleIO

#include "sep.h"

#include "source_base/global_variable.h"
#include "source_base/parallel_common.h"
#include "source_base/tool_title.h"
#include "source_io/output.h"

#include <fstream>
#include <sstream>
#include <string>

SepPot::SepPot()
{
}

SepPot::~SepPot()
{
    delete[] r;
    r = nullptr;
    delete[] rv;
    rv = nullptr;
}

int SepPot::read_sep(std::ifstream& ifs)
{
    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "Sep.Element")
        {
            iss >> label;
        }
        else if (key == "Sep.XcType")
        {
            iss >> xc_type;
        }
        else if (key == "Sep.Orbital")
        {
            iss >> orbital;
        }
        else if (key == "Sep.Points")
        {
            iss >> mesh;
            delete[] r;
            r = new double[mesh];
            delete[] rv;
            rv = new double[mesh];
        }
        else if (key == "Sep.StripAmount")
        {
            iss >> strip_elec;
        }
        else if (key == "<Sep.Potential")
        {
            double r_val, rv_val;
            int idx = 0;
            while (std::getline(ifs, line) && line != "Sep.Potential>")
            {
                std::istringstream data_line(line);
                if (data_line >> r_val >> rv_val)
                {
                    r[idx] = r_val;
                    rv[idx] = rv_val;
                    idx++;
                }
            }
            break;
        }
    }
    return 0;
}

void SepPot::print_sep_info(std::ofstream& ofs)
{
    ofs << "\n sep_vl:";
    ofs << "\n sep_info:";
    ofs << "\n label         " << label;
    ofs << "\n xc            " << xc_type;
    ofs << "\n orbital       " << orbital;
    ofs << "\n strip electron" << strip_elec;
}

void SepPot::print_sep_vsep(std::ofstream& ofs)
{
    ofs << "\n mesh  " << mesh;
    output::printr1_d(ofs, " r : ", r, mesh);
    output::printr1_d(ofs, " vsep : ", rv, mesh);
    ofs << "\n -----------------------------";
}

#ifdef __MPI

void SepPot::bcast_sep()
{
    ModuleBase::TITLE("SepPot", "bcast_sep");
    Parallel_Common::bcast_bool(is_enable);
    Parallel_Common::bcast_double(r_in);
    Parallel_Common::bcast_double(r_out);
    Parallel_Common::bcast_double(r_power);
    Parallel_Common::bcast_double(enhence_a);
    Parallel_Common::bcast_string(label);
    Parallel_Common::bcast_string(xc_type);
    Parallel_Common::bcast_string(orbital);
    Parallel_Common::bcast_int(strip_elec);
    Parallel_Common::bcast_int(mesh);

    if (GlobalV::MY_RANK != 0 && mesh > 0)
    {
        r = new double[mesh];
        rv = new double[mesh];
    }

    Parallel_Common::bcast_double(r, mesh);
    Parallel_Common::bcast_double(rv, mesh);

    return;
}
#endif // __MPI

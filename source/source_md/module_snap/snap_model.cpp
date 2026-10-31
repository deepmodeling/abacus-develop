#include "snap_model.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

namespace
{
bool next_data_line(std::istream& input, std::string& line)
{
    while (std::getline(input, line))
    {
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos)
        {
            line.erase(comment);
        }
        std::istringstream words(line);
        std::string first;
        if (words >> first)
        {
            return true;
        }
    }
    return false;
}

template <typename T>
bool read_parameter_value(std::istringstream& words, T& value, const std::string& keyword, std::string& error)
{
    if (!(words >> value) || words.fail())
    {
        error = "Invalid SNAP parameter value for keyword: " + keyword;
        return false;
    }
    return true;
}
}

namespace ModuleMD
{
namespace SNAP
{

bool ModelParameters::valid(std::string& error) const
{
    if (!std::isfinite(rcut) || rcut <= 0.0)
    {
        error = "SNAP cutoff must be finite and positive.";
        return false;
    }
    if (!std::isfinite(rcutfac) || rcutfac <= 0.0)
    {
        error = "SNAP rcutfac must be finite and positive.";
        return false;
    }
    if (!std::isfinite(rmin) || rmin < 0.0 || rmin >= rcut)
    {
        error = "SNAP rmin must satisfy 0 <= rmin < rcut.";
        return false;
    }
    if (!std::isfinite(rfac0) || rfac0 <= 0.0)
    {
        error = "SNAP rfac0 must be finite and positive.";
        return false;
    }
    if (twojmax < 0)
    {
        error = "SNAP twojmax must be non-negative.";
        return false;
    }
    if (!std::isfinite(wself))
    {
        error = "SNAP self weight must be finite.";
        return false;
    }
    if (elements.empty() || element_weights.size() != elements.size())
    {
        error = "SNAP element names and weights must be non-empty and have the same size.";
        return false;
    }
    if (element_radii.size() != elements.size())
    {
        error = "SNAP element names and radii must have the same size.";
        return false;
    }
    for (std::size_t i = 0; i < element_weights.size(); ++i)
    {
        if (!std::isfinite(element_weights[i]) || !std::isfinite(element_radii[i]) || element_radii[i] <= 0.0)
        {
            error = "SNAP element weights must be finite.";
            return false;
        }
    }
    return true;
}

bool ModelData::load(const std::string& coefficient_file,
                     const std::string& parameter_file,
                     const std::vector<std::string>& requested_elements,
                     std::string& error)
{
    parameters = ModelParameters();
    parameters.rfac0 = 0.99363;
    coefficients.clear();
    parameters.elements = requested_elements;
    if (parameters.elements.empty())
    {
        error = "SNAP requires at least one requested element.";
        return false;
    }

    std::ifstream coeff_stream(coefficient_file.c_str());
    if (!coeff_stream)
    {
        error = "Cannot open SNAP coefficient file: " + coefficient_file;
        return false;
    }

    std::string line;
    if (!next_data_line(coeff_stream, line))
    {
        error = "SNAP coefficient file is empty.";
        return false;
    }
    std::istringstream header(line);
    int block_count = 0;
    int coefficient_count = 0;
    if (!(header >> block_count >> coefficient_count) || block_count <= 0 || coefficient_count <= 0)
    {
        error = "Invalid SNAP coefficient file header.";
        return false;
    }

    std::vector<bool> found(parameters.elements.size(), false);
    parameters.element_weights.assign(parameters.elements.size(), 0.0);
    parameters.element_radii.assign(parameters.elements.size(), 0.0);
    coefficients.assign(parameters.elements.size(), std::vector<double>());
    for (int block = 0; block < block_count; ++block)
    {
        if (!next_data_line(coeff_stream, line))
        {
            error = "SNAP coefficient file ended before all element blocks were read.";
            return false;
        }
        std::istringstream element_line(line);
        std::string element;
        double radius = 0.0;
        double weight = 0.0;
        if (!(element_line >> element >> radius >> weight))
        {
            error = "Invalid SNAP element line.";
            return false;
        }
        std::size_t element_index = parameters.elements.size();
        for (std::size_t i = 0; i < parameters.elements.size(); ++i)
        {
            if (parameters.elements[i] == element)
            {
                element_index = i;
                break;
            }
        }
        std::vector<double> values(static_cast<std::size_t>(coefficient_count), 0.0);
        for (int coefficient = 0; coefficient < coefficient_count; ++coefficient)
        {
            if (!next_data_line(coeff_stream, line))
            {
                error = "SNAP coefficient file ended before all coefficients were read.";
                return false;
            }
            std::istringstream value_line(line);
            if (!(value_line >> values[static_cast<std::size_t>(coefficient)]))
            {
                error = "Invalid SNAP coefficient value.";
                return false;
            }
        }
        if (element_index < parameters.elements.size())
        {
            if (found[element_index])
            {
                error = "Duplicate SNAP element block: " + element;
                return false;
            }
            found[element_index] = true;
            parameters.element_radii[element_index] = radius;
            parameters.element_weights[element_index] = weight;
            coefficients[element_index] = values;
        }
    }
    for (std::size_t i = 0; i < found.size(); ++i)
    {
        if (!found[i])
        {
            error = "SNAP coefficient file does not contain requested element: " + parameters.elements[i];
            return false;
        }
    }

    std::ifstream parameter_stream(parameter_file.c_str());
    if (!parameter_stream)
    {
        error = "Cannot open SNAP parameter file: " + parameter_file;
        return false;
    }
    bool has_rcutfac = false;
    bool has_twojmax = false;
    while (next_data_line(parameter_stream, line))
    {
        std::istringstream words(line);
        std::string keyword;
        words >> keyword;
        if (keyword == "rcutfac")
        {
            if (!read_parameter_value(words, parameters.rcutfac, keyword, error))
            {
                return false;
            }
            has_rcutfac = true;
        }
        else if (keyword == "twojmax")
        {
            if (!read_parameter_value(words, parameters.twojmax, keyword, error))
            {
                return false;
            }
            has_twojmax = true;
        }
        else if (keyword == "rfac0")
        {
            if (!read_parameter_value(words, parameters.rfac0, keyword, error))
            {
                return false;
            }
        }
        else if (keyword == "rmin0")
        {
            if (!read_parameter_value(words, parameters.rmin, keyword, error))
            {
                return false;
            }
        }
        else if (keyword == "bnormflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.bnorm_flag = value != 0;
        }
        else if (keyword == "switchflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.switch_flag = value != 0;
        }
        else if (keyword == "bzeroflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.bzero_flag = value != 0;
        }
        else if (keyword == "quadraticflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.quadratic_flag = value != 0;
        }
        else if (keyword == "wselfallflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.wselfall_flag = value != 0;
        }
        else if (keyword == "chemflag")
        {
            int value = 0;
            if (!read_parameter_value(words, value, keyword, error))
            {
                return false;
            }
            parameters.chem_flag = value != 0;
        }
    }
    if (!has_rcutfac || !has_twojmax)
    {
        error = "SNAP parameter file must define rcutfac and twojmax.";
        return false;
    }
    double maximum_cutoff = 0.0;
    for (std::size_t i = 0; i < parameters.element_radii.size(); ++i)
    {
        maximum_cutoff = std::max(maximum_cutoff, 2.0 * parameters.element_radii[i] * parameters.rcutfac);
    }
    parameters.rcut = maximum_cutoff;
    return parameters.valid(error);
}

std::size_t ModelData::coefficient_count() const
{
    if (coefficients.empty())
    {
        return 0;
    }
    return coefficients[0].size();
}

} // namespace SNAP
} // namespace ModuleMD

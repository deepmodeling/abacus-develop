#ifndef SNAP_MODEL_H
#define SNAP_MODEL_H

#include <string>
#include <vector>

namespace ModuleMD
{
namespace SNAP
{

struct ModelParameters
{
    double rcut = 0.0;
    double rcutfac = 0.0;
    double rmin = 0.0;
    double rfac0 = 0.0;
    int twojmax = 0;
    bool bnorm_flag = false;
    bool chem_flag = false;
    bool switch_flag = true;
    bool bzero_flag = true;
    bool quadratic_flag = false;
    bool wselfall_flag = false;
    double wself = 1.0;
    std::vector<std::string> elements;
    std::vector<double> element_weights;
    std::vector<double> element_radii;

    bool valid(std::string& error) const;
};

struct ModelData
{
    ModelParameters parameters;
    std::vector<std::vector<double> > coefficients;

    bool load(const std::string& coefficient_file,
              const std::string& parameter_file,
              const std::vector<std::string>& requested_elements,
              std::string& error);

    std::size_t coefficient_count() const;
};

} // namespace SNAP
} // namespace ModuleMD

#endif // SNAP_MODEL_H

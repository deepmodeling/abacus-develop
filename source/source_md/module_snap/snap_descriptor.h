#ifndef SNAP_DESCRIPTOR_H
#define SNAP_DESCRIPTOR_H

#include "snap_cg.h"
#include "snap_model.h"

#include <complex>
#include <vector>

namespace ModuleMD
{
namespace SNAP
{

struct Neighbor
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double weight = 1.0;
    double cutoff_radius = 0.0;
    int element_index = 0;
    double sinner = 0.0;
    double dinner = 0.0;
};

struct UWithDerivative
{
    std::vector<std::complex<double>> derivatives;
};

class DescriptorEvaluator
{
  public:
    explicit DescriptorEvaluator(const ModelParameters &parameters);

    std::vector<std::complex<double>> compute_u(const Neighbor &neighbor) const;

    UWithDerivative compute_u_with_derivative(const Neighbor &neighbor, int direction) const;

    std::vector<std::complex<double>> accumulate_u(
        const std::vector<Neighbor> &neighbors, int central_element) const;

    std::vector<double> compute_bispectrum(const std::vector<std::complex<double>> &u) const;

    std::vector<double> compute_bispectrum_derivative(
        const std::vector<std::complex<double>> &u, const std::vector<std::complex<double>> &du) const;

    int u_size() const { return u_size_; }

    int bispectrum_size() const { return bispectrum_indices_.size(); }

    int channel_count() const { return parameters_.chem_flag ? static_cast<int>(parameters_.elements.size()) : 1; }

  private:
    struct SwitchValue
    {
        double value;
        double derivative;
    };
    struct BispectrumIndex
    {
        int j1;
        int j2;
        int j;
    };

    int u_index(int j, int ma, int mb) const;
    SwitchValue switching(double distance, const Neighbor &neighbor, double cutoff_radius) const;
    void init_root_factors();

    ModelParameters parameters_;
    ClebschGordanTable cg_table_;
    std::vector<double> root_factors_;
    std::vector<int> u_block_;
    std::vector<BispectrumIndex> bispectrum_indices_;
    int u_size_;
};

} // namespace SNAP
} // namespace ModuleMD

#endif // SNAP_DESCRIPTOR_H

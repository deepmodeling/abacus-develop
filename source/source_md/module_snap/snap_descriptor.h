#ifndef SNAP_DESCRIPTOR_H
#define SNAP_DESCRIPTOR_H

#include "snap_clebsch_gordan.h"
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
    std::size_t element_index = 0;
};

struct UWithDerivative
{
    std::vector<std::complex<double> > values;
    std::vector<std::complex<double> > derivatives;
};

class DescriptorEvaluator
{
  public:
    explicit DescriptorEvaluator(const ModelParameters& parameters);

    std::vector<std::complex<double> > compute_u(const Neighbor& neighbor) const;

    UWithDerivative compute_u_with_derivative(const Neighbor& neighbor, int direction) const;

    std::vector<std::complex<double> > accumulate_u(const std::vector<Neighbor>& neighbors) const;

    std::vector<std::complex<double> > accumulate_u(const std::vector<Neighbor>& neighbors,
                                                    std::size_t central_element) const;

    std::vector<double> compute_bispectrum(const std::vector<std::complex<double> >& u) const;

    std::vector<double> compute_bispectrum_derivative(const std::vector<std::complex<double> >& u,
                                                      const std::vector<std::complex<double> >& du) const;

    std::size_t u_size() const
    {
        return u_size_;
    }

    std::size_t bispectrum_size() const
    {
        return bispectrum_indices_.size();
    }

    std::size_t channel_count() const
    {
        return parameters_.chem_flag ? parameters_.elements.size() : 1;
    }

  private:
    struct BispectrumIndex
    {
        int j1;
        int j2;
        int j;
    };

    std::size_t u_index(int j, int ma, int mb) const;
    double cutoff(double distance, double cutoff_radius) const;
    void compute_rootpq();

    ModelParameters parameters_;
    ClebschGordanTable cg_table_;
    std::vector<double> rootpq_;
    std::vector<std::size_t> u_block_;
    std::vector<BispectrumIndex> bispectrum_indices_;
    std::size_t u_size_;
};

} // namespace SNAP
} // namespace ModuleMD

#endif // SNAP_DESCRIPTOR_H

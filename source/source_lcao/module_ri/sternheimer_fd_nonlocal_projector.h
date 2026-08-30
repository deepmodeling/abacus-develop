#ifndef STERNHEIMER_FD_NONLOCAL_PROJECTOR_H
#define STERNHEIMER_FD_NONLOCAL_PROJECTOR_H

#include <complex>
#include <vector>

namespace ModuleRI
{

class SternheimerFDNonlocalProjector
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using Matrix = std::vector<Vector>;

    struct ProjectorBlock
    {
        std::vector<Vector> projectors;
        Matrix d_matrix;
    };

    SternheimerFDNonlocalProjector(int grid_size, double volume_element, std::vector<ProjectorBlock> blocks);

    int grid_size() const;
    double volume_element() const;
    const std::vector<ProjectorBlock>& blocks() const;

    void apply(const Vector& psi, Vector& vpsi) const;
    void add_to(const Vector& psi, Vector& hpsi) const;
    void apply_batch(const Matrix& psi, Matrix& vpsi) const;
    void add_to_batch(const Matrix& psi, Matrix& hpsi) const;

  private:
    void validate_block(const ProjectorBlock& block) const;

    int grid_size_ = 0;
    double volume_element_ = 0.0;
    std::vector<ProjectorBlock> blocks_;
};

} // namespace ModuleRI

#endif

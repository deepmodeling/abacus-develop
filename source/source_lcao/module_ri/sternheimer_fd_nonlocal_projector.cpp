#include "source_lcao/module_ri/sternheimer_fd_nonlocal_projector.h"

#include <stdexcept>
#include <utility>

namespace ModuleRI
{

SternheimerFDNonlocalProjector::SternheimerFDNonlocalProjector(const int grid_size,
                                                               const double volume_element,
                                                               std::vector<ProjectorBlock> blocks)
    : grid_size_(grid_size), volume_element_(volume_element), blocks_(std::move(blocks))
{
    if (grid_size_ <= 0)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector requires a positive grid size.");
    }
    if (volume_element_ <= 0.0)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector requires a positive volume element.");
    }
    for (const ProjectorBlock& block : blocks_)
    {
        validate_block(block);
    }
}

int SternheimerFDNonlocalProjector::grid_size() const
{
    return grid_size_;
}

double SternheimerFDNonlocalProjector::volume_element() const
{
    return volume_element_;
}

const std::vector<SternheimerFDNonlocalProjector::ProjectorBlock>& SternheimerFDNonlocalProjector::blocks() const
{
    return blocks_;
}

void SternheimerFDNonlocalProjector::validate_block(const ProjectorBlock& block) const
{
    const int num_projectors = static_cast<int>(block.projectors.size());
    if (num_projectors <= 0)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector block requires at least one projector.");
    }
    if (static_cast<int>(block.d_matrix.size()) != num_projectors)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector D matrix row count does not match projectors.");
    }

    for (const Vector& projector : block.projectors)
    {
        if (static_cast<int>(projector.size()) != grid_size_)
        {
            throw std::invalid_argument("SternheimerFDNonlocalProjector projector size does not match the grid.");
        }
    }

    for (const Vector& row : block.d_matrix)
    {
        if (static_cast<int>(row.size()) != num_projectors)
        {
            throw std::invalid_argument("SternheimerFDNonlocalProjector D matrix must be square.");
        }
    }
}

void SternheimerFDNonlocalProjector::apply(const Vector& psi, Vector& vpsi) const
{
    if (static_cast<int>(psi.size()) != grid_size_)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector::apply input size does not match the grid.");
    }

    vpsi.assign(grid_size_, Complex(0.0, 0.0));
    add_to(psi, vpsi);
}

void SternheimerFDNonlocalProjector::add_to(const Vector& psi, Vector& hpsi) const
{
    if (static_cast<int>(psi.size()) != grid_size_ || static_cast<int>(hpsi.size()) != grid_size_)
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector::add_to vector sizes do not match the grid.");
    }

    for (const ProjectorBlock& block : blocks_)
    {
        const int num_projectors = static_cast<int>(block.projectors.size());
        Vector coefficients(num_projectors, Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
        for (int ip = 0; ip != num_projectors; ++ip)
        {
            const Vector& beta = block.projectors[ip];
            for (int ir = 0; ir != grid_size_; ++ir)
            {
                coefficients[ip] += volume_element_ * std::conj(beta[ir]) * psi[ir];
            }
        }

        Vector weighted_coefficients(num_projectors, Complex(0.0, 0.0));
        for (int ip = 0; ip != num_projectors; ++ip)
        {
            for (int jp = 0; jp != num_projectors; ++jp)
            {
                weighted_coefficients[ip] += block.d_matrix[ip][jp] * coefficients[jp];
            }
        }

#pragma omp parallel for schedule(static)
        for (int ir = 0; ir != grid_size_; ++ir)
        {
            Complex contribution(0.0, 0.0);
            for (int ip = 0; ip != num_projectors; ++ip)
            {
                contribution += block.projectors[ip][ir] * weighted_coefficients[ip];
            }
            hpsi[ir] += contribution;
        }
    }
}

} // namespace ModuleRI

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
    for (const ProjectorBlock& block: blocks_)
    {
        validate_block(block);
        SparseProjectorBlock sparse_block;
        sparse_block.projectors.reserve(block.projectors.size());
        sparse_block.d_matrix = block.d_matrix;
        for (const Vector& projector: block.projectors)
        {
            SparseProjector sparse_projector;
            for (int ir = 0; ir != grid_size_; ++ir)
            {
                const Complex value = projector[static_cast<std::size_t>(ir)];
                ++projector_value_count_;
                if (value != Complex(0.0, 0.0))
                {
                    sparse_projector.indices.push_back(ir);
                    sparse_projector.values.push_back(value);
                    ++nonzero_projector_value_count_;
                }
            }
            sparse_block.projectors.push_back(std::move(sparse_projector));
        }
        sparse_blocks_.push_back(std::move(sparse_block));
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

std::size_t SternheimerFDNonlocalProjector::projector_value_count() const
{
    return projector_value_count_;
}

std::size_t SternheimerFDNonlocalProjector::nonzero_projector_value_count() const
{
    return nonzero_projector_value_count_;
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

    for (const Vector& projector: block.projectors)
    {
        if (static_cast<int>(projector.size()) != grid_size_)
        {
            throw std::invalid_argument("SternheimerFDNonlocalProjector projector size does not match the grid.");
        }
    }

    for (const Vector& row: block.d_matrix)
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

    for (const SparseProjectorBlock& block: sparse_blocks_)
    {
        const int num_projectors = static_cast<int>(block.projectors.size());
        Vector coefficients(num_projectors, Complex(0.0, 0.0));
#pragma omp parallel for schedule(static)
        for (int ip = 0; ip != num_projectors; ++ip)
        {
            const SparseProjector& beta = block.projectors[static_cast<std::size_t>(ip)];
            for (std::size_t entry = 0; entry != beta.indices.size(); ++entry)
            {
                coefficients[static_cast<std::size_t>(ip)] += volume_element_ * std::conj(beta.values[entry])
                                                              * psi[static_cast<std::size_t>(beta.indices[entry])];
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

        for (int ip = 0; ip != num_projectors; ++ip)
        {
            const SparseProjector& beta = block.projectors[static_cast<std::size_t>(ip)];
            const Complex coefficient = weighted_coefficients[static_cast<std::size_t>(ip)];
            for (std::size_t entry = 0; entry != beta.indices.size(); ++entry)
            {
                hpsi[static_cast<std::size_t>(beta.indices[entry])] += beta.values[entry] * coefficient;
            }
        }
    }
}

void SternheimerFDNonlocalProjector::apply_batch(const Matrix& psi, Matrix& vpsi) const
{
    if (psi.empty())
    {
        vpsi.clear();
        return;
    }
    for (const Vector& column: psi)
    {
        if (static_cast<int>(column.size()) != grid_size_)
        {
            throw std::invalid_argument(
                "SternheimerFDNonlocalProjector::apply_batch input size does not match the grid.");
        }
    }

    vpsi.assign(psi.size(), Vector(static_cast<std::size_t>(grid_size_), Complex(0.0, 0.0)));
    add_to_batch(psi, vpsi);
}

void SternheimerFDNonlocalProjector::add_to_batch(const Matrix& psi, Matrix& hpsi) const
{
    if (psi.size() != hpsi.size())
    {
        throw std::invalid_argument("SternheimerFDNonlocalProjector::add_to_batch column counts do not match.");
    }
    for (std::size_t column = 0; column != psi.size(); ++column)
    {
        if (static_cast<int>(psi[column].size()) != grid_size_ || static_cast<int>(hpsi[column].size()) != grid_size_)
        {
            throw std::invalid_argument(
                "SternheimerFDNonlocalProjector::add_to_batch vector sizes do not match the grid.");
        }
    }

    const std::size_t batch_size = psi.size();
    for (const SparseProjectorBlock& block: sparse_blocks_)
    {
        const int num_projectors = static_cast<int>(block.projectors.size());
        Matrix coefficients(batch_size, Vector(static_cast<std::size_t>(num_projectors), Complex(0.0, 0.0)));
#pragma omp parallel for collapse(2) schedule(static)
        for (int ip = 0; ip != num_projectors; ++ip)
        {
            for (int column = 0; column != static_cast<int>(batch_size); ++column)
            {
                const SparseProjector& beta = block.projectors[static_cast<std::size_t>(ip)];
                Complex coefficient(0.0, 0.0);
                for (std::size_t entry = 0; entry != beta.indices.size(); ++entry)
                {
                    coefficient
                        += volume_element_ * std::conj(beta.values[entry])
                           * psi[static_cast<std::size_t>(column)][static_cast<std::size_t>(beta.indices[entry])];
                }
                coefficients[static_cast<std::size_t>(column)][static_cast<std::size_t>(ip)] = coefficient;
            }
        }

        Matrix weighted_coefficients(batch_size, Vector(static_cast<std::size_t>(num_projectors), Complex(0.0, 0.0)));
#pragma omp parallel for schedule(static)
        for (int column = 0; column != static_cast<int>(batch_size); ++column)
        {
            for (int ip = 0; ip != num_projectors; ++ip)
            {
                for (int jp = 0; jp != num_projectors; ++jp)
                {
                    weighted_coefficients[static_cast<std::size_t>(column)][static_cast<std::size_t>(ip)]
                        += block.d_matrix[static_cast<std::size_t>(ip)][static_cast<std::size_t>(jp)]
                           * coefficients[static_cast<std::size_t>(column)][static_cast<std::size_t>(jp)];
                }
            }
        }

#pragma omp parallel for schedule(static)
        for (int column = 0; column != static_cast<int>(batch_size); ++column)
        {
            for (int ip = 0; ip != num_projectors; ++ip)
            {
                const SparseProjector& beta = block.projectors[static_cast<std::size_t>(ip)];
                const Complex coefficient
                    = weighted_coefficients[static_cast<std::size_t>(column)][static_cast<std::size_t>(ip)];
                for (std::size_t entry = 0; entry != beta.indices.size(); ++entry)
                {
                    hpsi[static_cast<std::size_t>(column)][static_cast<std::size_t>(beta.indices[entry])]
                        += beta.values[entry] * coefficient;
                }
            }
        }
    }
}

} // namespace ModuleRI

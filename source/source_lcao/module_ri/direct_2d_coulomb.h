#ifndef DIRECT_2D_COULOMB_H
#define DIRECT_2D_COULOMB_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace Direct2dCoulomb
{
struct ReciprocalPoint
{
    std::array<int, 3> g_direct{};
    ModuleBase::Vector3<double> q_plus_g_cart;
    double k_bohr = 0.0;
};

struct QuadratureNode
{
    double coordinate_bohr = 0.0;
    double weight_bohr = 0.0;
};

struct ReciprocalLine2D
{
    std::array<int, 2> g_direct{};
    ModuleBase::Vector3<double> q_plus_g_parallel_cart;
    double k_parallel_bohr = 0.0;
    double kz_max_bohr = 0.0;
};

struct MixedFourierPoint
{
    ReciprocalPoint point;
    double kz_weight_bohr = 0.0;
};

struct InPlaneQuadraturePoint
{
    ModuleBase::Vector3<double> q_parallel_cart;
    double q_bohr = 0.0;
    double weight = 0.0;
};

struct RadialTable
{
    int l = 0;
    double dk = 0.0;
    const double* values = nullptr;
    std::size_t size = 0;
};

struct BasisFunction
{
    RadialTable radial;
    int m = 0;
    std::size_t global_index = 0;
    ModuleBase::Vector3<double> tau;
};

struct DenseMatrix
{
    std::size_t dimension = 0;
    std::vector<std::complex<double>> values;

    std::complex<double>& operator()(std::size_t row, std::size_t column)
    {
        return values.at(row * dimension + column);
    }

    const std::complex<double>& operator()(std::size_t row, std::size_t column) const
    {
        return values.at(row * dimension + column);
    }
};

struct Result
{
    DenseMatrix matrix;
    std::size_t regular_g_count = 0;
    double hermitian_residual = 0.0;
    double elapsed_seconds = 0.0;
};

struct MethodMetadata
{
    std::string method = "ewald";
    double ecut_ry = 0.0;
    int kz_order = 0;
    int gamma_order = 8;
    int nq = 0;
    std::size_t naux = 0;
    std::string source_revision;
};

struct MethodMetadata3D
{
    std::string method = "ewald";
    double ecut_ry = 0.0;
    int nq = 0;
    std::size_t naux = 0;
    std::string source_revision;
};

struct Metadata
{
    int iq = 0;
    ModuleBase::Vector3<double> q_direct;
    ModuleBase::Vector3<double> q_cart;
    double ecut_ry = 0.0;
    double alpha = 0.0;
    double chi = 0.0;
    int mpi_owner = 0;
    int mpi_size = 1;
    std::string build_version;
};

ModuleBase::Vector3<double> cartesian_to_direct(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& cartesian);

std::vector<ReciprocalPoint> enumerate_reciprocal_points(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& q_direct,
    double tpiba,
    double ecut_ry);

std::vector<QuadratureNode> gauss_legendre_interval(int order,
                                                    double lower_bohr,
                                                    double upper_bohr);

std::vector<ReciprocalLine2D> enumerate_reciprocal_lines(
    const ModuleBase::Matrix3& reciprocal,
    const ModuleBase::Vector3<double>& q_direct,
    double tpiba,
    double ecut_ry);

std::vector<MixedFourierPoint> expand_mixed_fourier_points(
    const std::vector<ReciprocalLine2D>& lines,
    double tpiba,
    int kz_order);

std::vector<InPlaneQuadraturePoint> gamma_plane_quadrature(
    const ModuleBase::Matrix3& reciprocal,
    double tpiba,
    const std::array<int, 3>& nmp,
    int plane_order,
    double massidda_chi);

std::vector<MixedFourierPoint> expand_gamma_zero_line(
    const std::vector<InPlaneQuadraturePoint>& plane_points,
    double tpiba,
    double ecut_ry,
    int kz_order);

std::vector<std::complex<double>> fourier_coefficients(
    const std::vector<BasisFunction>& basis,
    const ReciprocalPoint& point,
    double omega);

Result build_coulomb_matrix(const std::vector<BasisFunction>& basis,
                            const std::vector<ReciprocalPoint>& points,
                            double omega,
                            double alpha,
                            double chi,
                            bool is_gamma);

Result build_gram_matrix(const std::vector<BasisFunction>& basis,
                               const std::vector<MixedFourierPoint>& points,
                               double omega,
                               double lz_bohr,
                               double alpha);

std::vector<std::size_t> atom_global_indices(
    const std::vector<int>& atom_types,
    const std::vector<std::vector<std::size_t>>& type_traversal_to_local);

void write_matrix_market(const DenseMatrix& matrix, const std::string& path);
void write_reader_v1(const DenseMatrix& matrix,
                     const std::vector<int>& atom_naux,
                     int iq_one_based,
                     int mpi_rank,
                     const std::string& prefix);
std::string format_method_metadata(const MethodMetadata& metadata);
std::string format_3d_method_metadata(const MethodMetadata3D& metadata);
void write_metadata(const Result& result, const Metadata& metadata, const std::string& path);
} // namespace Direct2dCoulomb

#endif

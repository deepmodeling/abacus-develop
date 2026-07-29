#ifndef STERNHEIMER_WAVEFUNCTION_DIAGNOSTIC_H
#define STERNHEIMER_WAVEFUNCTION_DIAGNOSTIC_H

#include <array>
#include <complex>
#include <string>
#include <utility>
#include <vector>

namespace ModuleRI
{

class SternheimerWavefunctionDiagnostic
{
  public:
    using Complex = std::complex<double>;
    using Vector = std::vector<Complex>;
    using NamedVector = std::pair<std::string, Vector>;

    struct Selector
    {
        int iq = 0;
        int ik_full = -1;
        int ib = -1;
        int ifrequency = 0;
        int channel = -1;

        bool matches(int candidate_iq,
                     int candidate_ik_full,
                     int candidate_ib,
                     int candidate_ifrequency,
                     int candidate_channel) const;
    };

    struct Configuration
    {
        Selector selector;
        std::string output_filename;
    };

    struct Metadata
    {
        int nx = 0;
        int ny = 0;
        int nz = 0;
        int iq = 0;
        int ik_full = -1;
        int ib = -1;
        int ifrequency = 0;
        int channel = -1;
        std::array<double, 9> lattice{};
        std::array<double, 3> qpoint{};
        std::array<double, 3> source_kpoint{};
        std::array<double, 3> target_kpoint{};
        double omega_ha = 0.0;
        double omega_ry = 0.0;
        double volume_element = 0.0;
        double reference_eigenvalue_ry = 0.0;
        double weighted_occupation = 0.0;
        double rhs_norm = 0.0;
        double solver_relative_residual = 0.0;
        double equation_relative_residual = 0.0;
        Complex diagonal_branch_element{};
    };

    struct Record
    {
        Metadata metadata;
        std::vector<NamedVector> vectors;
    };

    static void write(const std::string& filename, const Record& record);
    static Record read(const std::string& filename);
    static Configuration parse_configuration(const std::string& specification);
};

} // namespace ModuleRI

#endif

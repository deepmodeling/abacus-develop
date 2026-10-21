#include "source_pw/module_pwdft/dftu_base_io.h"

#include "source_base/global_function.h"
#include "source_base/module_external/lapack_connector.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_pwdft/dftu_base.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <ostream>
#include <vector>

namespace DFTU_BASE
{

struct SocCompactBlocks
{
    std::vector<double> up_up;
    std::vector<double> up_down_real;
    std::vector<double> up_down_imag;
    std::vector<double> down_down;
};

// Retained for the follow-up occupation-eigensystem output refactor. The
// compact writer does not currently call these helpers: this PR only removes
// the legacy text format and separates the writer from the file orchestration.
namespace retained_eigenvalue_helpers
{

inline void jacobi_rotate(std::vector<std::vector<double>>& matrix,
                          const int p,
                          const int q,
                          const int size)
{
    if (std::abs(matrix[p][q]) <= 1.0e-10)
    {
        return;
    }

    const double ratio = (matrix[q][q] - matrix[p][p]) / (2.0 * matrix[p][q]);
    const double tangent = ratio >= 0.0
                               ? 1.0 / (ratio + std::sqrt(1.0 + ratio * ratio))
                               : -1.0 / (-ratio + std::sqrt(1.0 + ratio * ratio));
    const double cosine = 1.0 / std::sqrt(1.0 + tangent * tangent);
    const double sine = tangent * cosine;

    matrix[p][p] -= tangent * matrix[p][q];
    matrix[q][q] += tangent * matrix[p][q];
    matrix[p][q] = matrix[q][p] = 0.0;

    for (int k = 0; k < size; ++k)
    {
        if (k == p || k == q)
        {
            continue;
        }
        const double matrix_kp = cosine * matrix[k][p] - sine * matrix[k][q];
        const double matrix_kq = sine * matrix[k][p] + cosine * matrix[k][q];
        matrix[k][p] = matrix[p][k] = matrix_kp;
        matrix[k][q] = matrix[q][k] = matrix_kq;
    }
}

inline std::vector<double> calculate_eigenvalues(std::vector<std::vector<double>>& matrix,
                                                 const int size)
{
    std::vector<double> eigenvalues(size);
    while (true)
    {
        int p = 0;
        int q = 1;
        for (int row = 0; row < size; ++row)
        {
            for (int column = row + 1; column < size; ++column)
            {
                if (std::abs(matrix[row][column]) > std::abs(matrix[p][q]))
                {
                    p = row;
                    q = column;
                }
            }
        }

        if (std::abs(matrix[p][q]) < 1.0e-10)
        {
            for (int i = 0; i < size; ++i)
            {
                eigenvalues[i] = matrix[i][i];
            }
            break;
        }
        jacobi_rotate(matrix, p, q, size);
    }
    return eigenvalues;
}

inline std::vector<double> calculate_spinor_eigenvalues(const SocCompactBlocks& blocks,
                                                        const int orbital_size)
{
    const int spinor_size = 2 * orbital_size;
    ModuleBase::ComplexMatrix spinor_matrix(spinor_size, spinor_size);
    for (int row = 0; row < orbital_size; ++row)
    {
        for (int column = 0; column < orbital_size; ++column)
        {
            const int index = row * orbital_size + column;
            spinor_matrix(row, column) = blocks.up_up[index];
            spinor_matrix(row, column + orbital_size)
                = std::complex<double>(blocks.up_down_real[index],
                                       blocks.up_down_imag[index]);
            spinor_matrix(row + orbital_size, column)
                = std::complex<double>(blocks.up_down_real[index],
                                       -blocks.up_down_imag[index]);
            spinor_matrix(row + orbital_size, column + orbital_size)
                = blocks.down_down[index];
        }
    }

    std::vector<double> eigenvalues(spinor_size, 0.0);
    std::vector<std::complex<double>> work(1, std::complex<double>(0.0, 0.0));
    std::vector<double> rwork(std::max(1, 3 * spinor_size - 2), 0.0);
    int info = 0;
    int lwork = -1;
    LapackConnector::zheev('N',
                           'U',
                           spinor_size,
                           spinor_matrix,
                           spinor_size,
                           eigenvalues.data(),
                           work.data(),
                           lwork,
                           rwork.data(),
                           &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev workspace query failed");
    }

    lwork = std::max(1, static_cast<int>(work[0].real()));
    work.resize(lwork);
    LapackConnector::zheev('N',
                           'U',
                           spinor_size,
                           spinor_matrix,
                           spinor_size,
                           eigenvalues.data(),
                           work.data(),
                           lwork,
                           rwork.data(),
                           &info);
    if (info != 0)
    {
        ModuleBase::WARNING_QUIT("calculate_spinor_eigenvalues",
                                 "LAPACK zheev diagonalization failed");
    }
    return eigenvalues;
}

} // namespace retained_eigenvalue_helpers

static SocCompactBlocks extract_soc_compact_blocks(const Plus_U_Base& dftu,
                                                   const int iat,
                                                   const int l,
                                                   const OccmatSocLayout layout)
{
    const int size = 2 * l + 1;
    const int block_size = size * size;
    SocCompactBlocks blocks;
    blocks.up_up.resize(block_size, 0.0);
    blocks.up_down_real.resize(block_size, 0.0);
    blocks.up_down_imag.resize(block_size, 0.0);
    blocks.down_down.resize(block_size, 0.0);

    const OccupationMatrix& occmat = dftu.occmat();
    if (layout == SOC_LAYOUT_PAULI)
    {
        const ModuleBase::matrix& occ = occmat.mat(iat, l, 0);
        for (int index = 0; index < block_size; ++index)
        {
            const double charge = occ.c[0 * block_size + index];
            const double spin_x = occ.c[1 * block_size + index];
            const double spin_y = occ.c[2 * block_size + index];
            const double spin_z = occ.c[3 * block_size + index];
            blocks.up_up[index] = 0.5 * (charge + spin_z);
            blocks.up_down_real[index] = 0.5 * spin_x;
            blocks.up_down_imag[index] = 0.5 * spin_y;
            blocks.down_down[index] = 0.5 * (charge - spin_z);
        }
    }
    else
    {
        for (int row = 0; row < size; ++row)
        {
            for (int column = 0; column < size; ++column)
            {
                const int index = row * size + column;
                const double up_up = occmat.get(iat, l, 0, row, column);
                const double up_down = occmat.get(iat, l, 0, row, column + size);
                const double down_up = occmat.get(iat, l, 0, row + size, column);
                const double down_down = occmat.get(iat, l, 0, row + size, column + size);
                blocks.up_up[index] = up_up;
                blocks.up_down_real[index] = 0.5 * (up_down + down_up);
                blocks.down_down[index] = down_down;
            }
        }
    }
    return blocks;
}

static void write_matrix(std::ostream& os, const double* matrix, const int size)
{
    for (int row = 0; row < size; ++row)
    {
        for (int column = 0; column < size; ++column)
        {
            os << " " << std::setw(10) << matrix[row * size + column];
        }
        os << std::endl;
    }
}

static void write_collinear_atom(const Plus_U_Base& dftu,
                                 const UnitCell& ucell,
                                 std::ostream& os,
                                 const int type,
                                 const int iat,
                                 const int l)
{
    const int size = 2 * l + 1;
    double electron_count[2] = {0.0, 0.0};
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int m = 0; m < size; ++m)
        {
            electron_count[spin] += dftu.occmat().get(iat, l, spin, m, m);
        }
    }

    os << "\n " << ucell.atoms[type].label;
    os << " Atom " << iat + 1;
    os << " L " << l;
    os << " mag " << std::fixed << std::setprecision(8)
       << electron_count[0] - electron_count[1] << std::endl;

    os << std::fixed << std::setprecision(8);
    for (int spin = 0; spin < 2; ++spin)
    {
        os << " spin " << spin + 1;
        os << " nelec " << electron_count[spin] << std::endl;
        write_matrix(os, dftu.occmat().mat(iat, l, spin).c, size);
    }
}

static void write_soc_atom(const Plus_U_Base& dftu,
                           const UnitCell& ucell,
                           std::ostream& os,
                           const int type,
                           const int iat,
                           const int l,
                           const OccmatSocLayout layout)
{
    const int size = 2 * l + 1;
    const SocCompactBlocks blocks = extract_soc_compact_blocks(dftu, iat, l, layout);
    double electron_up = 0.0;
    double electron_down = 0.0;
    double magnetization_x = 0.0;
    double magnetization_y = 0.0;
    double magnetization_z = 0.0;
    for (int m = 0; m < size; ++m)
    {
        const int index = m * size + m;
        electron_up += blocks.up_up[index];
        electron_down += blocks.down_down[index];
        magnetization_x += 2.0 * blocks.up_down_real[index];
        magnetization_y += 2.0 * blocks.up_down_imag[index];
        magnetization_z += blocks.up_up[index] - blocks.down_down[index];
    }

    os << "\n " << ucell.atoms[type].label;
    os << " Atom " << iat + 1;
    os << " L " << l;
    os << " mag " << std::fixed << std::setprecision(8)
       << magnetization_x << " " << magnetization_y << " " << magnetization_z
       << std::endl;
    os << std::fixed << std::setprecision(8);

    os << " spin 1 nelec " << electron_up << std::endl;
    write_matrix(os, blocks.up_up.data(), size);
    os << " spin 12 re" << std::endl;
    write_matrix(os, blocks.up_down_real.data(), size);
    os << " spin 12 im" << std::endl;
    write_matrix(os, blocks.up_down_imag.data(), size);
    os << " spin 2 nelec " << electron_down << std::endl;
    write_matrix(os, blocks.down_down.data(), size);
}

void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ostream& os,
                   const OccmatSocLayout soc_layout)
{
    ModuleBase::TITLE("DFTU_BASE", "write_occup_m");
    const int nspin = dftu.occmat().nspin();
    for (int type = 0; type < ucell.ntype; ++type)
    {
        if (!dftu.has_l_channel(type))
        {
            continue;
        }
        const int l = dftu.get_l_channel(type);
        for (int atom = 0; atom < ucell.atoms[type].na; ++atom)
        {
            const int iat = ucell.itia2iat(type, atom);
            if (nspin == 4)
            {
                write_soc_atom(dftu, ucell, os, type, iat, l, soc_layout);
            }
            else
            {
                write_collinear_atom(dftu, ucell, os, type, iat, l);
            }
        }
    }
}

} // namespace DFTU_BASE

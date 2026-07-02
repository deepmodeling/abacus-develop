#include "source_lcao/module_ri/sternheimer_abacus_dft_reference.h"

#include "source_base/matrix.h"
#include "source_estate/elecstate.h"

#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

std::vector<double> copy_matrix_row(const ModuleBase::matrix& values,
                                    const int row_index,
                                    const int num_columns,
                                    const char* quantity_name)
{
    if (row_index < 0 || row_index >= values.nr)
    {
        throw std::invalid_argument(std::string("Sternheimer DFT reference ") + quantity_name
                                    + " k index is out of range.");
    }

    if (num_columns <= 0 || num_columns > values.nc)
    {
        throw std::invalid_argument(std::string("Sternheimer DFT reference ") + quantity_name
                                    + " band count is out of range.");
    }

    std::vector<double> row(num_columns);
    for (int ib = 0; ib != num_columns; ++ib)
    {
        row[ib] = values(row_index, ib);
    }
    return row;
}

} // namespace

std::vector<double> copy_sternheimer_dft_eigenvalues(const elecstate::ElecState& elec_state,
                                                     const int k_index,
                                                     const int num_bands)
{
    return copy_matrix_row(elec_state.ekb, k_index, num_bands, "eigenvalue");
}

std::vector<double> copy_sternheimer_dft_occupations(const elecstate::ElecState& elec_state,
                                                     const int k_index,
                                                     const int num_bands)
{
    return copy_matrix_row(elec_state.wg, k_index, num_bands, "occupation");
}

} // namespace ModuleRI

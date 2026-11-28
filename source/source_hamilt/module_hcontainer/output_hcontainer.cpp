#include "output_hcontainer.h"

#include "source_base/module_out/sparse_matrix.h"

#include <fstream>

namespace hamilt
{

/**
 * @brief Constructor of Output_HContainer
 * @attention ofs should be open outside of this interface
 */
template <typename T>
Output_HContainer<T>::Output_HContainer(
    hamilt::HContainer<T> *hcontainer,
    std::ostream &ofs,
    double sparse_threshold,
    int precision)
    : _hcontainer(hcontainer), _ofs(ofs), _sparse_threshold(sparse_threshold), _precision(precision)
{
    if (this->_sparse_threshold == -1)
    {
        this->_sparse_threshold = 1e-10;
    }
    if (this->_precision == -1)
    {
        this->_precision = 8;
    }
}

template <typename T>
void Output_HContainer<T>::write(bool write_empty)
{
    ModuleIO::write_csr_banner(_ofs);

    int size_for_loop_R = this->_hcontainer->size_R_loop();
    int rx = 0;
    int ry = 0;
    int rz = 0;
    int R_range[2] = {0, 0};
    // find the range of R
    for (int iR = 0; iR < size_for_loop_R; iR++)
    {
        this->_hcontainer->loop_R(iR, rx, ry, rz);
        int max_R = std::max({rx, ry, rz});
        int min_R = std::min({rx, ry, rz});
        if (max_R > R_range[1])
        {
            R_range[1] = max_R;
        }
        if (min_R < R_range[0])
        {
            R_range[0] = min_R;
        }
    }
    // write in order of R
    for (int ix = R_range[0]; ix <= R_range[1]; ix++)
    {
        for (int iy = R_range[0]; iy <= R_range[1]; iy++)
        {
            for (int iz = R_range[0]; iz <= R_range[1]; iz++)
            {
                if (this->_hcontainer->find_R(ix, iy, iz) != -1)
                {
                    this->write_single_R(ix, iy, iz);
                }
                else if (write_empty)
                {
                    _ofs << " " << ix << " " << iy << " " << iz << " 0" << std::endl;
                }
            }
        }
    }
}

template <typename T>
void Output_HContainer<T>::write(int rx_in, int ry_in, int rz_in)
{
    int size_for_loop_R = this->_hcontainer->size_R_loop();
    int rx = 0;
    int ry = 0;
    int rz = 0;
    int find_R = 0;
    for (int iR = 0; iR < size_for_loop_R; iR++)
    {
        this->_hcontainer->loop_R(iR, rx, ry, rz);
        if (rx == rx_in && ry == ry_in && rz == rz_in)
        {
            find_R += 1;
            this->write_single_R(rx, ry, rz);
            break;
        }
    }
    if (find_R == 0)
    {
        ModuleBase::WARNING_QUIT("Output_HContainer::write", "Cannot find the R vector from the HContainer.");
    }
}

template <typename T>
void Output_HContainer<T>::write_single_R(int rx, int ry, int rz)
{
    // The atom-pair traversal below only needs the atom-pair matrices and the
    // global basis size, both of which are available for the flattened
    // BasisFlattenedIndexing layout as well. Reject only when neither layout
    // provides a valid global basis.
    if (this->_hcontainer->get_paraV() == nullptr && this->_hcontainer->get_basis_index() == nullptr)
    {
        ModuleBase::WARNING_QUIT("Output_HContainer::write_single_R",
                                 "neither paraV nor basis_index is set! Unable to write the matrix.");
    }
    this->_hcontainer->fix_R(rx, ry, rz);

    ModuleIO::SparseMatrix<T> sparse_matrix =
        ModuleIO::SparseMatrix<T>(_hcontainer->get_nbasis(), _hcontainer->get_nbasis());
    assert(_hcontainer->get_nbasis() > 0);

    sparse_matrix.setSparseThreshold(this->_sparse_threshold);

    // For the flattened layout the atom-pair window is given in local
    // block-cyclic coordinates; translate each local index to its global
    // position through the distribution map. For the paraV/serial layout the
    // window is already in global coordinates.
    const bool flattened = (this->_hcontainer->get_basis_index() != nullptr);
    const Parallel_Blkcyclic *dist = flattened ? this->_hcontainer->get_matrix_dist() : nullptr;

    for (int iap = 0; iap < this->_hcontainer->size_atom_pairs(); ++iap)
    {
        const hamilt::AtomPair<T>& atom_pair = this->_hcontainer->get_atom_pair(iap);
        const int r_index = atom_pair.find_R(rx, ry, rz);
        if (r_index < 0)
        {
            continue;
        }
        std::tuple<std::vector<int>, T*> tmp_matrix_info = atom_pair.get_matrix_values(r_index);
        int *tmp_index = std::get<0>(tmp_matrix_info).data();
        T *tmp_data = std::get<1>(tmp_matrix_info);
        const int row_begin = tmp_index[0];
        const int row_size = tmp_index[1];
        const int col_begin = tmp_index[2];
        const int col_size = tmp_index[3];
        for (int irow = 0; irow < row_size; ++irow)
        {
            const int g_row = flattened ? dist->irglb(row_begin + irow) : (row_begin + irow);
            for (int icol = 0; icol < col_size; ++icol)
            {
                const int g_col = flattened ? dist->icglb(col_begin + icol) : (col_begin + icol);
                sparse_matrix.insert(g_row, g_col, *tmp_data);
                tmp_data++;
            }
        }
    }

    if (sparse_matrix.getNNZ() != 0)
    {
        _ofs << " " << rx << " " << ry << " " << rz << " " << sparse_matrix.getNNZ() << std::endl;
        sparse_matrix.printToCSR(_ofs, _precision);
    }
    else
    {
        // Write R-block header with 0 nonzero elements and empty CSR data
        // to keep nR in file header consistent with actual R-block count
        _ofs << " " << rx << " " << ry << " " << rz << " 0" << std::endl;
        _ofs << " # CSR values" << std::endl;
        _ofs << std::endl;
        _ofs << " # CSR column indices" << std::endl;
        _ofs << std::endl;
        _ofs << " # CSR row pointers" << std::endl;
        int nbasis = _hcontainer->get_nbasis();
        for (int i = 0; i < nbasis + 1; i++)
        {
            _ofs << " 0";
        }
        _ofs << std::endl;
    }
    this->_hcontainer->unfix_R();
}

template class Output_HContainer<double>;
template class Output_HContainer<std::complex<double>>;

} // namespace hamilt

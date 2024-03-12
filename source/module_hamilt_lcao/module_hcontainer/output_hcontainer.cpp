#include "output_hcontainer.h"

#include <fstream>

#include "module_io/sparse_matrix.h"

namespace hamilt
{

/**
 * @brief Constructor of Output_HContainer
 * @attention ofs should be open outside of this interface
 */
template <typename T>
Output_HContainer<T>::Output_HContainer(hamilt::HContainer<T>* hcontainer,
                                        const Parallel_Orbitals* ParaV,
                                        const UnitCell& ucell,
                                        std::ostream& ofs,
                                        double sparse_threshold,
                                        int precision)
    : _hcontainer(hcontainer),
      _ParaV(ParaV),
      _ucell(ucell),
      _ofs(ofs),
      _sparse_threshold(sparse_threshold),
      _precision(precision)
{
}

template <typename T>
void Output_HContainer<T>::write()
{
    int size_for_loop_R = this->_hcontainer->size_R_loop();
    int rx, ry, rz;
    for (int iR = 0; iR < size_for_loop_R; iR++)
    {
        this->_hcontainer->loop_R(iR, rx, ry, rz);
        this->write_single_R(rx, ry, rz);
    }
}

template <typename T>
void Output_HContainer<T>::write(int rx_in, int ry_in, int rz_in)
{
    int size_for_loop_R = this->_hcontainer->size_R_loop();
    int rx, ry, rz;
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
        ModuleBase::WARNING_QUIT("Output_HContainer::write", "Cannot find the R vector from the HContaine");
    }
}

template <typename T>
void Output_HContainer<T>::write_single_R(int rx, int ry, int rz)
{
    this->_hcontainer->fix_R(rx, ry, rz);
    ModuleIO::SparseMatrix<T> sparse_matrix = ModuleIO::SparseMatrix<T>(this->_ParaV->nrow, this->_ParaV->ncol);
    sparse_matrix.setSparseThreshold(this->_sparse_threshold);
    for (int iap = 0; iap < this->_hcontainer->size_atom_pairs(); ++iap)
    {
        auto atom_pair = this->_hcontainer->get_atom_pair(iap);
        auto tmp_matrix_info = atom_pair.get_matrix_values();
        int* tmp_index = std::get<0>(tmp_matrix_info).data();
        T* tmp_data = std::get<1>(tmp_matrix_info);
        for (int irow = tmp_index[0]; irow < tmp_index[0] + tmp_index[1]; ++irow)
        {
            for (int icol = tmp_index[2]; icol < tmp_index[2] + tmp_index[3]; ++icol)
            {
                sparse_matrix.insert(irow, icol, *tmp_data);
                tmp_data++;
                // to do: consider 2D block-cyclic distribution
            }
        }
    }
    if (sparse_matrix.getNNZ() != 0)
    {
        _ofs << rx << " " << ry << " " << rz << " " << sparse_matrix.getNNZ() << std::endl;
        sparse_matrix.printToCSR(_ofs, _precision);
    }
    this->_hcontainer->unfix_R();
}

// explicit instantiation of template class with double type
template class Output_HContainer<double>;
// to do: explicit instantiation of template class with std::complex<double> type
// template class Output_HContainer<std::complex<double>>;

} // namespace hamilt
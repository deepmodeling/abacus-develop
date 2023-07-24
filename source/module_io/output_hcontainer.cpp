#include "module_io/output_hcontainer.h"

#include <fstream>

#include "module_io/sparse_matrix.h"

namespace ModuleIO
{

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
    int nonzero = 0;
    for (int iap = 0; iap < this->_hcontainer->size_atom_pairs(); ++iap)
    {
        auto atom_pair = this->_hcontainer->get_atom_pair(iap);
        int iat1 = atom_pair.get_atom_i();
        int iat2 = atom_pair.get_atom_j();
        int T1 = _ucell.iat2it[iat1];
        int T2 = _ucell.iat2it[iat2];
        int I1 = _ucell.iat2ia[iat1];
        int I2 = _ucell.iat2ia[iat2];
        const int start1 = _ucell.itiaiw2iwt(T1, I1, 0);
        const int start2 = _ucell.itiaiw2iwt(T2, I2, 0);
        int size1 = _ucell.atoms[T1].nw;
        int size2 = _ucell.atoms[T2].nw;
        for (int iw1 = 0; iw1 < size1; ++iw1)
        {
            const int global_index1 = start1 + iw1;
            for (int iw2 = 0; iw2 < size2; ++iw2)
            {
                const int global_index2 = start2 + iw2;

                T tmp_matrix_value = atom_pair.get_matrix_value(global_index1, global_index2);

                if (std::abs(tmp_matrix_value) > _sparse_threshold)
                {
                    nonzero++;
                    sparse_matrix.addValue(global_index1, global_index2, tmp_matrix_value);
                    // to do: consider 2D block-cyclic distribution
                }
            }
        }
    }
    if (nonzero != 0)
    {
        _ofs << rx << " " << ry << " " << rz << " " << nonzero << std::endl;
        sparse_matrix.convertToCSR(_sparse_threshold);
        sparse_matrix.printCSR(_ofs, _precision);
    }
    this->_hcontainer->unfix_R();
}

// explicit instantiation of template class with double type
template class Output_HContainer<double>;
// to do: explicit instantiation of template class with std::complex<double> type
// template class Output_HContainer<std::complex<double>>;

} // namespace ModuleIO
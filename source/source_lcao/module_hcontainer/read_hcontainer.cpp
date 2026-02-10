#include "read_hcontainer.h"

#include "source_io/module_output/sparse_matrix.h"
#include "source_io/module_output/csr_reader.h"
#include "hcontainer_funcs.h"

#include <fstream>

namespace hamilt
{

/**
 * @brief Constructor of Read_HContainer
 * @attention ifs should be open outside of this interface
 */
template <typename T>
Read_HContainer<T>::Read_HContainer(hamilt::HContainer<T>* hcontainer,
                                        const std::string& filename,
                                        const int nlocal,
                                        const UnitCell* ucell)
    : _hcontainer(hcontainer), _filename(filename), _nlocal(nlocal), _ucell(ucell)
{
}

template <typename T>
void Read_HContainer<T>::read()
{
    // build atom index of col and row
    std::vector<int> atom_index_row;
    std::vector<int> atom_index_col;
    int natom = this->_ucell->nat;
    Parallel_Orbitals pv_serial;
    pv_serial.set_serial(this->_nlocal, this->_nlocal);
    pv_serial.set_atomic_trace(this->_ucell->get_iat2iwt(), this->_ucell->nat, this->_nlocal);
    for (int iat = 0; iat < natom; ++iat)
    {
        int row_size = pv_serial.get_row_size(iat);
        int col_size = pv_serial.get_col_size(iat);
        for (int i = 0; i < row_size; ++i)
        {
            atom_index_row.push_back(iat);
        }
        for (int j = 0; j < col_size; ++j)
        {
            atom_index_col.push_back(iat);
        }
    }
    //
    hamilt::HContainer<T> hcontainer_serial(&pv_serial);

#ifdef __MPI
    if(GlobalV::MY_RANK == 0)
    {
#endif
    ModuleIO::csrFileReader<T> csr(this->_filename);
    int step = csr.getStep();
    int matrix_dimension = csr.getMatrixDimension();
    int r_number = csr.getNumberOfR();

    //construct serial hcontainer firstly
    // prepare atom index mapping from csr row/col to atom index
    for (int i = 0; i < r_number; i++)
    {
        std::vector<int> RCoord = csr.getRCoordinate(i);
        ModuleIO::SparseMatrix<T> sparse_matrix = csr.getMatrix(i);
        for (const auto& element: sparse_matrix.getElements())
        {
            int row = element.first.first;
            int col = element.first.second;
            T value = element.second;
           

            //insert into hcontainer
            int atom_i = atom_index_row[row];
            int atom_j = atom_index_col[col];
            auto* ij_pair = hcontainer_serial.find_pair(atom_i, atom_j);
            if(ij_pair == nullptr)
            {
                //insert new pair
                hamilt::AtomPair<T> new_pair(atom_i, atom_j, RCoord[0], RCoord[1], RCoord[2], &pv_serial);
                hcontainer_serial.insert_pair(new_pair);
            }
            else
            {
                if(ij_pair->find_R(RCoord[0], RCoord[1], RCoord[2]) == -1)
                {
                    //insert new R
                    hamilt::AtomPair<T> new_pair(atom_i, atom_j, RCoord[0], RCoord[1], RCoord[2], &pv_serial);
                    hcontainer_serial.insert_pair(new_pair);
                }
            }
        }
    }
    hcontainer_serial.allocate(nullptr, true);
    // second loop, add values into hcontainer
    for (int i = 0; i < r_number; i++)
    {
        std::vector<int> RCoord = csr.getRCoordinate(i);
        ModuleIO::SparseMatrix<T> sparse_matrix = csr.getMatrix(i);
        for (const auto& element: sparse_matrix.getElements())
        {
            int row = element.first.first;
            int col = element.first.second;
            T value = element.second;

            //insert into hcontainer
            int atom_i = atom_index_row[row];
            int atom_j = atom_index_col[col];
            auto* matrix = hcontainer_serial.find_matrix(atom_i, atom_j, RCoord[0], RCoord[1], RCoord[2]);
            matrix->add_element(row - pv_serial.atom_begin_row[atom_i],
                                col - pv_serial.atom_begin_col[atom_j],
                                value);
        }
    }
#ifdef __MPI
}
    // thirdly, distribute hcontainer_serial to parallel hcontainer
    // send <IJR>s from serial_rank to all ranks
    int my_rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> para_ijrs;
    if (my_rank == 0)
    {
        para_ijrs = hcontainer_serial.get_ijr_info();
        this->_hcontainer->insert_ijrs(&para_ijrs);
        this->_hcontainer->allocate();
    }
    if (my_rank != 0)
    {
        std::vector<int> tmp_ijrs;
        MPI_Status status;
        long tmp_size = 0;
        MPI_Recv(&tmp_size, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
        tmp_ijrs.resize(tmp_size);
        MPI_Recv(tmp_ijrs.data(),
                    tmp_ijrs.size(),
                    MPI_INT,
                    0,
                    1,
                    MPI_COMM_WORLD,
                    &status);
        this->_hcontainer->insert_ijrs(&tmp_ijrs);
        this->_hcontainer->allocate();
    }
    else
    {
        for (int i = 1; i < size; ++i)
        {
            long tmp_size = para_ijrs.size();
            MPI_Send(&tmp_size, 1, MPI_LONG, i, 0, MPI_COMM_WORLD);
            MPI_Send(para_ijrs.data(), para_ijrs.size(), MPI_INT, i, 1, MPI_COMM_WORLD);
        }
    }
    // gather values from serial_rank to Parallels
    transferSerial2Parallels(hcontainer_serial, this->_hcontainer, 0);
#else
    std::vector<int> para_ijrs = hcontainer_serial.get_ijr_info();
    this->_hcontainer->insert_ijrs(&para_ijrs);
    this->_hcontainer->allocate();
    this->_hcontainer->add(hcontainer_serial);
#endif

}

template class Read_HContainer<double>;
template class Read_HContainer<std::complex<double>>;

} // namespace hamilt

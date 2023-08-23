#include "density_matrix.h"

#include "module_base/libm/libm.h"

namespace elecstate
{

//----------------------------------------------------
// density matrix class
//----------------------------------------------------

// destructor
template <typename T>
DensityMatrix<T>::~DensityMatrix()
{
    // do nothing
}

// use unitcell to initialize atom_pairs
template <typename T>
DensityMatrix<T>::DensityMatrix(const K_Vectors* kv_in, const Parallel_Orbitals* paraV_in)
{
    this->_kv = kv_in;
    this->_paraV = paraV_in;
    this->_DMR = new hamilt::HContainer<T>(this->_paraV);
    this->_DMK.reserve(this->_kv->nks);
    ModuleBase::ComplexMatrix zero_DMK(this->_paraV->nrow, this->_paraV->ncol);
    for (int ik = 0; ik < this->_kv->nks; ik++)
    {
        this->_DMK.push_back(zero_DMK);
    }
}

// initialize density matrix DMR from UnitCell (mainly used in UnitTest)
template <typename T>
void DensityMatrix<T>::init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell)
{
    // set up a HContainer
    for (int iat1 = 0; iat1 < ucell->nat; iat1++)
    {
        auto tau1 = ucell->get_tau(iat1);
        int T1, I1;
        ucell->iat2iait(iat1, &I1, &T1);
        AdjacentAtomInfo adjs;
        GridD_in->Find_atom(*ucell, tau1, T1, I1, &adjs);
        // std::cout << "adjs.adj_num: " <<adjs.adj_num << std::endl;
        for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
        {
            const int T2 = adjs.ntype[ad];
            const int I2 = adjs.natom[ad];
            int iat2 = ucell->itia2iat(T2, I2);
            if (this->_paraV->get_row_size(iat1) <= 0 || this->_paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            ModuleBase::Vector3<int>& R_index = adjs.box[ad];
            // std::cout << "R_index: " << R_index.x << " " << R_index.y << " " << R_index.z << std::endl;
            hamilt::AtomPair<T> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, this->_paraV);
            this->_DMR->insert_pair(tmp);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    this->_DMR->allocate(true);
}

// initialize density matrix DMR from another HContainer
template <typename T>
void DensityMatrix<T>::init_DMR(const hamilt::HContainer<T>& DMR_in)
{
    // set up a HContainer
    this->_DMR = new hamilt::HContainer<T>(this->_paraV);
    // synchronize shape
    this->_DMR->shape_synchron(DMR_in);
}

// get _DMR
template <typename T>
hamilt::HContainer<T>* DensityMatrix<T>::get_DMR_pointer()
{
    return this->_DMR;
}

// set _DMK element
template <typename T>
void DensityMatrix<T>::set_DMK(const int ik, const int i, const int j, const T value)
{
    this->_DMK[ik](i, j) = value;
}

// get _DMK nks, nrow, ncol
template <typename T>
int DensityMatrix<T>::get_DMK_nks() const
{
    return this->_DMK.size();
}

template <typename T>
int DensityMatrix<T>::get_DMK_nrow() const
{
    return this->_DMK[0].nr;
}

template <typename T>
int DensityMatrix<T>::get_DMK_ncol() const
{
    return this->_DMK[0].nc;
}

// calculate DMR from DMK
template <typename T>
void DensityMatrix<T>::cal_DMR()
{
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    for (int i = 0; i < this->_DMR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<T>& tmp_ap = this->_DMR->get_atom_pair(i);
        int iat1 = tmp_ap.get_atom_i();
        int iat2 = tmp_ap.get_atom_j();
        // get global indexes of whole matrix for each atom in this process
        int row_ap = this->_paraV->atom_begin_row[iat1];
        int col_ap = this->_paraV->atom_begin_col[iat2];
        if (row_ap == -1 || col_ap == -1)
        {
            throw std::string("Atom-pair not belong this process");
        }
        for (int ir = 0; ir < tmp_ap.get_R_size(); ++ir)
        {
            const int* r_index = tmp_ap.get_R_index(ir);
            hamilt::BaseMatrix<T>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
            if (tmp_matrix == nullptr)
            {
                std::cout << "tmp_matrix is nullptr" << std::endl;
                continue;
            }
#endif
            std::complex<T> tmp_res;
            // loop over k-points
            for (int ik = 0; ik < this->_kv->nks; ++ik)
            {
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                const double arg = (this->_kv->kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                std::complex<double> kphase = std::complex<double>(cosp, sinp);
                // set DMR element
                for (int i = 0; i < this->_paraV->get_row_size(iat1); ++i)
                {
                    for (int j = 0; j < this->_paraV->get_col_size(iat2); ++j)
                    {
                        // since DMK is column-major, we need to transpose it
                        tmp_res = kphase * this->_DMK[ik](row_ap + i, col_ap + j);
                        tmp_matrix->add_element(i, j, tmp_res.real());
                    }
                }
            }
        }
    }
}

// calculate DMR from DMK using blas
template <typename T>
void DensityMatrix<T>::cal_DMR_blas()
{
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    for (int i = 0; i < this->_DMR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<T>& tmp_ap = this->_DMR->get_atom_pair(i);
        int iat1 = tmp_ap.get_atom_i();
        int iat2 = tmp_ap.get_atom_j();
        // get global indexes of whole matrix for each atom in this process
        int row_ap = this->_paraV->atom_begin_row[iat1];
        int col_ap = this->_paraV->atom_begin_col[iat2];
        if (row_ap == -1 || col_ap == -1)
        {
            throw std::string("Atom-pair not belong this process");
        }
        for (int ir = 0; ir < tmp_ap.get_R_size(); ++ir)
        {
            const int* r_index = tmp_ap.get_R_index(ir);
            hamilt::BaseMatrix<T>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
            if (tmp_matrix == nullptr)
            {
                std::cout << "tmp_matrix is nullptr" << std::endl;
                continue;
            }
#endif
            // loop over k-points
            for (int ik = 0; ik < this->_kv->nks; ++ik)
            {
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                const double arg = (this->_kv->kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                std::complex<double> kphase = std::complex<double>(cosp, sinp);
                // set DMR element
                T* tmp_DMR = tmp_matrix->get_pointer();
                std::complex<double>* tmp_DMK = this->_DMK[ik].c;
                T* DMK_real_pointer = nullptr;
                T* DMK_imag_pointer = nullptr;
                tmp_DMK += row_ap * this->_paraV->ncol + col_ap;
                for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
                {
                    DMK_real_pointer = (T*)tmp_DMK;
                    DMK_imag_pointer = DMK_real_pointer + 1;
                    BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                        kphase.real(),
                                        DMK_real_pointer,
                                        2,
                                        tmp_DMR,
                                        1);
                    BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                        kphase.imag(),
                                        DMK_imag_pointer,
                                        2,
                                        tmp_DMR,
                                        1);
                    tmp_DMK += this->_paraV->get_col_size(iat2);
                    tmp_DMR += this->_paraV->get_col_size(iat2);
                }
            }
        }
    }
}

// read *.dmk into density matrix dm(k)
template <typename T>
void DensityMatrix<T>::read_DMK(const std::string directory, const int ik)
{
    // read
    std::string fn;
    fn = directory + std::to_string(ik) + ".dmk";
    //
    bool quit_abacus = false;

    std::ifstream ifs;

    ifs.open(fn.c_str());
    if (!ifs)
    {
        quit_abacus = true;
    }
    else
    {
        // if the number is not match,
        // quit the program or not.
        bool quit = false;

        ModuleBase::CHECK_DOUBLE(ifs, this->_kv->kvec_d[ik].x, quit);
        ModuleBase::CHECK_DOUBLE(ifs, this->_kv->kvec_d[ik].y, quit);
        ModuleBase::CHECK_DOUBLE(ifs, this->_kv->kvec_d[ik].z, quit);
        ModuleBase::CHECK_INT(ifs, this->_paraV->nrow);
        ModuleBase::CHECK_INT(ifs, this->_paraV->ncol);
    } // If file exist, read in data.
    // Finish reading the first part of density matrix.

    for (int i = 0; i < this->_paraV->nrow; ++i)
    {
        for (int j = 0; j < this->_paraV->ncol; ++j)
        {
            ifs >> this->_DMK[ik](i, j);
        }
    }
    ifs.close();
}

// output density matrix dm(k) into *.dmk
template <typename T>
void DensityMatrix<T>::write_DMK(const std::string directory, const int ik)
{

    // write
    std::string fn;
    fn = directory + std::to_string(ik) + ".dmk";
    std::ofstream ofs;
    ofs.open(fn.c_str());
    if (!ofs)
    {
        ModuleBase::WARNING("elecstate::write_dmk", "Can't create DENSITY MATRIX File!");
    }
    ofs << this->_kv->kvec_d[ik].x << " " << this->_kv->kvec_d[ik].y << " " << this->_kv->kvec_d[ik].z << std::endl;
    ofs << "\n  " << this->_paraV->nrow << " " << this->_paraV->ncol << std::endl;

    ofs << std::setprecision(3);
    ofs << std::scientific;

    for (int i = 0; i < this->_paraV->nrow; ++i)
    {
        for (int j = 0; j < this->_paraV->ncol; ++j)
        {
            if (j % 8 == 0)
                ofs << "\n";
            ofs << " " << this->_DMK[ik](i, j).real();
            // ofs << " " << DM[is][i][j];
        }
    }

    ofs.close();
}

// read *.dmk into density matrix dm(k)
template <typename T>
void DensityMatrix<T>::set_DMK_files(const std::string directory)
{
    // read
    for (int ik = 0; ik < this->_kv->nks; ++ik)
    {
        this->read_DMK(directory, ik);
    }
}

// write density matrix dm(k) into *.dmk with all k-points
template <typename T>
void DensityMatrix<T>::output_DMK(const std::string directory)
{
    // write
    for (int ik = 0; ik < this->_kv->nks; ++ik)
    {
        this->write_DMK(directory, ik);
    }
}

// get a matrix element of density matrix dm(k)
template <typename T>
T DensityMatrix<T>::get_DMK(const int ik, const int i, const int j) const
{
    return this->_DMK[ik](i, j).real();
}

// T of HContainer can be double or complex<double>
template class DensityMatrix<double>;

} // namespace elecstate

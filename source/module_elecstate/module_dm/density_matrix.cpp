#include "density_matrix.h"

#include "module_base/libm/libm.h"

namespace elecstate
{

//----------------------------------------------------
// density matrix class
//----------------------------------------------------

// destructor
template <typename TK, typename TR>
DensityMatrix<TK, TR>::~DensityMatrix()
{
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
}

// constructor
template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const K_Vectors* kv_in, const Parallel_Orbitals* paraV_in, const int nspin)
{
    this->_kv = kv_in;
    this->_paraV = paraV_in;
    // set this->_nspin
    if (nspin == 1 || nspin == 4)
    {
        this->_nspin = 1;
    }
    else if (nspin == 2)
    {
        this->_nspin = 2;
#ifdef __DEBUG
        assert(kv_in->nks % 2 == 0);
#endif
    }
    else
    {
        throw std::string("nspin must be 1, 2 or 4");
    }
    // set this->_nks, which is real number of k-points
    this->_nks = kv_in->nks / this->_nspin;
    // reserve memory for _DMK
    this->_DMK.reserve(this->_kv->nks);
    std::vector<TK> tmp_DMK(this->_paraV->nrow * this->_paraV->ncol);
    for (int ik = 0; ik < this->_kv->nks; ik++)
    {
        this->_DMK.push_back(tmp_DMK);
    }
}

// initialize density matrix DMR from UnitCell (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell)
{
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // construct a new DMR
    hamilt::HContainer<TR>* tmp_DMR;
    tmp_DMR = new hamilt::HContainer<TR>(this->_paraV);
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
            hamilt::AtomPair<TR> tmp_ap(iat1, iat2, R_index.x, R_index.y, R_index.z, this->_paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    tmp_DMR->allocate(true);
    this->_DMR.push_back(tmp_DMR);
    // add another DMR if nspin==2
    if (this->_nspin == 2)
    {
        hamilt::HContainer<TR>* tmp_DMR1;
        tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
        this->_DMR.push_back(tmp_DMR1);
    }
}

// initialize density matrix DMR from another HContainer (mainly used)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(const hamilt::HContainer<TR>& DMR_in)
{
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // set up a HContainer using another one
    for (int is = 0; is < this->_nspin; ++is) // loop over spin
    {
        hamilt::HContainer<TR>* tmp_DMR;
        tmp_DMR = new hamilt::HContainer<TR>(DMR_in);
        // zero.out
        tmp_DMR->set_zero();
        this->_DMR.push_back(tmp_DMR);
    }
}

// get _DMR pointer
template <typename TK, typename TR>
hamilt::HContainer<TR>* DensityMatrix<TK, TR>::get_DMR_pointer(const int ispin) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    return this->_DMR[ispin - 1];
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK(const int ispin, const int ik, const int i, const int j, const TK value)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->ncol + j] = value;
}

// get a matrix element of density matrix dm(k)
template <typename TK, typename TR>
TK DensityMatrix<TK, TR>::get_DMK(const int ispin, const int ik, const int i, const int j) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    return this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->ncol + j];
}

// get _DMK nks, nrow, ncol
template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nks() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_kv->nks;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nrow() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_paraV->nrow;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_ncol() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_paraV->ncol;
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR()
{
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nks * (is - 1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = this->_DMR[is - 1];
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        for (int i = 0; i < tmp_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<double>& tmp_ap = tmp_DMR->get_atom_pair(i);
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
                hamilt::BaseMatrix<double>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
                if (tmp_matrix == nullptr)
                {
                    std::cout << "tmp_matrix is nullptr" << std::endl;
                    continue;
                }
#endif
                // loop over k-points
                for (int ik = 0; ik < this->_nks; ++ik)
                {
                    // cal k_phase
                    // if TK==std::complex<double>, kphase is e^{ikR}
                    const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                    const double arg = (this->_kv->kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    std::complex<double> kphase = std::complex<double>(cosp, sinp);
                    // set DMR element
                    double* tmp_DMR_pointer = tmp_matrix->get_pointer();
                    std::complex<double>* tmp_DMK_pointer = this->_DMK[ik + ik_begin].data();
                    double* DMK_real_pointer = nullptr;
                    double* DMK_imag_pointer = nullptr;
                    tmp_DMK_pointer += row_ap * this->_paraV->ncol + col_ap;
                    for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
                    {
                        DMK_real_pointer = (double*)tmp_DMK_pointer;
                        DMK_imag_pointer = DMK_real_pointer + 1;
                        BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                            kphase.real(),
                                            DMK_real_pointer,
                                            2,
                                            tmp_DMR_pointer,
                                            1);
                        BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                            kphase.imag(),
                                            DMK_imag_pointer,
                                            2,
                                            tmp_DMR_pointer,
                                            1);
                        tmp_DMK_pointer += this->_paraV->get_col_size(iat2);
                        tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
                    }
                }
            }
        }
    }
}

// calculate DMR from DMK using blas for gamma-only calculation
template <>
void DensityMatrix<double, double>::cal_DMR()
{
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nks * (is - 1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = this->_DMR[is - 1];
#ifdef __DEBUG
        assert(tmp_DMR->is_gamma_only() == true);
        assert(this->_nks == 1);
#endif
        // #ifdef _OPENMP
        // #pragma omp parallel for
        // #endif
        for (int i = 0; i < tmp_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<double>& tmp_ap = tmp_DMR->get_atom_pair(i);
            int iat1 = tmp_ap.get_atom_i();
            int iat2 = tmp_ap.get_atom_j();
            // get global indexes of whole matrix for each atom in this process
            int row_ap = this->_paraV->atom_begin_row[iat1];
            int col_ap = this->_paraV->atom_begin_col[iat2];
            if (row_ap == -1 || col_ap == -1)
            {
                throw std::string("Atom-pair not belong this process");
            }
            // R index
            const int* r_index = tmp_ap.get_R_index(0);
#ifdef __DEBUG
            assert(tmp_ap.get_R_size() == 1);
            assert(r_index[0] == 0 && r_index[1] == 0 && r_index[2] == 0);
#endif
            hamilt::BaseMatrix<double>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
            if (tmp_matrix == nullptr)
            {
                std::cout << "tmp_matrix is nullptr" << std::endl;
                continue;
            }
#endif
            // k index
            double kphase = 1;
            // set DMR element
            double* tmp_DMR_pointer = tmp_matrix->get_pointer();
            double* tmp_DMK_pointer = this->_DMK[0 + ik_begin].data();
            tmp_DMK_pointer += row_ap * this->_paraV->ncol + col_ap;
            for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
            {
                BlasConnector::axpy(this->_paraV->get_col_size(iat2), kphase, tmp_DMK_pointer, 1, tmp_DMR_pointer, 1);
                tmp_DMK_pointer += this->_paraV->get_col_size(iat2);
                tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
            }
        }
    }
}

// read *.dmk into density matrix dm(k)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::read_DMK(const std::string directory, const int ispin, const int ik)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // read
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
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
            ifs >> this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->ncol + j];
        }
    }
    ifs.close();
}

// output density matrix dm(k) into *.dmk
template <>
void DensityMatrix<double, double>::write_DMK(const std::string directory, const int ispin, const int ik)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
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
            ofs << " " << this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->ncol + j];
        }
    }

    ofs.close();
}

template <>
void DensityMatrix<std::complex<double>, double>::write_DMK(const std::string directory, const int ispin, const int ik)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // write
    std::string fn;
    fn = directory + "SPIN" + std::to_string(ispin) + "_" + std::to_string(ik) + ".dmk";
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
            ofs << " " << this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->ncol + j].real();
        }
    }

    ofs.close();
}

// T of HContainer can be double or complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case

} // namespace elecstate

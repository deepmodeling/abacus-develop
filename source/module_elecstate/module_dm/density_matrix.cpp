#include "density_matrix.h"

#include "module_base/libm/libm.h"
#include "module_base/tool_title.h"
#include "module_base/memory.h"
#include "module_base/timer.h"

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
}

// constructor for multi-k
template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const K_Vectors* kv_in, const Parallel_Orbitals* paraV_in, const int nspin)
{
    ModuleBase::TITLE("DensityMatrix", "DensityMatrix-MK");
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
    // allocate memory for _DMK
    this->_DMK.resize(this->_kv->nks);
    for (int ik = 0; ik < this->_kv->nks; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}

// constructor for Gamma-Only
template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in, const int nspin)
{
    ModuleBase::TITLE("DensityMatrix", "DensityMatrix-GO");
    this->_paraV = paraV_in;
    // set this->_nspin
    if (nspin == 1 || nspin == 4)
    {
        this->_nspin = 1;
    }
    else if (nspin == 2)
    {
        this->_nspin = 2;
    }
    else
    {
        throw std::string("nspin must be 1, 2 or 4");
    }
    // set this->_nks, which is real number of k-points
    this->_nks = 1;
    // allocate memory for _DMK
    this->_DMK.resize(_nspin);
    for (int ik = 0; ik < this->_nspin; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}

// initialize density matrix DMR from UnitCell (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(Grid_Driver* GridD_in, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
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
    if(std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->_DMR.push_back(tmp_DMR);
    // add another DMR if nspin==2
    if (this->_nspin == 2)
    {
        hamilt::HContainer<TR>* tmp_DMR1;
        tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
        this->_DMR.push_back(tmp_DMR1);
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

// initialize density matrix DMR from UnitCell and RA (mainly used in UnitTest)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(Record_adj& ra, const UnitCell* ucell)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
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
        for (int ad = 0; ad < ra.na_each[iat1]; ++ad)
        {
            const int T2 = ra.info[iat1][ad][3];
            const int I2 = ra.info[iat1][ad][4];
            int iat2 = ucell->itia2iat(T2, I2);
            if (this->_paraV->get_row_size(iat1) <= 0 || this->_paraV->get_col_size(iat2) <= 0)
            {
                continue;
            }
            hamilt::AtomPair<TR> tmp_ap(iat1, iat2, ra.info[iat1][ad][0], ra.info[iat1][ad][1], ra.info[iat1][ad][2], this->_paraV);
            tmp_DMR->insert_pair(tmp_ap);
        }
    }
    // allocate the memory of BaseMatrix in SR, and set the new values to zero
    if(std::is_same<TK, double>::value)
    {
        tmp_DMR->fix_gamma();
    }
    tmp_DMR->allocate(nullptr, true);
    this->_DMR.push_back(tmp_DMR);
    // add another DMR if nspin==2
    if (this->_nspin == 2)
    {
        hamilt::HContainer<TR>* tmp_DMR1;
        tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
        this->_DMR.push_back(tmp_DMR1);
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

// initialize density matrix DMR from another HContainer (mainly used)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(const hamilt::HContainer<TR>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
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
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
}

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::init_DMR(const hamilt::HContainer<TRShift>& DMR_in)
{
    ModuleBase::TITLE("DensityMatrix", "init_DMR");
    // ensure _DMR is empty
    for (auto& it: this->_DMR)
    {
        delete it;
    }
    this->_DMR.clear();
    // set up a HContainer using another one
    int size_ap = DMR_in.size_atom_pairs();
    if(size_ap > 0)
    {   
        const Parallel_Orbitals* paraV_ = DMR_in.get_atom_pair(0).get_paraV();
        hamilt::HContainer<TR>* tmp_DMR = new hamilt::HContainer<TR>(paraV_);
        for(int iap=0;iap<size_ap;iap++)
        {
            const int iat1 = DMR_in.get_atom_pair(iap).get_atom_i();
            const int iat2 = DMR_in.get_atom_pair(iap).get_atom_j();
            for(int ir = 0;ir<DMR_in.get_atom_pair(iap).get_R_size();ir++)
            {
                const int* R_index = DMR_in.get_atom_pair(iap).get_R_index(ir);
                hamilt::AtomPair<TR> tmp_ap(iat1, iat2, R_index[0], R_index[1], R_index[2], paraV_);
                tmp_DMR->insert_pair(tmp_ap);
            }
        }
        tmp_DMR->allocate(nullptr, true);
        this->_DMR.push_back(tmp_DMR);
        if(this->_nspin == 2)
        {
            hamilt::HContainer<TR>* tmp_DMR1 = new hamilt::HContainer<TR>(*tmp_DMR);
            this->_DMR.push_back(tmp_DMR1);
        }
    }
    ModuleBase::Memory::record("DensityMatrix::DMR", this->_DMR.size() * this->_DMR[0]->get_memory_size());
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

// get _DMK[ik] pointer
template <typename TK, typename TR>
TK* DensityMatrix<TK, TR>::get_DMK_pointer(const int ik) const
{
#ifdef __DEBUG
    assert(ik < this->_nks * this->_nspin);
#endif
    return const_cast<TK*>(this->_DMK[ik].data());
}

// get _DMK[ik] vector
template <typename TK, typename TR>
std::vector<std::vector<TK>> DensityMatrix<TK, TR>::get_DMK_vector() const
{
    return this->_DMK;
}

// get _paraV pointer
template <typename TK, typename TR>
const Parallel_Orbitals* DensityMatrix<TK, TR>::get_paraV_pointer() const
{
    return this->_paraV;
}

// get _kv pointer
template <typename TK, typename TR>
const K_Vectors* DensityMatrix<TK, TR>::get_kv_pointer() const
{
    return this->_kv;
}

// set DMK using a pointer
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_pointer(const int ik, TK* DMK_in)
{
#ifdef __DEBUG
    assert(ik < this->_nks * this->_nspin);
#endif
    this->_DMK[ik].assign(DMK_in, DMK_in + this->_paraV->nrow * this->_paraV->ncol);
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK(const int ispin, const int ik, const int i, const int j, const TK value)
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
    assert(ik >= 0 && ik < this->_nks);
#endif
    // consider transpose col=>row
    this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->nrow + j] = value;
}

// set _DMK element
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::set_DMK_zero()
{
    for (int ik = 0; ik < _nspin * _nks; ik++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->_DMK[ik].data(), this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
}

// get a matrix element of density matrix dm(k)
template <typename TK, typename TR>
TK DensityMatrix<TK, TR>::get_DMK(const int ispin, const int ik, const int i, const int j) const
{
#ifdef __DEBUG
    assert(ispin > 0 && ispin <= this->_nspin);
#endif
    // consider transpose col=>row
    return this->_DMK[ik + this->_nks * (ispin - 1)][i * this->_paraV->nrow + j];
}

// get _DMK nks, nrow, ncol
template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_nks() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
    assert(this->_kv != nullptr);
#endif
    return this->_kv->nks;
}

template <typename TK, typename TR>
int DensityMatrix<TK, TR>::get_DMK_size() const
{
#ifdef __DEBUG
    assert(this->_DMK.size() != 0);
#endif
    return this->_DMK.size();
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

template <typename TK, typename TR>
void DensityMatrix<TK, TR>::save_DMR()
{
    ModuleBase::TITLE("DensityMatrix", "save_DMR");
    ModuleBase::timer::tick("DensityMatrix", "save_DMR");
    
    const int nnr = this->_DMR[0]->get_nnr();
    // allocate if _DMR_save is empty 
    if(_DMR_save.size() == 0)
    {
        _DMR_save.resize(this->_DMR.size());
    }
    // resize if _DMR_save[is].size is not equal to _DMR.size
    for(int is = 0; is < _DMR_save.size(); is++)
    {
        if(_DMR_save[is].size() != nnr)
        {
            _DMR_save[is].resize(nnr);
        }
    }
    // save _DMR to _DMR_save
    for(int is=0;is<this->_DMR.size();is++)
    {
        TR* DMR_pointer = this->_DMR[is]->get_wrapper();
        TR* DMR_save_pointer = _DMR_save[is].data();
        // set to zero
        ModuleBase::GlobalFunc::ZEROS(DMR_save_pointer, nnr);
        for(int i=0;i<nnr;i++)
        {
            DMR_save_pointer[i] = DMR_pointer[i];
        }
    }

    ModuleBase::timer::tick("DensityMatrix", "save_DMR");
}

// calculate DMR from DMK using add_element
template <typename TK, typename TR>
void DensityMatrix<TK,TR>::cal_DMR_test()
{
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nks*(is-1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<TR>* tmp_DMR = this->_DMR[is-1];
        // set zero since this function is called in every scf step
        tmp_DMR->set_zero();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int iap = 0; iap < tmp_DMR->size_atom_pairs(); ++iap)
        {
            hamilt::AtomPair<double>& tmp_ap = tmp_DMR->get_atom_pair(iap);
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
                hamilt::BaseMatrix<TR>* tmp_matrix = tmp_ap.find_matrix(r_index[0], r_index[1], r_index[2]);
#ifdef __DEBUG
            if (tmp_matrix == nullptr)
            {
                std::cout << "tmp_matrix is nullptr" << std::endl;
                continue;
            }
#endif
                std::complex<TR> tmp_res;
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
                    for (int i = 0; i < this->_paraV->get_row_size(iat1); ++i)
                    {
                        for (int j = 0; j < this->_paraV->get_col_size(iat2); ++j)
                        {
                            // since DMK is column-major, we need to transpose it col=>row
                            tmp_res = kphase * this->_DMK[ik_begin+ik][(col_ap+j)*this->_paraV->nrow+row_ap+i];
                            tmp_matrix->add_element(i, j, tmp_res.real());
                        }
                    }
                }
            }
        }
    }
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR()
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR");
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    int ld_hk = this->_paraV->nrow;
    int ld_hk2 = 2 * ld_hk;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nks * (is - 1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        tmp_DMR->set_zero();
#ifdef _OPENMP
        #pragma omp parallel for
#endif
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
                if(GlobalV::NSPIN !=4 )
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
                    // jump DMK to fill DMR
                    // DMR is row-major, DMK is column-major
                    tmp_DMK_pointer += col_ap * this->_paraV->nrow + row_ap;
                    for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
                    {
                        DMK_real_pointer = (double*)tmp_DMK_pointer;
                        DMK_imag_pointer = DMK_real_pointer + 1;
                        BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                            kphase.real(),
                                            DMK_real_pointer,
                                            ld_hk2,
                                            tmp_DMR_pointer,
                                            1);
                        // "-" since i^2 = -1
                        BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                            -kphase.imag(),
                                            DMK_imag_pointer,
                                            ld_hk2,
                                            tmp_DMR_pointer,
                                            1);
                        tmp_DMK_pointer += 1;
                        tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
                    }
                }
                // treat DMR as pauli matrix when NSPIN=4
                if(GlobalV::NSPIN==4)
                {
                    std::vector<std::complex<double>> tmp_DMR(this->_paraV->get_col_size() * this->_paraV->get_row_size(), std::complex<double>(0.0, 0.0));
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
                        std::complex<double>* tmp_DMR_pointer = tmp_DMR.data();
                        std::complex<double>* tmp_DMK_pointer = this->_DMK[ik + ik_begin].data();
                        double* DMK_real_pointer = nullptr;
                        double* DMK_imag_pointer = nullptr;
                        // jump DMK to fill DMR
                        // DMR is row-major, DMK is column-major
                        tmp_DMK_pointer += col_ap * this->_paraV->nrow + row_ap;
                        for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
                        {
                            BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                                kphase,
                                                tmp_DMK_pointer,
                                                ld_hk,
                                                tmp_DMR_pointer,
                                                1);
                            tmp_DMK_pointer += 1;
                            tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
                        }
                    }
                    int npol = 2;
                    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
                    std::vector<int> step_trace(npol * npol, 0);
                    for (int is = 0; is < npol; is++)
                    {
                        for (int is2 = 0; is2 < npol; is2++)
                        {
                            step_trace[is * npol + is2] = this->_paraV->get_col_size(iat2) * is + is2;
                            //step_trace[is + is2 * npol] = this->_paraV->get_col_size(iat2) * is + is2;
                        }
                    }
                    std::complex<double> tmp[4];
                    double* target_DMR = tmp_matrix->get_pointer();
                    std::complex<double>* tmp_DMR_pointer = tmp_DMR.data();
                    for(int irow=0;irow<this->_paraV->get_row_size(iat1);irow += 2)
                    {
                        for(int icol=0;icol<this->_paraV->get_col_size(iat2);icol += 2)
                        {
                            // catch the 4 spin component value of one orbital pair
                            tmp[0] = tmp_DMR_pointer[icol + step_trace[0]];
                            tmp[1] = tmp_DMR_pointer[icol + step_trace[1]];
                            tmp[2] = tmp_DMR_pointer[icol + step_trace[2]];
                            tmp[3] = tmp_DMR_pointer[icol + step_trace[3]];
                            // transfer to Pauli matrix and save the real part
                            // save them back to the tmp_matrix
                            target_DMR[icol + step_trace[0]] = tmp[0].real() + tmp[3].real();
                            target_DMR[icol + step_trace[1]] = tmp[1].real() + tmp[2].real();
                            target_DMR[icol + step_trace[2]] = - tmp[1].imag() + tmp[2].imag();// (i * (rho_updown - rho_downup)).real()
                            target_DMR[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();
                        }
                        tmp_DMR_pointer += this->_paraV->get_col_size(iat2) * 2;
                        target_DMR += this->_paraV->get_col_size(iat2) * 2;
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
}


// calculate DMR from DMK using blas for gamma-only calculation
template <>
void DensityMatrix<double, double>::cal_DMR()
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR");
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nks * (is - 1); // jump this->_nks for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        tmp_DMR->set_zero();
        
#ifdef __DEBUG
        //assert(tmp_DMR->is_gamma_only() == true);
        assert(this->_nks == 1);
#endif
#ifdef _OPENMP
        #pragma omp parallel for
#endif
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
            // transpose DMK col=>row
            tmp_DMK_pointer += col_ap * this->_paraV->nrow + row_ap;
            for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
            {
                BlasConnector::axpy(this->_paraV->get_col_size(iat2), kphase, tmp_DMK_pointer, ld_hk, tmp_DMR_pointer, 1);
                tmp_DMK_pointer += 1;
                tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
}

// merge density matrix DMR with different spin
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::sum_DMR_spin()
{
    ModuleBase::TITLE("DensityMatrix", "sum_DMR_spin");
    if (this->_nspin == 1)
    {
        return;
    }
    ModuleBase::timer::tick("DensityMatrix", "sum_DMR_spin");
    if (this->_nspin == 2)
    {
        hamilt::HContainer<double>* tmp_DMR_up = this->_DMR[0];
        hamilt::HContainer<double>* tmp_DMR_down = this->_DMR[1];
        for (int i = 0; i < tmp_DMR_up->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR>& tmp_ap_up = tmp_DMR_up->get_atom_pair(i);
            hamilt::AtomPair<TR>& tmp_ap_down = tmp_DMR_down->get_atom_pair(i);
            for (int ir = 0; ir < tmp_ap_up.get_R_size(); ++ir)
            {
                const int* r_index = tmp_ap_up.get_R_index(ir);
                hamilt::BaseMatrix<double>* tmp_matrix_up = tmp_ap_up.find_matrix(r_index[0], r_index[1], r_index[2]);
                hamilt::BaseMatrix<double>* tmp_matrix_down = tmp_ap_down.find_matrix(r_index[0], r_index[1], r_index[2]);
                TR* ptr_up = tmp_matrix_up->get_pointer();
                TR* ptr_down = tmp_matrix_down->get_pointer();
                for (int i = 0; i < tmp_ap_up.get_size(); ++i)
                {
                    ptr_up[i] += ptr_down[i];
                }
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "sum_DMR_spin");
}

// read *.dmk into density matrix dm(k)
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::read_DMK(const std::string directory, const int ispin, const int ik)
{
    ModuleBase::TITLE("DensityMatrix", "read_DMK");
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
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
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
    ModuleBase::TITLE("DensityMatrix", "write_DMK");
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
//template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace elecstate

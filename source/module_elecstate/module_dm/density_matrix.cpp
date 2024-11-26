#include "density_matrix.h"

#include "module_parameter/parameter.h"
#include "module_base/libm/libm.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_cell/klist.h"

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
    delete[] this->dmr_tmp_;
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in, const int nspin, const std::vector<ModuleBase::Vector3<double>>& kvec_d, const int nk)
    : _paraV(paraV_in), _nspin(nspin), _kvec_d(kvec_d), _nk((nk > 0 && nk <= _kvec_d.size()) ? nk : _kvec_d.size())
{
    ModuleBase::TITLE("DensityMatrix", "DensityMatrix-MK");
    const int nks = _nk * _nspin;
    this->_DMK.resize(nks);
    for (int ik = 0; ik < nks; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}

template <typename TK, typename TR>
DensityMatrix<TK, TR>::DensityMatrix(const Parallel_Orbitals* paraV_in, const int nspin) :_paraV(paraV_in), _nspin(nspin), _kvec_d({ ModuleBase::Vector3<double>(0,0,0) }), _nk(1)
{
    ModuleBase::TITLE("DensityMatrix", "DensityMatrix-GO");
    this->_DMK.resize(_nspin);
    for (int ik = 0; ik < this->_nspin; ik++)
    {
        this->_DMK[ik].resize(this->_paraV->get_row_size() * this->_paraV->get_col_size());
    }
    ModuleBase::Memory::record("DensityMatrix::DMK", this->_DMK.size() * this->_DMK[0].size() * sizeof(TK));
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR(const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR");

    // To check whether DMR has been initialized
#ifdef __DEBUG
    assert(!this->_DMR.empty() && "DMR has not been initialized!");
#endif

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    int ld_hk = this->_paraV->nrow;
    int ld_hk2 = 2 * ld_hk;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
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
            std::vector<std::complex<double>> tmp_DMR;
            if (PARAM.inp.nspin == 4)
            {
                tmp_DMR.resize(tmp_ap.get_size());
            }
            for (int ir = 0; ir < tmp_ap.get_R_size(); ++ir)
            {
                const ModuleBase::Vector3<int> r_index = tmp_ap.get_R_index(ir);
                hamilt::BaseMatrix<double>* tmp_matrix = tmp_ap.find_matrix(r_index);
#ifdef __DEBUG
                if (tmp_matrix == nullptr)
                {
                    std::cout << "tmp_matrix is nullptr" << std::endl;
                    continue;
                }
#endif
                // loop over k-points
                if (PARAM.inp.nspin != 4)
                {
                    for (int ik = 0; ik < this->_nk; ++ik)
                    {
                        if(ik_in >= 0 && ik_in != ik) { continue;
}
                        // cal k_phase
                        // if TK==std::complex<double>, kphase is e^{ikR}
                        const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                        const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
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
                }

                // treat DMR as pauli matrix when NSPIN=4
                if (PARAM.inp.nspin == 4)
                {
                    tmp_DMR.assign(tmp_ap.get_size(), std::complex<double>(0.0, 0.0));
                    for (int ik = 0; ik < this->_nk; ++ik)
                    {
                        if(ik_in >= 0 && ik_in != ik) { continue;
}
                        // cal k_phase
                        // if TK==std::complex<double>, kphase is e^{ikR}
                        const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                        const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
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
                        for (int mu = 0; mu < tmp_ap.get_row_size(); ++mu)
                        {
                            BlasConnector::axpy(tmp_ap.get_col_size(),
                                                kphase,
                                                tmp_DMK_pointer,
                                                ld_hk,
                                                tmp_DMR_pointer,
                                                1);
                            tmp_DMK_pointer += 1;
                            tmp_DMR_pointer += tmp_ap.get_col_size();
                        }
                    }
                    int npol = 2;
                    // step_trace = 0 for NSPIN=1,2; ={0, 1, local_col, local_col+1} for NSPIN=4
                    int step_trace[4];
                    for (int is = 0; is < npol; is++)
                    {
                        for (int is2 = 0; is2 < npol; is2++)
                        {
                            step_trace[is * npol + is2] = tmp_ap.get_col_size() * is + is2;
                        }
                    }
                    std::complex<double> tmp[4];
                    double* target_DMR = tmp_matrix->get_pointer();
                    std::complex<double>* tmp_DMR_pointer = tmp_DMR.data();
                    for (int irow = 0; irow < tmp_ap.get_row_size(); irow += 2)
                    {
                        for (int icol = 0; icol < tmp_ap.get_col_size(); icol += 2)
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
                            target_DMR[icol + step_trace[2]]
                                = -tmp[1].imag() + tmp[2].imag(); // (i * (rho_updown - rho_downup)).real()
                            target_DMR[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();
                        }
                        tmp_DMR_pointer += tmp_ap.get_col_size() * 2;
                        target_DMR += tmp_ap.get_col_size() * 2;
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<double, double>::cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out)const{}
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out)const
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR_full");

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR_full");
    int ld_hk = this->_paraV->nrow;
    int ld_hk2 = 2 * ld_hk;
    hamilt::HContainer<std::complex<double>>* tmp_DMR = dmR_out;
    // set zero since this function is called in every scf step
    tmp_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < tmp_DMR->size_atom_pairs(); ++i)
    {
        auto& tmp_ap = tmp_DMR->get_atom_pair(i);
        int iat1 = tmp_ap.get_atom_i();
        int iat2 = tmp_ap.get_atom_j();
        // get global indexes of whole matrix for each atom in this process
        int row_ap = this->_paraV->atom_begin_row[iat1];
        int col_ap = this->_paraV->atom_begin_col[iat2];
        for (int ir = 0; ir < tmp_ap.get_R_size(); ++ir)
        {
            const ModuleBase::Vector3<int> r_index = tmp_ap.get_R_index(ir);
            auto* tmp_matrix = tmp_ap.find_matrix(r_index);
#ifdef __DEBUG
            if (tmp_matrix == nullptr)
            {
                std::cout << "tmp_matrix is nullptr" << std::endl;
                continue;
            }
#endif
            // loop over k-points
            // calculate full matrix for complex density matrix
            for (int ik = 0; ik < this->_nk; ++ik)
            {
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                std::complex<double> kphase = std::complex<double>(cosp, sinp);
                // set DMR element
                std::complex<double>* tmp_DMR_pointer = tmp_matrix->get_pointer();
                const std::complex<double>* tmp_DMK_pointer = this->_DMK[ik].data();
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
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR_full");
}

// calculate DMR from DMK using blas for gamma-only calculation
template <>
void DensityMatrix<double, double>::cal_DMR(const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR");

    assert(ik_in == -1 || ik_in == 0);

    // To check whether DMR has been initialized
#ifdef __DEBUG
    assert(!this->_DMR.empty() && "DMR has not been initialized!");
#endif

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<double>* tmp_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        tmp_DMR->set_zero();

#ifdef __DEBUG
        // assert(tmp_DMR->is_gamma_only() == true);
        assert(this->_nk == 1);
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
            const ModuleBase::Vector3<int> r_index = tmp_ap.get_R_index(0);
#ifdef __DEBUG
            assert(tmp_ap.get_R_size() == 1);
            assert(r_index.x == 0 && r_index.y == 0 && r_index.z == 0);
#endif
            hamilt::BaseMatrix<double>* tmp_matrix = tmp_ap.find_matrix(r_index);
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
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
}

// switch_dmr
template <typename TK, typename TR>
void DensityMatrix<TK, TR>::switch_dmr(const int mode)
{
    ModuleBase::TITLE("DensityMatrix", "switch_dmr");
    if (this->_nspin != 2)
    {
        return;
    }
    else
    {
        ModuleBase::timer::tick("DensityMatrix", "switch_dmr");
        switch(mode)
        {
        case 0:
            // switch to original density matrix
            if (this->dmr_tmp_ != nullptr && this->dmr_origin_.size() != 0) 
            {
                this->_DMR[0]->allocate(this->dmr_origin_.data(), false);
                delete[] this->dmr_tmp_;
                this->dmr_tmp_ = nullptr;
            }
            // else: do nothing
            break;
        case 1:
            // switch to total magnetization density matrix, dmr_up + dmr_down
            if(this->dmr_tmp_ == nullptr)
            {
                const size_t size = this->_DMR[0]->get_nnr();
                this->dmr_tmp_ = new TR[size];
                this->dmr_origin_.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin_[i] = this->_DMR[0]->get_wrapper()[i];
                    this->dmr_tmp_[i] = this->dmr_origin_[i] + this->_DMR[1]->get_wrapper()[i];
                }
                this->_DMR[0]->allocate(this->dmr_tmp_, false);
            }
            else
            {
                const size_t size = this->_DMR[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp_[i] = this->dmr_origin_[i] + this->_DMR[1]->get_wrapper()[i];
                }
            }
            break;
        case 2:
            // switch to magnetization density matrix, dmr_up - dmr_down
            if(this->dmr_tmp_ == nullptr)
            {
                const size_t size = this->_DMR[0]->get_nnr();
                this->dmr_tmp_ = new TR[size];
                this->dmr_origin_.resize(size);
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_origin_[i] = this->_DMR[0]->get_wrapper()[i];
                    this->dmr_tmp_[i] = this->dmr_origin_[i] - this->_DMR[1]->get_wrapper()[i];
                }
                this->_DMR[0]->allocate(this->dmr_tmp_, false);
            }
            else
            {
                const size_t size = this->_DMR[0]->get_nnr();
                for (int i = 0; i < size; ++i)
                {
                    this->dmr_tmp_[i] = this->dmr_origin_[i] - this->_DMR[1]->get_wrapper()[i];
                }
            }
            break;
        default:
            throw std::string("Unknown mode in switch_dmr");
        }
        ModuleBase::timer::tick("DensityMatrix", "switch_dmr");
    }
}

// T of HContainer can be double or complex<double>
template class DensityMatrix<double, double>;               // Gamma-Only case
template class DensityMatrix<std::complex<double>, double>; // Multi-k case
template class DensityMatrix<std::complex<double>, std::complex<double>>; // For EXX in future

} // namespace elecstate

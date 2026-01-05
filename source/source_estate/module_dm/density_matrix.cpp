#include "density_matrix.h"

#include "source_io/module_parameter/parameter.h"
#include "source_base/libm/libm.h"
#include "source_base/memory.h"
#include "source_base/timer.h"
#include "source_base/tool_title.h"
#include "source_cell/klist.h"

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
    ModuleBase::TITLE("DensityMatrix", "resize_DMK");
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
    ModuleBase::TITLE("DensityMatrix", "resize_gamma");
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
    using TR = double;

    // To check whether DMR has been initialized
    assert(!this->_DMR.empty() && "DMR has not been initialized!");

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    const int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        const int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<TR>*const target_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        target_DMR->set_zero();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            // get global indexes of whole matrix for each atom in this process
            const int row_ap = this->_paraV->atom_begin_row[iat1];
            const int col_ap = this->_paraV->atom_begin_col[iat2];
            const int row_size = this->_paraV->get_row_size(iat1);
            const int col_size = this->_paraV->get_col_size(iat2);
            const int mat_size = row_size * col_size;
            const int R_size = target_ap.get_R_size();
            assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");

            // calculate kphase and target_mat_ptr
            std::vector<std::vector<std::complex<double>>> kphase_vec(this->_nk, std::vector<std::complex<double>>(R_size));
            std::vector<TR*> target_DMR_mat_vec(R_size);
            for(int iR = 0; iR < R_size; ++iR)
            {
                const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(iR);
                hamilt::BaseMatrix<TR>*const target_mat = target_ap.find_matrix(R_index);
                #ifdef __DEBUG
                if (target_mat == nullptr)
                {
                    std::cout << "target_mat is nullptr" << std::endl;
                    continue;
                }
                #endif
                target_DMR_mat_vec[iR] = target_mat->get_pointer();
                for(int ik = 0; ik < this->_nk; ++ik)
                {
                    if(ik_in >= 0 && ik_in != ik) { continue; }
                    // cal k_phase
                    // if TK==std::complex<double>, kphase is e^{ikR}
                    const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
                    const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    kphase_vec[ik][iR] = std::complex<double>(cosp, sinp);
                }
            }

            std::vector<std::complex<double>> DMK_mat_trans(mat_size);
            std::vector<std::complex<double>> tmp_DMR( (PARAM.inp.nspin==4) ? mat_size*R_size : 0);
            for(int ik = 0; ik < this->_nk; ++ik)
            {
                if(ik_in >= 0 && ik_in != ik) { continue; }

                // copy column-major DMK to row-major DMK_mat_trans (for the purpose of computational efficiency)
                const std::complex<double>*const DMK_mat_ptr
                    = this->_DMK[ik + ik_begin].data()
                      + col_ap * this->_paraV->nrow + row_ap;
                for(int icol = 0; icol < col_size; ++icol) {
                    for(int irow = 0; irow < row_size; ++irow) {
                        DMK_mat_trans[irow * col_size + icol] = DMK_mat_ptr[icol * ld_hk + irow];
                }}

                // if nspin != 4, fill DMR
                // if nspin == 4, fill tmp_DMR
                for(int iR = 0; iR < R_size; ++iR)
                {
                    // (kr+i*ki) * (Dr+i*Di) = (kr*Dr-ki*Di) + i*(kr*Di+ki*Dr)
                    const std::complex<double> kphase = kphase_vec[ik][iR];
                    if(PARAM.inp.nspin != 4)                // only save real kr*Dr-ki*Di
                    {
                        TR* target_DMR_mat = target_DMR_mat_vec[iR];
                        for(int i = 0; i < mat_size; i++)
                        {
                            target_DMR_mat[i]
                                += kphase.real() * DMK_mat_trans[i].real() 
                                 - kphase.imag() * DMK_mat_trans[i].imag();
                        }
                    } else if(PARAM.inp.nspin == 4)
                    {
                        BlasConnector::axpy(mat_size,
                                            kphase,
                                            DMK_mat_trans.data(),
                                            1,
                                            &tmp_DMR[iR * mat_size],
                                            1);
                    }
                }
            }

            // if nspin == 4
            // copy tmp_DMR to fill target_DMR
            if(PARAM.inp.nspin == 4)
            {
                // step_trace ={0, 1, local_col, local_col+1} for NSPIN=4
                int step_trace[4]{};
                constexpr int npol = 2;
                for (int is = 0; is < npol; is++) {
                    for (int is2 = 0; is2 < npol; is2++) {
                        step_trace[is * npol + is2] = target_ap.get_col_size() * is + is2;
                }}

                std::complex<double> tmp[4]{};
                for(int iR = 0; iR < R_size; ++iR)
                {
                    const std::complex<double>* tmp_DMR_mat = &tmp_DMR[iR * mat_size];
                    TR* target_DMR_mat = target_DMR_mat_vec[iR];
                    for (int irow = 0; irow < row_size; irow += 2)
                    {
                        for (int icol = 0; icol < col_size; icol += 2)
                        {
                            // catch the 4 spin component value of one orbital pair
                            tmp[0] = tmp_DMR_mat[icol + step_trace[0]];
                            tmp[1] = tmp_DMR_mat[icol + step_trace[1]];
                            tmp[2] = tmp_DMR_mat[icol + step_trace[2]];
                            tmp[3] = tmp_DMR_mat[icol + step_trace[3]];
                            // transfer to Pauli matrix and save the real part
                            // save them back to the target_mat
                            target_DMR_mat[icol + step_trace[0]] = tmp[0].real() + tmp[3].real();  // (rho_upup + rho_downdown).real()
                            target_DMR_mat[icol + step_trace[1]] = tmp[1].real() + tmp[2].real();  // (rho_updown + rho_downup).real()
                            target_DMR_mat[icol + step_trace[2]] = -tmp[1].imag() + tmp[2].imag(); // (i * (rho_updown - rho_downup)).real()
                            target_DMR_mat[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();  // (rho_upup - rho_downdown).real()
                        }
                        tmp_DMR_mat += col_size * 2;
                        target_DMR_mat += col_size * 2;
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR_td(const UnitCell& ucell, const ModuleBase::Vector3<double> At, const int ik_in)
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR_td");
    using TR = double;
    // To check whether DMR has been initialized
    assert(!this->_DMR.empty() && "DMR has not been initialized!");

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR_td");
    const int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        const int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<TR>*const target_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        target_DMR->set_zero();
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<TR>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            // get global indexes of whole matrix for each atom in this process
            const int row_ap = this->_paraV->atom_begin_row[iat1];
            const int col_ap = this->_paraV->atom_begin_col[iat2];
            const int row_size = this->_paraV->get_row_size(iat1);
            const int col_size = this->_paraV->get_col_size(iat2);
            const int mat_size = row_size * col_size;
            const int R_size = target_ap.get_R_size();
            assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");

            // calculate kphase and target_mat_ptr
            std::vector<std::vector<std::complex<double>>> kphase_vec(this->_nk, std::vector<std::complex<double>>(R_size));
            std::vector<TR*> target_DMR_mat_vec(R_size);
            for(int iR = 0; iR < R_size; ++iR)
            {
                const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(iR);
                hamilt::BaseMatrix<TR>*const target_mat = target_ap.find_matrix(R_index);
                #ifdef __DEBUG
                if (target_mat == nullptr)
                {
                    std::cout << "target_mat is nullptr" << std::endl;
                    continue;
                }
                #endif
                target_DMR_mat_vec[iR] = target_mat->get_pointer();
                double arg_td = 0.0;
                //cal tddft phase for hybrid gauge
                ModuleBase::Vector3<double> dtau = ucell.cal_dtau(iat1, iat2, R_index);
                arg_td = At * dtau * ucell.lat0;
                for(int ik = 0; ik < this->_nk; ++ik)
                {
                    if(ik_in >= 0 && ik_in != ik) { continue; }
                    // cal k_phase
                    // if TK==std::complex<double>, kphase is e^{ikR}
                    const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
                    const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI + arg_td;
                    double sinp, cosp;
                    ModuleBase::libm::sincos(arg, &sinp, &cosp);
                    kphase_vec[ik][iR] = std::complex<double>(cosp, sinp);
                }
            }

            std::vector<std::complex<double>> DMK_mat_trans(mat_size);
            std::vector<std::complex<double>> tmp_DMR( (PARAM.inp.nspin==4) ? mat_size*R_size : 0);
            for(int ik = 0; ik < this->_nk; ++ik)
            {
                if(ik_in >= 0 && ik_in != ik) { continue; }

                // copy column-major DMK to row-major DMK_mat_trans (for the purpose of computational efficiency)
                const std::complex<double>*const DMK_mat_ptr
                    = this->_DMK[ik + ik_begin].data()
                      + col_ap * this->_paraV->nrow + row_ap;
                for(int icol = 0; icol < col_size; ++icol) {
                    for(int irow = 0; irow < row_size; ++irow) {
                        DMK_mat_trans[irow * col_size + icol] = DMK_mat_ptr[icol * ld_hk + irow];
                }}

                // if nspin != 4, fill DMR
                // if nspin == 4, fill tmp_DMR
                for(int iR = 0; iR < R_size; ++iR)
                {
                    // (kr+i*ki) * (Dr+i*Di) = (kr*Dr-ki*Di) + i*(kr*Di+ki*Dr)
                    const std::complex<double> kphase = kphase_vec[ik][iR];
                    if(PARAM.inp.nspin != 4)                // only save real kr*Dr-ki*Di
                    {
                        TR* target_DMR_mat = target_DMR_mat_vec[iR];
                        for(int i = 0; i < mat_size; i++)
                        {
                            target_DMR_mat[i]
                                += kphase.real() * DMK_mat_trans[i].real() 
                                 - kphase.imag() * DMK_mat_trans[i].imag();
                        }
                    } else if(PARAM.inp.nspin == 4)
                    {
                        BlasConnector::axpy(mat_size,
                                            kphase,
                                            DMK_mat_trans.data(),
                                            1,
                                            &tmp_DMR[iR * mat_size],
                                            1);
                    }
                }
            }

            // if nspin == 4
            // copy tmp_DMR to fill target_DMR
            if(PARAM.inp.nspin == 4)
            {
                // step_trace ={0, 1, local_col, local_col+1} for NSPIN=4
                int step_trace[4]{};
                constexpr int npol = 2;
                for (int is = 0; is < npol; is++) {
                    for (int is2 = 0; is2 < npol; is2++) {
                        step_trace[is * npol + is2] = target_ap.get_col_size() * is + is2;
                }}

                std::complex<double> tmp[4]{};
                for(int iR = 0; iR < R_size; ++iR)
                {
                    const std::complex<double>* tmp_DMR_mat = &tmp_DMR[iR * mat_size];
                    TR* target_DMR_mat = target_DMR_mat_vec[iR];
                    for (int irow = 0; irow < row_size; irow += 2)
                    {
                        for (int icol = 0; icol < col_size; icol += 2)
                        {
                            // catch the 4 spin component value of one orbital pair
                            tmp[0] = tmp_DMR_mat[icol + step_trace[0]];
                            tmp[1] = tmp_DMR_mat[icol + step_trace[1]];
                            tmp[2] = tmp_DMR_mat[icol + step_trace[2]];
                            tmp[3] = tmp_DMR_mat[icol + step_trace[3]];
                            // transfer to Pauli matrix and save the real part
                            // save them back to the target_mat
                            target_DMR_mat[icol + step_trace[0]] = tmp[0].real() + tmp[3].real();  // (rho_upup + rho_downdown).real()
                            target_DMR_mat[icol + step_trace[1]] = tmp[1].real() + tmp[2].real();  // (rho_updown + rho_downup).real()
                            target_DMR_mat[icol + step_trace[2]] = -tmp[1].imag() + tmp[2].imag(); // (i * (rho_updown - rho_downup)).real()
                            target_DMR_mat[icol + step_trace[3]] = tmp[0].real() - tmp[3].real();  // (rho_upup - rho_downdown).real()
                        }
                        tmp_DMR_mat += col_size * 2;
                        target_DMR_mat += col_size * 2;
                    }
                }
            }
        }
    }
    ModuleBase::timer::tick("DensityMatrix", "cal_DMR_td");
}

// calculate DMR from DMK using blas for multi-k calculation
template <>
void DensityMatrix<double, double>::cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out)const{}
template <>
void DensityMatrix<std::complex<double>, double>::cal_DMR_full(hamilt::HContainer<std::complex<double>>* dmR_out)const
{
    ModuleBase::TITLE("DensityMatrix", "cal_DMR_full");

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR_full");
    const int ld_hk = this->_paraV->nrow;
    hamilt::HContainer<std::complex<double>>* target_DMR = dmR_out;
    // set zero since this function is called in every scf step
    target_DMR->set_zero();
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
    {
        hamilt::AtomPair<std::complex<double>>& target_ap = target_DMR->get_atom_pair(i);
        const int iat1 = target_ap.get_atom_i();
        const int iat2 = target_ap.get_atom_j();
        // get global indexes of whole matrix for each atom in this process
        const int row_ap = this->_paraV->atom_begin_row[iat1];
        const int col_ap = this->_paraV->atom_begin_col[iat2];
        const int row_size = this->_paraV->get_row_size(iat1);
        const int col_size = this->_paraV->get_col_size(iat2);
        const int mat_size = row_size * col_size;
        const int R_size = target_ap.get_R_size();

        // calculate kphase and target_mat_ptr
        std::vector<std::vector<std::complex<double>>> kphase_vec(this->_nk, std::vector<std::complex<double>>(R_size));
        std::vector<std::complex<double>*> target_DMR_mat_vec(R_size);
        for(int iR = 0; iR < R_size; ++iR)
        {
            const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(iR);
            hamilt::BaseMatrix<std::complex<double>>*const target_mat = target_ap.find_matrix(R_index);
            #ifdef __DEBUG
            if (target_mat == nullptr)
            {
                std::cout << "target_mat is nullptr" << std::endl;
                continue;
            }
            #endif
            target_DMR_mat_vec[iR] = target_mat->get_pointer();
            for(int ik = 0; ik < this->_nk; ++ik)
            {
                // cal k_phase
                // if TK==std::complex<double>, kphase is e^{ikR}
                const ModuleBase::Vector3<double> dR(R_index[0], R_index[1], R_index[2]);
                const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                double sinp, cosp;
                ModuleBase::libm::sincos(arg, &sinp, &cosp);
                kphase_vec[ik][iR] = std::complex<double>(cosp, sinp);
            }
        }

        std::vector<std::complex<double>> DMK_mat_trans(mat_size);
        for(int ik = 0; ik < this->_nk; ++ik)
        {
            // copy column-major DMK to row-major DMK_mat_trans (for the purpose of computational efficiency)
            const std::complex<double>*const DMK_mat_ptr
                = this->_DMK[ik].data()
                  + col_ap * this->_paraV->nrow + row_ap;
            for(int icol = 0; icol < col_size; ++icol) {
                for(int irow = 0; irow < row_size; ++irow) {
                    DMK_mat_trans[irow * col_size + icol] = DMK_mat_ptr[icol * ld_hk + irow];
            }}

            for(int iR = 0; iR < R_size; ++iR)
            {
                const std::complex<double> kphase = kphase_vec[ik][iR];
                std::complex<double>* target_DMR_mat = target_DMR_mat_vec[iR];
                BlasConnector::axpy(mat_size,
                                    kphase,
                                    DMK_mat_trans.data(),
                                    1,
                                    target_DMR_mat,
                                    1);
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
    assert(this->_nk == 1);

    // To check whether DMR has been initialized
    assert(!this->_DMR.empty() && "DMR has not been initialized!");

    ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    const int ld_hk = this->_paraV->nrow;
    for (int is = 1; is <= this->_nspin; ++is)
    {
        const int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
        hamilt::HContainer<double>*const target_DMR = this->_DMR[is - 1];
        // set zero since this function is called in every scf step
        target_DMR->set_zero();

        // assert(target_DMR->is_gamma_only() == true);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int i = 0; i < target_DMR->size_atom_pairs(); ++i)
        {
            hamilt::AtomPair<double>& target_ap = target_DMR->get_atom_pair(i);
            const int iat1 = target_ap.get_atom_i();
            const int iat2 = target_ap.get_atom_j();
            // get global indexes of whole matrix for each atom in this process
            const int row_ap = this->_paraV->atom_begin_row[iat1];
            const int col_ap = this->_paraV->atom_begin_col[iat2];
            assert(row_ap != -1 && col_ap != -1 && "Atom-pair not belong this process");
            const int row_size = this->_paraV->get_row_size(iat1);
            const int col_size = this->_paraV->get_col_size(iat2);
            const int R_size = target_ap.get_R_size();
            assert(R_size == 1);
            const ModuleBase::Vector3<int> R_index = target_ap.get_R_index(0);
            assert(R_index.x == 0 && R_index.y == 0 && R_index.z == 0);
            hamilt::BaseMatrix<double>*const target_mat = target_ap.find_matrix(R_index);
            #ifdef __DEBUG
            if (target_mat == nullptr)
            {
                std::cout << "target_mat is nullptr" << std::endl;
                continue;
            }
            #endif
            // k index
            constexpr double kphase = 1;
            // transpose DMK col=>row
            const double* DMK_ptr
                = this->_DMK[0 + ik_begin].data()
                  + col_ap * this->_paraV->nrow + row_ap;
            // set DMR element
            double* target_DMR_ptr = target_mat->get_pointer();
            for (int mu = 0; mu < row_size; ++mu)
            {
                BlasConnector::axpy(col_size,
                                    kphase,
                                    DMK_ptr,
                                    ld_hk,
                                    target_DMR_ptr,
                                    1);
                DMK_ptr += 1;
                target_DMR_ptr += col_size;
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

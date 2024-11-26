#include "module_elecstate/module_dm/density_matrix.h"
#include "module_base/timer.h"
#include "module_parameter/parameter.h"
#include "module_base/libm/libm.h"
namespace elecstate
{
    template<>
    void DensityMatrix<std::complex<double>, std::complex<double>>::cal_DMR(int ik_in)
    {
        ModuleBase::TITLE("DensityMatrix", "cal_DMR");
        ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
        for (int is = 1; is <= this->_nspin; ++is)
        {
            int ik_begin = this->_nk * (is - 1); // jump this->_nk for spin_down if nspin==2
            hamilt::HContainer<std::complex<double>>* tmp_DMR = this->_DMR[is - 1];
            // set zero since this function is called in every scf step
            tmp_DMR->set_zero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < tmp_DMR->size_atom_pairs(); ++i)
            {
                hamilt::AtomPair<std::complex<double>>& tmp_ap = tmp_DMR->get_atom_pair(i);
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
                    const ModuleBase::Vector3<int> r_index = tmp_ap.get_R_index(ir);
                    hamilt::BaseMatrix<std::complex<double>>* tmp_matrix = tmp_ap.find_matrix(r_index);
#ifdef __DEBUG
                    if (tmp_matrix == nullptr)
                    {
                        std::cout << "tmp_matrix is nullptr" << std::endl;
                        continue;
                    }
#endif
                    // loop over k-points
                    if (PARAM.inp.nspin != 4) {
                        for (int ik = 0; ik < this->_nk; ++ik)
                        {
                            if (ik_in >= 0 && ik_in != ik) continue;
                            // cal k_phase
                            // if TK==std::complex<double>, kphase is e^{ikR}
                            const ModuleBase::Vector3<double> dR(r_index[0], r_index[1], r_index[2]);
                            const double arg = (this->_kvec_d[ik] * dR) * ModuleBase::TWO_PI;
                            double sinp, cosp;
                            ModuleBase::libm::sincos(arg, &sinp, &cosp);
                            std::complex<double> kphase = std::complex<double>(cosp, sinp);
                            // set DMR element
                            std::complex<double>* tmp_DMR_pointer = tmp_matrix->get_pointer();
                            std::complex<double>* tmp_DMK_pointer = this->_DMK[ik + ik_begin].data();
                            // jump DMK to fill DMR
                            // DMR is row-major, DMK is column-major
                            tmp_DMK_pointer += col_ap * this->_paraV->nrow + row_ap;
                            for (int mu = 0; mu < this->_paraV->get_row_size(iat1); ++mu)
                            {
                                BlasConnector::axpy(this->_paraV->get_col_size(iat2),
                                    kphase,
                                    tmp_DMK_pointer,
                                    this->_paraV->get_row_size(),
                                    tmp_DMR_pointer,
                                    1);
                                tmp_DMK_pointer += 1;
                                tmp_DMR_pointer += this->_paraV->get_col_size(iat2);
                            }
                        }
}
                    // treat DMR as pauli matrix when NSPIN=4
                    if (PARAM.inp.nspin == 4) {
                        throw std::runtime_error("complex DM(R) with NSPIN=4 is not implemented yet");
}
                }
            }
        }
        ModuleBase::timer::tick("DensityMatrix", "cal_DMR");
    }
    // template class DensityMatrix<std::complex<double>, std::complex<double>>;
}
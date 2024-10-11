#include "esolver_ks_lcao_tddft.h"

#include "module_io/cal_r_overlap_R.h"
#include "module_io/dipole_io.h"
#include "module_io/rho_io.h"
#include "module_io/td_current_io.h"
#include "module_io/write_HS.h"
#include "module_io/write_HS_R.h"
#include "module_io/write_wfc_nao.h"

//--------------temporary----------------------------
#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/scalapack_connector.h"
#include "module_base/lapack_connector.h"
#include "module_elecstate/module_charge/symmetry_rho.h"
#include "module_elecstate/occupy.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_domain.h" // need divide_HS_in_frag
#include "module_hamilt_lcao/module_tddft/evolve_elec.h"
#include "module_hamilt_lcao/module_tddft/td_velocity.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_io/print_info.h"

//-----HSolver ElecState Hamilt--------
#include "module_elecstate/elecstate_lcao.h"
#include "module_elecstate/elecstate_lcao_tddft.h"
#include "module_hamilt_lcao/hamilt_lcaodft/hamilt_lcao.h"
#include "module_hsolver/hsolver_lcao.h"
#include "module_parameter/parameter.h"
#include "module_psi/psi.h"

//-----force& stress-------------------
#include "module_hamilt_lcao/hamilt_lcaodft/FORCE_STRESS.h"

//---------------------------------------------------

namespace ModuleESolver
{

// use the original formula (Hamiltonian matrix) to calculate energy density
// matrix
void ESolver_KS_LCAO_TDDFT::cal_edm_tddft()
{
    // mohan add 2024-03-27
    const int nlocal = PARAM.globalv.nlocal;
    assert(nlocal >= 0);

    dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)
        ->get_DM()
        ->EDMK.resize(kv.get_nks());
    for (int ik = 0; ik < kv.get_nks(); ++ik) {
        std::complex<double>* tmp_dmk
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->get_DMK_pointer(ik);

        ModuleBase::ComplexMatrix& tmp_edmk
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->EDMK[ik];

        const Parallel_Orbitals* tmp_pv
            = dynamic_cast<elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)->get_DM()->get_paraV_pointer();

#ifdef __MPI

        // mohan add 2024-03-27
        //! be careful, the type of nloc is 'long'
        //! whether the long type is safe, needs more discussion
        const long nloc = this->pv.nloc;
        const int ncol = this->pv.ncol;
        const int nrow = this->pv.nrow;

        tmp_edmk.create(ncol, nrow);
        complex<double>* Htmp = new complex<double>[nloc];
        complex<double>* Sinv = new complex<double>[nloc];
        complex<double>* tmp1 = new complex<double>[nloc];
        complex<double>* tmp2 = new complex<double>[nloc];
        complex<double>* tmp3 = new complex<double>[nloc];
        complex<double>* tmp4 = new complex<double>[nloc];

        ModuleBase::GlobalFunc::ZEROS(Htmp, nloc);
        ModuleBase::GlobalFunc::ZEROS(Sinv, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp1, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp2, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp3, nloc);
        ModuleBase::GlobalFunc::ZEROS(tmp4, nloc);

        const int inc = 1;

        hamilt::MatrixBlock<complex<double>> h_mat;
        hamilt::MatrixBlock<complex<double>> s_mat;

        p_hamilt->matrix(h_mat, s_mat);
        zcopy_(&nloc, h_mat.p, &inc, Htmp, &inc);
        zcopy_(&nloc, s_mat.p, &inc, Sinv, &inc);

        vector<int> ipiv(nloc, 0);
        int info = 0;
        const int one_int = 1;

        pzgetrf_(&nlocal, &nlocal, Sinv, &one_int, &one_int, this->pv.desc, ipiv.data(), &info);

        int lwork = -1;
        int liwork = -1;

        // if lwork == -1, then the size of work is (at least) of length 1.
        std::vector<std::complex<double>> work(1, 0);

        // if liwork = -1, then the size of iwork is (at least) of length 1.
        std::vector<int> iwork(1, 0);

        pzgetri_(&nlocal,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->pv.desc,
                 ipiv.data(),
                 work.data(),
                 &lwork,
                 iwork.data(),
                 &liwork,
                 &info);

        lwork = work[0].real();
        work.resize(lwork, 0);
        liwork = iwork[0];
        iwork.resize(liwork, 0);

        pzgetri_(&nlocal,
                 Sinv,
                 &one_int,
                 &one_int,
                 this->pv.desc,
                 ipiv.data(),
                 work.data(),
                 &lwork,
                 iwork.data(),
                 &liwork,
                 &info);

        const char N_char = 'N';
        const char T_char = 'T';
        const complex<double> one_float = {1.0, 0.0};
        const complex<double> zero_float = {0.0, 0.0};
        const complex<double> half_float = {0.5, 0.0};

        pzgemm_(&N_char,
                &N_char,
                &nlocal,
                &nlocal,
                &nlocal,
                &one_float,
                tmp_dmk,
                &one_int,
                &one_int,
                this->pv.desc,
                Htmp,
                &one_int,
                &one_int,
                this->pv.desc,
                &zero_float,
                tmp1,
                &one_int,
                &one_int,
                this->pv.desc);

        pzgemm_(&N_char,
                &N_char,
                &nlocal,
                &nlocal,
                &nlocal,
                &one_float,
                tmp1,
                &one_int,
                &one_int,
                this->pv.desc,
                Sinv,
                &one_int,
                &one_int,
                this->pv.desc,
                &zero_float,
                tmp2,
                &one_int,
                &one_int,
                this->pv.desc);

        pzgemm_(&N_char,
                &N_char,
                &nlocal,
                &nlocal,
                &nlocal,
                &one_float,
                Sinv,
                &one_int,
                &one_int,
                this->pv.desc,
                Htmp,
                &one_int,
                &one_int,
                this->pv.desc,
                &zero_float,
                tmp3,
                &one_int,
                &one_int,
                this->pv.desc);

        pzgemm_(&N_char,
                &N_char,
                &nlocal,
                &nlocal,
                &nlocal,
                &one_float,
                tmp3,
                &one_int,
                &one_int,
                this->pv.desc,
                tmp_dmk,
                &one_int,
                &one_int,
                this->pv.desc,
                &zero_float,
                tmp4,
                &one_int,
                &one_int,
                this->pv.desc);

        pzgeadd_(&N_char,
                 &nlocal,
                 &nlocal,
                 &half_float,
                 tmp2,
                 &one_int,
                 &one_int,
                 this->pv.desc,
                 &half_float,
                 tmp4,
                 &one_int,
                 &one_int,
                 this->pv.desc);

        zcopy_(&nloc, tmp4, &inc, tmp_edmk.c, &inc);

        delete[] Htmp;
        delete[] Sinv;
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp3;
        delete[] tmp4;
#else
        // for serial version
        tmp_edmk.create(this->pv.ncol, this->pv.nrow);
        ModuleBase::ComplexMatrix Sinv(nlocal, nlocal);
        ModuleBase::ComplexMatrix Htmp(nlocal, nlocal);

        hamilt::MatrixBlock<complex<double>> h_mat;
        hamilt::MatrixBlock<complex<double>> s_mat;

        p_hamilt->matrix(h_mat, s_mat);
        // cout<<"hmat "<<h_mat.p[0]<<endl;
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                Htmp(i, j) = h_mat.p[i * nlocal + j];
                Sinv(i, j) = s_mat.p[i * nlocal + j];
            }
        }
        int INFO = 0;

        int lwork = 3 * nlocal - 1; // tmp
        std::complex<double>* work = new std::complex<double>[lwork];
        ModuleBase::GlobalFunc::ZEROS(work, lwork);

        int IPIV[nlocal];

        LapackConnector::zgetrf(nlocal, nlocal, Sinv, nlocal, IPIV, &INFO);
        LapackConnector::zgetri(nlocal, Sinv, nlocal, IPIV, work, lwork, &INFO);
        // I just use ModuleBase::ComplexMatrix temporarily, and will change it
        // to complex<double>*
        ModuleBase::ComplexMatrix tmp_dmk_base(nlocal, nlocal);
        for (int i = 0; i < nlocal; i++)
        {
            for (int j = 0; j < nlocal; j++)
            {
                tmp_dmk_base(i, j) = tmp_dmk[i * nlocal + j];
            }
        }
        tmp_edmk = 0.5 * (Sinv * Htmp * tmp_dmk_base + tmp_dmk_base * Htmp * Sinv);
        delete[] work;
#endif
    }
    return;
}
} // namespace ModuleESolver

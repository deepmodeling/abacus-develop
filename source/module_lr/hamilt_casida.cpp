#include "hamilt_casida.h"
#include "module_lr/utils/lr_util_print.h"
namespace LR
{
    template<typename T>
    std::vector<T> HamiltLR<T>::matrix()const
    {
        ModuleBase::TITLE("HamiltLR", "matrix");
        const int no = this->nocc[0];
        const int nv = this->nvirt[0];
        const auto& px = this->pX[0];
        const int ldim = nk * px.get_local_size();
        int npairs = no * nv;
        std::vector<T> Amat_full(this->nk * npairs * this->nk * npairs, 0.0);
        for (int ik = 0;ik < this->nk;++ik) {
            for (int j = 0;j < no;++j) {
                for (int b = 0;b < nv;++b)
                {//calculate A^{ai} for each bj
                    int bj = j * nv + b;    //global
                    int kbj = ik * npairs + bj; //global
                    psi::Psi<T> X_bj(1, 1, this->nk * px.get_local_size(), this->nk * px.get_local_size(), true); // k1-first, like in iterative solver
                    X_bj.zero_out();
                    // X_bj(0, 0, lj * px.get_row_size() + lb) = this->one();
                    int lj = px.global2local_col(j);
                    int lb = px.global2local_row(b);
                    if (px.in_this_processor(b, j)) { X_bj(0, 0, ik * px.get_local_size() + lj * px.get_row_size() + lb) = this->one(); }
                    psi::Psi<T> A_aibj(1, 
                                       1, 
                                       this->nk * px.get_local_size(),
                                       this->nk * px.get_local_size(),
                                       true); // k1-first
                    A_aibj.zero_out();

                    this->cal_dm_trans(0, X_bj.get_pointer());
                    hamilt::Operator<T>* node(this->ops);
                    while (node != nullptr)
                    {   // act() on and return the k1-first type of psi
                        node->act(1, ldim, /*npol=*/1, X_bj.get_pointer(), A_aibj.get_pointer());
                        node = (hamilt::Operator<T>*)(node->next_op);
                    }
                    // reduce ai for a fixed bj
                    A_aibj.fix_kb(0, 0);
#ifdef __MPI
                    for (int ik_ai = 0;ik_ai < this->nk;++ik_ai) {
                        LR_Util::gather_2d_to_full(px, &A_aibj.get_pointer()[ik_ai * px.get_local_size()],
                            Amat_full.data() + kbj * this->nk * npairs /*col, bj*/ + ik_ai * npairs/*row, ai*/,
                            false, nv, no);
}
#endif
                }
}
}
        // output Amat
        std::cout << "Full A matrix: (elements < 1e-10 is set to 0)" << std::endl;
        LR_Util::print_value(Amat_full.data(), nk * npairs, nk * npairs);
        return Amat_full;
    }

    template<> double HamiltLR<double>::one() const { return 1.0; }
    template<> std::complex<double> HamiltLR<std::complex<double>>::one() const { return std::complex<double>(1.0, 0.0); }

    template class HamiltLR<double>;
    template class HamiltLR<std::complex<double>>;
}
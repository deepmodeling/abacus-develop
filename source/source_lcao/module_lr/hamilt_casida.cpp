#include "hamilt_casida.h"
#include <fstream>
#include <cstdlib>
#include "source_lcao/module_lr/utils/lr_util_print.h"
namespace LR
{
    template<typename T>
    std::vector<T> HamiltLR<T>::matrix()const
    {
        ModuleBase::TITLE("HamiltLR", "matrix");
        const int no = this->nocc[0];
        const int nv = this->nvirt[0];
        const auto& px = this->pX[0];
        const int ldim = nk * px.size();
        int npairs = no * nv;
        std::vector<T> Amat_full(this->nk * npairs * this->nk * npairs, 0.0);
        for (int ik = 0;ik < this->nk;++ik) {
            for (int j = 0;j < no;++j) {
                for (int b = 0;b < nv;++b)
                {//calculate A^{ai} for each bj
                    int bj = j * nv + b;    //global
                    int kbj = ik * npairs + bj; //global
                    // Keep valid backing addresses even when this rank owns no pairs.
                    // ldim remains zero: all ranks still execute the operator collectives.
                    const int storage_size = std::max(ldim, 1);
                    std::vector<T> X_bj(storage_size, T(0));
                    std::vector<T> A_aibj(storage_size, T(0));
                    const int lj = px.icloc(j);
                    const int lb = px.irloc(b);
                    if (px.in_this_processor(b, j))
                    {
                        X_bj[ik * px.size() + lj * px.nrloc() + lb] = this->one();
                    }

                    this->cal_dm_trans(0, X_bj.data());
                    hamilt::Operator<T>* node(this->ops);
                    while (node != nullptr)
                    {   // act() on and return the k1-first type of psi
                        node->act(1, ldim, /*npol=*/1, X_bj.data(), A_aibj.data());
                        node = (hamilt::Operator<T>*)(node->next_op);
                    }
                    // reduce ai for a fixed bj
#ifdef __MPI
                    for (int ik_ai = 0;ik_ai < this->nk;++ik_ai) {
                        LR_Util::gather_2d_to_full(px, &A_aibj.data()[ik_ai * px.size()],
                            Amat_full.data() + kbj * this->nk * npairs /*col, bj*/ + ik_ai * npairs/*row, ai*/,
                            false, nv, no);
                    }
#else
                    std::memcpy(Amat_full.data() + kbj * this->nk * npairs, A_aibj.data(), this->nk * npairs * sizeof(T));
#endif
                }
}
}
        return Amat_full;
    }

    template<> double HamiltLR<double>::one() const { return 1.0; }
    template<> std::complex<double> HamiltLR<std::complex<double>>::one() const { return std::complex<double>(1.0, 0.0); }

    template class HamiltLR<double>;
    template class HamiltLR<std::complex<double>>;
}
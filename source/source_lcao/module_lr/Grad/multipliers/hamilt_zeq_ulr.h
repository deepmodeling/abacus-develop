#pragma once
#include "source_hamilt/hamilt.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include <cassert>
#include <cstring>
#include <vector>

namespace LR
{
    /// @brief Common skeleton of the open-shell (spin-unrestricted) Z-vector operators.
    ///
    /// Both sides of the Z-vector equation have the same block structure as the open-shell
    /// Casida Hamiltonian `HamiltULR`: the vector is the concatenation [up-block | down-block]
    /// and the operator is a 2x2 array of spin blocks
    ///     ops[(s_out << 1) + s_in],
    /// because every kernel action carries a free outer spin and a summed inner spin,
    ///     $K_{pq\sigma}[D]=\sum_{\kappa\lambda\sigma'}K_{pq\sigma,\kappa\lambda\sigma'}D_{\kappa\lambda\sigma'}$.
    /// The exchange part is $\propto\delta_{\sigma\sigma'}$, so EXX operators only ever go on
    /// the diagonal blocks 0 and 3.
    ///
    /// Derived classes fill `ops` in their constructor and implement `set_dm`, which rebuilds
    /// whatever density matrices the operators read from the `s_in` block of X.
    template<typename T>
    class ZeqULR
    {
    public:
        ZeqULR(const std::vector<int>& nocc_in,
            const std::vector<int>& nvirt_in,
            const std::vector<Parallel_2D>& pX_in,
            const int nk_in)
            : nocc(nocc_in), nvirt(nvirt_in), pX(pX_in), nk(nk_in),
            ldim(nk_in* (pX_in[0].get_local_size() + pX_in[1].get_local_size())),
            gdim(nk_in* (nocc_in[0] * nvirt_in[0] + nocc_in[1] * nvirt_in[1])),
            ops(4, nullptr)
        {
        }
        virtual ~ZeqULR() { for (auto& op : this->ops) { delete op; } }

        void hPsi(const T* const psi_in, T* const hpsi, const int ld_psi, const int nband) const
        {
            assert(ld_psi == this->ldim);
            const std::vector<int> ldim_is = { nk * pX[0].get_local_size(), nk * pX[1].get_local_size() };
            for (int ib = 0;ib < nband;++ib)
            {
                const int offset_band = ib * ld_psi;
                // Terms that are not bilinear in a (out, in) spin pair -- currently only the
                // $g^{xc}$ term of the right-hand side, which is quadratic in $D^X$ and needs
                // both transition-density channels on the grid at once.
                this->band_extra(psi_in + offset_band, hpsi + offset_band);
                for (int is_in : {0, 1})
                {
                    const int offset_in = offset_band + is_in * ldim_is[0];
                    this->set_dm(is_in, psi_in + offset_in);
                    for (int is_out : {0, 1})
                    {
                        const int offset_out = offset_band + is_out * ldim_is[0];
                        hamilt::Operator<T>* node(this->ops[(is_out << 1) + is_in]);
                        while (node != nullptr)
                        {
                            node->act(/*nband=*/1, ldim_is[is_in], /*npol=*/1,
                                psi_in + offset_in, hpsi + offset_out);
                            node = (hamilt::Operator<T>*)(node->next_op);
                        }
                    }
                }
            }
        }

        /// @brief The full (replicated) matrix, column by column. Only used by the LAPACK solver.
        std::vector<T> matrix() const
        {
            ModuleBase::TITLE("ZeqULR", "matrix");
            const std::vector<int> npairs = { nocc[0] * nvirt[0], nocc[1] * nvirt[1] };
            const std::vector<int> ldim_is = { nk * pX[0].get_local_size(), nk * pX[1].get_local_size() };
            const std::vector<int> gdim_is = { nk * npairs[0], nk * npairs[1] };
            std::vector<T> mat_full(static_cast<std::size_t>(gdim) * gdim, T(0));
            for (int is_in : {0, 1})
            {
                const auto& px = this->pX[is_in];
                const int loffset_in = is_in * ldim_is[0];
                const int goffset_in = is_in * gdim_is[0];
                for (int ik_in = 0;ik_in < nk;++ik_in)
                {
                    for (int j = 0;j < nocc[is_in];++j)
                    {
                        for (int b = 0;b < nvirt[is_in];++b)
                        {
                            const int gcol = goffset_in + ik_in * npairs[is_in] + j * nvirt[is_in] + b;
                            std::vector<T> X_col(this->ldim, T(0));
                            const int lj = px.global2local_col(j);
                            const int lb = px.global2local_row(b);
                            if (px.in_this_processor(b, j))
                            {
                                X_col[loffset_in + ik_in * px.get_local_size() + lj * px.get_row_size() + lb] = T(1);
                            }
                            this->set_dm(is_in, X_col.data() + loffset_in);
                            std::vector<T> col(this->ldim, T(0));
                            for (int is_out : {0, 1})
                            {
                                const int loffset_out = is_out * ldim_is[0];
                                const int goffset_out = is_out * gdim_is[0];
                                const auto& pax = this->pX[is_out];
                                hamilt::Operator<T>* node(this->ops[(is_out << 1) + is_in]);
                                while (node != nullptr)
                                {
                                    node->act(1, ldim_is[is_in], /*npol=*/1,
                                        X_col.data() + loffset_in, col.data() + loffset_out);
                                    node = (hamilt::Operator<T>*)(node->next_op);
                                }
#ifdef __MPI
                                for (int ik_out = 0;ik_out < this->nk;++ik_out)
                                {
                                    LR_Util::gather_2d_to_full(pax,
                                        col.data() + loffset_out + ik_out * pax.get_local_size(),
                                        mat_full.data() + static_cast<std::size_t>(gcol) * gdim
                                        + goffset_out + ik_out * npairs[is_out],
                                        false, nvirt[is_out], nocc[is_out]);
                                }
#else
                                std::memcpy(mat_full.data() + static_cast<std::size_t>(gcol) * gdim + goffset_out,
                                    col.data() + loffset_out, gdim_is[is_out] * sizeof(T));
#endif
                            }
                        }
                    }
                }
            }
            return mat_full;
        }

        const std::vector<int>& nocc;
        const std::vector<int>& nvirt;
        const std::vector<Parallel_2D>& pX;
        const int nk = 1;
        const int ldim = 1;
        const int gdim = 1;

    protected:
        /// @brief rebuild the density matrices the operators read, from the `is` block of X
        virtual void set_dm(const int is, const T* const X) const = 0;
        /// @brief optional per-band contribution that does not fit the 2x2 block structure.
        /// NOTE `matrix()` deliberately does NOT call this: it exists for the right-hand side,
        /// which is not a linear operator, while `matrix()` is only ever used for the LHS.
        virtual void band_extra(const T* const X, T* const out) const {}
        std::vector<hamilt::Operator<T>*> ops;
    };
}

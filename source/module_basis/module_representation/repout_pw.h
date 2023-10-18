#ifndef REPOUT_PW_H
#define REPOUT_PW_H

#include "module_basis/module_representation/repout.h"
#include "module_base/macros.h"
/*
    RepOut_PW: output representation of plane wave basis
    (note that in class Representation, the intermediate representation is pw, so 
    this class is the simplest case that only need to directly output, copy psig to psi_out,
    but needed to multiply by psi_in.

    For psi initialization, the psi_in would be identity matrix,
    For lcao to transform to pw, the psi_in would be the one in lcao basis,
    )
*/
template<typename T, typename Device>
class RepOut_PW : public RepOut<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
        constexpr static const Device * ctx = {};
    public:
        RepOut_PW();
        ~RepOut_PW();
        /// @brief wrapped gemm function, in parameter list the number of row is directly the C++ context, calculates AB = C
        /// @param mat1 A
        /// @param mat2 B
        /// @param mat3 C
        /// @param nrow1 number of rows of A
        /// @param ncol1 number of columns of A
        /// @param nrow2 number of rows of B
        /// @param ncol2 number of columns of B
        /// @param nrow3 number of rows of C
        /// @param ncol3 number of columns of C
        void matmul(
            T* mat1, T* mat2, T* mat3,
            const int nrow1, const int ncol1, 
            const int nrow2, const int ncol2,
            const int nrow3, const int ncol3
        );

        void project(psi::Psi<T, Device>* psi_in,
                     psi::Psi<T, Device>* psig,
                     psi::Psi<T, Device>* psi_out) override;

};

#endif // REPOUT_PW_H
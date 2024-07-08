#ifndef ABACUS_HS_MATRIX_K_HPP
#define ABACUS_HS_MATRIX_K_HPP

#include "module_basis/module_ao/parallel_orbitals.h"

#include <vector>
namespace hamilt
{
    template <typename TK>
    class HS_Matrix_K
    {
        public:
            HS_Matrix_K(const Parallel_Orbitals* paraV, bool no_s=false){
                hk.resize(paraV->nloc);
                if(!no_s) 
                {
                    sk.resize(paraV->nloc);
                }
                this->pv = paraV;
            }
            TK* get_hk() {return hk.data();}
            TK* get_sk() {return sk.data();}
            int get_size() {return hk.size();}
            void set_zero_hk() {hk.assign(hk.size(), 0);}
            void set_zero_sk() {sk.assign(sk.size(), 0);}
            const Parallel_Orbitals* get_pv() const {return this->pv;}
        private:
            std::vector<TK> hk;
            std::vector<TK> sk;
            const Parallel_Orbitals* pv = nullptr;
    };
}

#endif
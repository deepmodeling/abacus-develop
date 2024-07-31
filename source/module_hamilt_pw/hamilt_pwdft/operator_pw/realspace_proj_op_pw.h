#include <string>
#include <vector>
#include "operator_pw.h"
#include "module_base/macros.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis_k.h"

namespace hamilt {

    #ifndef DFTUTEMPLATE_H
    #define DFTUTEMPLATE_H
    template<class T> class DFTU: public T {};
    #endif
    template<typename T, typename Device>
    class DFTU<OperatorPW<T, Device>> : public OperatorPW<T, Device>
    {
        private:
            using Real = typename GetTypeReal<T>::type;
        public:
            // various constructors, support different types of projectors

            DFTU(const std::vector<int>& isk,
                 const UnitCell* ucell_in,      //< the dependency on UnitCell
                 const ModulePW::PW_Basis_K* pw_basis);

            DFTU(const std::vector<int> isk,
                 const std::vector<int>& l_hubbard,
                 const std::vector<double>& u,
                 const std::vector<double>& rgrid,
                 const std::vector<std::vector<double>>& projs,
                 const std::vector<int>& natom,                         //< UnitCell::nat
                 const std::vector<ModuleBase::Vector3<double>*>& tau,  //< UnitCell::...
                 const double& omega,                                   //< UnitCell::omega
                 const double& tpiba,                                   //< UnitCell::tpiba
                 const std::vector<ModuleBase::Vector3<double>>& q,     //< PW_Basis::getgpluskcar
                 const double& dq = 0.01,                               //< GlobalV::DQ
                 const double& nq = 10000);                             //< GlobalV::NQX

            virtual void act(const int nbands,
                             const int nbasis,
                             const int npol,
                             const T* tmpsi_in,
                             T* tmhpsi,
                             const int ngk = 0) const override;

            static void read_abacus_orb(const std::string& forb, 
                                        int& lmax, 
                                        std::vector<int>& nzeta,
                                        int& nr, 
                                        double& dr,
                                        std::vector<std::vector<double>>& radials);
                                    
    };
}
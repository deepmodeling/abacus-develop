#ifndef HYDROGEN_RADIALS_H_
#define HYDROGEN_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"
#include "module_base/assoc_laguerre.h"

class HydrogenRadials : public RadialSet
{
    public:
        HydrogenRadials() {}
        
        HydrogenRadials& operator=(const HydrogenRadials& rhs);
        HydrogenRadials* clone() const { return new HydrogenRadials(*this); } // covariant return type
        
        ~HydrogenRadials() {}

        void build(const int itype = 0,
                   const double charge = 1.0,
                   const int nmax = 0,
                   const double rcut = 10.0,
                   const double dr = 0.01,
                   const int rank = 0,
                   std::ofstream* ptr_log = nullptr        
        );
    private:
        void generate_hydrogen_radials(const double charge,
                                       const int nmax,
                                       const double rcut,
                                       const double dr,
                                       const int rank,
                                       std::ofstream* ptr_log);
        Assoc_Laguerre assoc_laguerre_;
};
#endif // HYDROGEN_RADIALS_H_
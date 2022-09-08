#ifndef __NONLOCALPW
#define __NONLOCALPW

#include "operator_pw.h"

#include "module_cell/unitcell_pseudo.h"

#include "src_pw/VNL_in_pw.h"

namespace hamilt
{

#ifndef __NONLOCALTEMPLATE
#define __NONLOCALTEMPLATE

template<class T> class Nonlocal : public T
{};

#endif

template<>
class Nonlocal<OperatorPW> : public OperatorPW
{
    public:
    Nonlocal(
        const int* isk_in,
        const pseudopot_cell_vnl* ppcell_in,
        const UnitCell_pseudo* ucell_in
    );

    virtual ~Nonlocal(){};

    virtual void init(const int ik_in)override;

    virtual void act
    (
        const psi::Psi<std::complex<double>> *psi_in, 
        const int n_npwx, 
        const std::complex<double>* tmpsi_in, 
        std::complex<double>* tmhpsi
    )const override;

    private:
    void add_nonlocal_pp(std::complex<double> *hpsi_in, const std::complex<double> *becp, const int m) const;

    mutable int max_npw = 0;

    mutable int npw = 0;

    mutable int npol = 0;

    const int* isk = nullptr;

    const pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell_pseudo* ucell = nullptr;
};

} // namespace hamilt

#endif
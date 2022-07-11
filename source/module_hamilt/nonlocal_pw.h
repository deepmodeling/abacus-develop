#ifndef __NONLOCALPW
#define __NONLOCALPW

#include "operator.h"

#include "module_cell/unitcell_pseudo.h"

#include "src_pw/VNL_in_pw.h"

namespace hamilt
{

class NonlocalPW : public Operator
{
    public:
    NonlocalPW(
        int max_npw_in,
        int npol_in,
        const int* ngk_in,
        const int* isk_in,
        const pseudopot_cell_vnl* ppcell_in,
        const UnitCell_pseudo* ucell_in
    );

    void init(const int ik)override;

    void act(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const override;

    private:
    void add_nonlocal_pp(std::complex<double> *hpsi_in, const std::complex<double> *becp, const int m) const;

    int max_npw = 0;

    int npol = 0;
    
    const int* ngk = nullptr;

    const int* isk = nullptr;

    const pseudopot_cell_vnl* ppcell = nullptr;

    const UnitCell_pseudo* ucell = nullptr;
};

} // namespace hamilt

#endif
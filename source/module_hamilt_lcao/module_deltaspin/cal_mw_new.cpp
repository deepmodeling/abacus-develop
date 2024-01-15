#include "spin_constrain.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_elecstate/elecstate_lcao.h"
#include "module_basis/module_ao/ORB_gen_tables.h"

template <>
void SpinConstrain<std::complex<double>, psi::DEVICE_CPU>::cal_MW_new()
{
    std::cout << "hello world" << std::endl;
    const ORB_gen_tables& uot = ORB_gen_tables::get_const_instance();
    for (int ipol = 0; ipol < this->npol_; ++ipol)
    {
        hamilt::HContainer<double>* dmR_current = dynamic_cast<const elecstate::ElecStateLCAO<std::complex<double>>*>(this->pelec)
            ->get_DM()->get_DMR_pointer(ipol+1);
    }
}
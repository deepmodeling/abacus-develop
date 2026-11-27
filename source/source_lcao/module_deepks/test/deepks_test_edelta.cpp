#include "deepks_test_runner.h"

#include "source_io/module_parameter/parameter.h"
#include "source_lcao/module_deepks/deepks_basic.h"

#include <complex>
#include <iomanip>

template <typename T>
void test_deepks<T>::check_edelta(std::vector<torch::Tensor>& descriptor)
{
    DeePKS_domain::load_model("model.ptg", ld.model_deepks);
    ld.allocate_V_delta(ucell.nat, kv.get_nkstot());
    if (PARAM.inp.deepks_equiv)
    {
        DeePKS_domain::cal_edelta_gedm_equiv(ucell.nat,
                                             this->ld.deepks_param,
                                             descriptor,
                                             this->ld.model_deepks,
                                             this->ld.gedm,
                                             this->ld.E_delta,
                                             0);
    }
    else
    {
        DeePKS_domain::cal_edelta_gedm(ucell.nat,
                                       this->ld.deepks_param,
                                       this->ld.model_deepks,
                                       this->ld.E_delta,
                                       descriptor,
                                       this->ld.pdm,
                                       this->ld.gedm);
    }

    std::ofstream ofs("E_delta.dat");
    ofs << std::setprecision(10) << this->ld.E_delta << std::endl;
    ofs.close();
    this->assert_file_matches_reference("E_delta.dat", "E_delta_ref.dat");
}

template void test_deepks<double>::check_edelta(std::vector<torch::Tensor>& descriptor);
template void test_deepks<std::complex<double>>::check_edelta(std::vector<torch::Tensor>& descriptor);

template <typename T>
void run_deepks_unit_edelta(test_deepks<T>& test)
{
    std::vector<torch::Tensor> descriptor;
    DeepksTestRunner::build_edelta(test, descriptor);
}

template void run_deepks_unit_edelta<double>(test_deepks<double>& test);
template void run_deepks_unit_edelta<std::complex<double>>(test_deepks<std::complex<double>>& test);

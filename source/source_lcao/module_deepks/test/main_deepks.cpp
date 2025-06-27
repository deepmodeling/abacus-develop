#include "LCAO_deepks_test.h"
#ifdef __MPI
#include <mpi.h>
#endif

int calculate();

template <typename T>
void run_tests(test_deepks<T>& test);

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    int status = calculate();
#ifdef __MPI
    MPI_Finalize();
#endif

    if (status > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int calculate()
{
    std::ifstream ifs("INPUT");
    char word[80];
    bool gamma_only_local;
    ifs >> word;
    ifs >> gamma_only_local;
    ifs.close();

    if (gamma_only_local)
    {
        test_deepks<double> test;
        run_tests(test);
        return test.failed_check;
    }
    else
    {
        test_deepks<std::complex<double>> test;
        run_tests(test);
        return test.failed_check;
    }
}

template <typename T>
void run_tests(test_deepks<T>& test)
{
    test.preparation();

    test.check_dstable();
    test.check_phialpha();

    test.check_pdm();

    std::vector<torch::Tensor> descriptor;
    test.check_descriptor(descriptor);

    torch::Tensor gdmx;
    test.check_gdmx(gdmx);
    test.check_gvx(gdmx);

    torch::Tensor gdmepsl;
    test.check_gdmepsl(gdmepsl);
    test.check_gvepsl(gdmepsl);

    test.check_orbpre();

    test.check_vdpre();
    test.check_vdrpre();

    test.check_edelta(descriptor);
    test.check_e_deltabands();
    test.check_f_delta_and_stress_delta();
    test.check_o_delta();

    std::cout << " [  ------  ] Total checks : " << test.total_check << std::endl;
    if (test.failed_check > 0)
    {
        std::cout << "\e[1;31m [  FAILED  ]\e[0m Failed checks : " << test.failed_check << std::endl;
    }
    else
    {
        std::cout << "\e[1;32m [  PASS    ]\e[0m All checks passed!" << std::endl;
    }
}

template void run_tests(test_deepks<double>& test);
template void run_tests(test_deepks<std::complex<double>>& test);
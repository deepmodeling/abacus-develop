#include "LCAO_deepks_test.h"
#ifdef __MPI
#include <mpi.h>
#endif

int calculate();

int main(int argc, char **argv)
{
#ifdef __MPI
	MPI_Init(&argc,&argv);
#endif
    int status = calculate();
#ifdef __MPI
	MPI_Finalize();
#endif

	if(status>0)
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
	test_deepks test;

	test.preparation();

	test.check_dstable();
	test.check_psialpha();

	test.check_pdm();
	test.check_gdmx();

	test.check_descriptor();
	test.check_gvx();

	test.check_edelta();
	test.check_v_delta();
	test.check_e_deltabands();
	test.check_f_delta();

	std::cout << " [  ------  ] Total checks : " << test.total_check <<std::endl;
	if(test.failed_check>0)
	{
		std::cout << "\e[1;31m [  FAILED  ]\e[0m Failed checks : " << test.failed_check <<std::endl;
	}
	else
	{
		std::cout << "\e[1;32m [  PASS    ]\e[0m All checks passed!" << std::endl;
	}

    return test.failed_check;
}
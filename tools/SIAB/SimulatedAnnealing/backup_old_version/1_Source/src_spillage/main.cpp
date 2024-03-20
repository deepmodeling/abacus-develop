#include "tools.h"
#include "read_INPUT.h"
#include "ReadData.h"
#include "MultiZeta.h"
#include "Plot_Psi.h"
#include "../src_parallel/parallel_global.h"
#include "../src_pw/memory_calculation.h"

int main(int argc, char **argv)
{
	cout << " ***************** " << endl;
	cout << "  WELCOME TO SIA! " << endl;
	cout << " ***************** " << endl;
	cout << " Generating Numerical Orbitals via minimizing spillage program start." << endl;
	// (1) recording start time.
    time_t time_start = time(NULL);
	timer::start();
    cout <<" Start  Time : " << ctime(&time_start);

	// (2) if parallel, prepare.
	cout << " First read in the parameters from INPUT." << endl;
	Parallel_Global::read_mpi_parameters(argc,argv);
	
	// (2) init class input.
	input.init();

	if(USEPW)
	{
		cout << " Init plane wave:" << endl;
		PW.init();
	}

    mz.set_levels();
    mz.cal_lmax_nmax();

	if(USEPW && CALSPI==1)
	{
		PW.table(); // mohan add 2010-06-16
	}

    input.Coef.init(
        NTYPE,
		LMAXUSED,
		NMAXUSED,
        ECUT,
        ECUT_JLQ,
        RCUT,
        TOLERENCE );

    ofs_running << "\n Coefficients init done." << endl;

    mz.metro.init();

    ofs_running << "\n Metropolis  init done." << endl;


    // mohan add 2009-09-11
    // mohan remove 2009-09-25
    // mohan get back 2009-10-10
    mz.metro.Psi2.init(input.ifs, NE,
                          LMAXUSED, RCUT, NCHIUSED,
                          TOLERENCE);


	// (3) calculate spillgae
	if( CALSPI==1 )
	{
		mz.init();
		// output "full" information.
		if( input.output && MY_RANK==0)
		{
			// output orbitals for ploting.
			Plot_Psi plot;
			plot.radial_wave_function();
		}
	}
	else if( CALSPI == 0)
	{
		Memory mem;
		mem.calculation();
	}

	// (4) output time.
	timer::finish();
    time_t time_finish = time(NULL);
    cout <<"\n Start  Time : " << ctime(&time_start);
    cout <<" Finish Time : " << ctime(&time_finish);
	cout <<" Total  Time : " << time_finish - time_start << " Seconds" << endl;

#ifdef __MPI
	MPI_Finalize();
#endif
	return 0;
}

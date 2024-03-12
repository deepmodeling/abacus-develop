#include "Metropolis.h"
#include "../src_parallel/parallel_common.h"
#include "Calculate_C4.h"
#include "../src_tools/Simulated_Annealing.h"
#include "Simulated_Annealing_Orbital.h"

Metropolis::Metropolis()
{
	prob_old = 0.0;
	prob_new = 0.0;
	temperature = 0.0;
	random_c4number = 0;
	output_each_istep = 0;
	kin_last = new double[1];
	kin_temp = new double[1];
}


Metropolis::~Metropolis()
{
	delete[] kin_last;	
	delete[] kin_temp;
}

// read in Metropolis parameters.
void Metropolis::init(void)
{
    if(test==1) TITLE("Metropolis","init");

	bool begin = false;

	if(MY_RANK==0)
	{
    	if(SCAN_BEGIN(input.ifs, "<METROPOLIS>"))
    	{
			begin = true;
		}
	}

#ifdef __MPI
	Parallel_Common::bcast_bool(begin);	
#endif

	ofs_running << "\n begin read in Metropolis parameter." << endl;
        
	if(begin)
	{
		if(MY_RANK==0)
		{	
        	READ_VALUE( input.ifs, this->spi_ini_temp);
        	READ_VALUE( input.ifs, this->spi_rate);
        	READ_VALUE( input.ifs, this->spi_ntemp);
        	READ_VALUE( input.ifs, this->spi_nsteps);

        	READ_VALUE( input.ifs, this->kin_ini_temp);
        	READ_VALUE( input.ifs, this->kin_rate);
        	READ_VALUE( input.ifs, this->kin_ntemp);
        	READ_VALUE( input.ifs, this->kin_nsteps);

			READ_VALUE( input.ifs, this->delta_kappa);
			READ_VALUE( input.ifs, this->output_each_istep);

			READ_VALUE( input.ifs, ACCEPTANCE_STEPS ); // mohan add 2009-10-31
			READ_VALUE( input.ifs, ACCEPTANCE_HIGH );
			READ_VALUE( input.ifs, ACCEPTANCE_LOW );

			READ_VALUE( input.ifs, KINETIC_MAX ); // mohan add 2009-10-31
			READ_VALUE( input.ifs, KINETIC_DR); //mohan add 2010-04-12
			READ_VALUE( input.ifs, OPTIMIZE_METHOD); // mohan add 2010-04-14
            //SCAN_END(input.ifs, "</METROPOLIS>");
		}
		
#ifdef __MPI
		MPI_Barrier(MPI_COMM_WORLD);

		Parallel_Common::bcast_double(this->spi_ini_temp);
		Parallel_Common::bcast_double(this->spi_rate);
		Parallel_Common::bcast_int(this->spi_ntemp);
		Parallel_Common::bcast_int(this->spi_nsteps);

		Parallel_Common::bcast_double(this->kin_ini_temp);
		Parallel_Common::bcast_double(this->kin_rate);
		Parallel_Common::bcast_int(this->kin_ntemp);
		Parallel_Common::bcast_int(this->kin_nsteps);

		Parallel_Common::bcast_double(this->delta_kappa);
		Parallel_Common::bcast_int(this->output_each_istep);

		Parallel_Common::bcast_int(ACCEPTANCE_STEPS );
		Parallel_Common::bcast_double(ACCEPTANCE_HIGH );
		Parallel_Common::bcast_double(ACCEPTANCE_LOW );

		Parallel_Common::bcast_double(KINETIC_MAX );
		Parallel_Common::bcast_double(KINETIC_DR );
		Parallel_Common::bcast_int(OPTIMIZE_METHOD);
#endif
		
    }

	this->temperature = this->spi_ini_temp;
   
   	// check whether the parameters are in the range.
	assert(spi_ini_temp >= 0.0);
	assert(spi_rate >= 0.0);
    assert(spi_ntemp >= 0);
	assert(spi_nsteps >= 0);

	assert(kin_ini_temp >= 0.0);
	assert(kin_rate >= 0.0);
    assert(kin_ntemp >= 0);
	assert(kin_nsteps >= 0);

	assert(delta_kappa >= 0.0);
	assert(output_each_istep>0);
	
    return;
}

// be called in MultiZeta::init(), reset the start temperature.
void Metropolis::reset_temperature(const int &istep)
{
	temperature = this->spi_ini_temp;
	for(int i=0; i<istep; i++)
	{
		temperature *= 0.8;
		// I think 0.8 is more reasonable than 0.1, mohan note, 2009-08-20
	}
	return;
}

void Metropolis::move_various_temperature( SpillageStep &Level)
{
	timer::tick("Metropolis","move_various_t");

	// there are two parts of this algorithm,
	// first part is about spillage, 
	// the second part is about kinetic energy.
	
	Simulated_Annealing_Orbital_Spillage sap(this,Level);
	sap.set_temperature_num(this->spi_ntemp)
	   .set_step_num(this->spi_nsteps);	
	try
	{
		sap.cal_all_temperature();
	}
	catch(File_Finish_Metropolis()) { return; }
	
	Simulated_Annealing_Orbital_Kinetic sak(this,Level);
	sak.set_temperature_num(this->kin_ntemp)
	   .set_step_num(this->kin_nsteps);	
	try
	{
		sak.cal_all_temperature();
	}
	catch(File_Finish_Metropolis()) { return; }
	
	timer::tick("Metropolis","move_various_t");
	return;	
}
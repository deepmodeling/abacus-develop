#include "Simulated_Annealing.h"
#include "mathzone.h"
#include "Random.h"

void Simulated_Annealing::cal_all_temperature(void)
{
	preprocess_all_temperature();
	for( size_t itemperature(0); itemperature != this->temperature_num; ++itemperature)
	{
		cal_each_temperature(itemperature);
		if(quit_loop_all_temperature()){ break; }
		update_temperature(itemperature);
	}
	reprocess_all_temperature();
	return;
}

void Simulated_Annealing::cal_each_temperature( size_t itemperature)
{
	preprocess_each_temperature(itemperature);
	for( size_t istep(0); istep != this->step_num; ++istep)
	{
		cal_each_step(itemperature, istep);
		if(quit_loop_each_temperature(itemperature)){ break; }
	}
	reprocess_each_temperature(itemperature);
	return;
}

void Simulated_Annealing::cal_each_step( size_t itemperature, size_t istep)
{
	preprocess_each_step(itemperature, istep);
	move_variable(itemperature, istep);
	double new_function;
	bool accept_flag(false), reject_flag(false);
	judge_accept_reject(new_function, itemperature, istep, accept_flag, reject_flag);
	process_accept_reject(new_function, itemperature, istep, accept_flag, reject_flag);
	reprocess_each_step(itemperature, istep);
	return;
}

void Simulated_Annealing::judge_accept_reject( double &new_function, size_t itemperature, size_t istep, bool &accept_flag, bool &reject_flag)
{
	if(must_accept_before(new_function, itemperature, istep)) { accept_flag = true; }
	else if(must_reject_before(new_function, itemperature, istep)) { reject_flag = true; }
	if(accept_flag || reject_flag){ return; }
	else
	{
		new_function = cal_new_function(itemperature, istep);
		if(must_accept_after(new_function, itemperature, istep)) { accept_flag = true; }
		else if(must_reject_after(new_function, itemperature, istep)) { reject_flag = true; }
	}

	if(accept_flag)
	{
		if(reject_flag) { throw Accept_Reject_Conflict(); }
		else { return; }
	}
	else if(reject_flag) { return; }
	else
	{
		if(new_function<=old_function) { accept_flag = true; }
		else
		{
			const double p = Mathzone::Boltzmann_dist(new_function - old_function, this->temperature);
			const double r = Random::between0and1();
			if (r < p) { accept_flag = true; }
			else { reject_flag = true; }		
		}
	}	
}

void Simulated_Annealing::process_accept_reject( double &new_function, size_t itemperature, size_t istep, bool &accept_flag, bool &reject_flag)
{
	// accept and reject may conflict without checking before
	if(accept_flag)
	{
		accept_process( new_function, itemperature, istep); 
	}
	if(reject_flag)
	{
		reject_process( new_function, itemperature, istep);		
	}
}

void Simulated_Annealing::accept_process( double new_function, size_t itemperature, size_t istep)
{
	this->old_function = new_function;
}

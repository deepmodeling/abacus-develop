#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include<cstddef>
#include<exception>
using std::exception;

class Accept_Reject_Conflict: public exception{};

class Simulated_Annealing
{
	public:
	
	virtual void cal_all_temperature(void);
	virtual void cal_each_temperature( size_t itemperature);
	virtual void cal_each_step( size_t itemperature, size_t istep);
	virtual void judge_accept_reject( double &new_function, size_t itemperature, size_t istep, bool &accept_flag, bool &reject_flag);
	virtual void process_accept_reject( double &new_function, size_t itemperature, size_t istep, bool &accept_flag, bool &reject_flag);
	virtual void accept_process( double new_function, size_t itemperature, size_t istep);

	virtual void move_variable(size_t itemperature, size_t istep)=0;
	virtual double cal_new_function(size_t itemperature, size_t istep)=0;
	
	virtual void preprocess_all_temperature(){}
	virtual bool quit_loop_all_temperature(){return false;}
	virtual void update_temperature(size_t itemperature){}
	virtual void reprocess_all_temperature(){}
	
	virtual void preprocess_each_temperature(size_t itemperature){}
	virtual bool quit_loop_each_temperature(size_t itemperature){return false;}
	virtual void reprocess_each_temperature(size_t itemperature){}
	
	virtual void preprocess_each_step(size_t itemperature, size_t istep){}
	virtual void reprocess_each_step(size_t itemperature, size_t iste){}
	
	virtual bool must_accept_before(double new_function, size_t itemperature, size_t istep){return false;}
	virtual bool must_reject_before(double new_function, size_t itemperature, size_t istep){return false;}
	virtual bool must_accept_after(double new_function, size_t itemperature, size_t istep){return false;}
	virtual bool must_reject_after(double new_function, size_t itemperature, size_t istep){return false;}
	
	virtual void reject_process(double new_function, size_t itemperature, size_t istep){}
	
	protected:
	size_t temperature_num;
	size_t step_num;
	double temperature;
	double old_function;
	
	public:
	Simulated_Annealing &set_temperature_num(size_t temperature_num_in){temperature_num = temperature_num_in; return *this;}
	Simulated_Annealing &set_step_num(size_t step_num_in){step_num = step_num_in; return *this;}
	Simulated_Annealing &set_temperature(double temperature_in){temperature = temperature_in; return *this;}
	Simulated_Annealing &set_old_function(double old_function_in){old_function = old_function_in; return *this;}
};

#endif
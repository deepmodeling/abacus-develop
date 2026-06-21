#ifndef SIMULATED_ANNEALING_ORBITAL_H
#define SIMULATED_ANNEALING_ORBITAL_H

class Metropolis;

class Metropolis::Simulated_Annealing_Orbital: public Simulated_Annealing
{
	public:
	Simulated_Annealing_Orbital( Metropolis *metropolis_in, SpillageStep &Level_in)
		:Simulated_Annealing(),
		 metropolis(metropolis_in),
		 Level(Level_in){}
	protected:
	virtual void preprocess_all_temperature()
	{
		metropolis->preprocess_all_temperature(Level);
	}
	virtual void preprocess_each_temperature(size_t itemperature)
	{
		metropolis->Psi2.ofso << "spillage_ntemp = " << itemperature << " kappa = " << metropolis->kappa << endl;	
	}
	virtual void update_temperature(size_t itemperature)
	{
		metropolis->update_t();		// change temperature			
	}
	virtual void preprocess_each_step( size_t itemperature, size_t istep)
	{
		metropolis->small_jql();
	}
	virtual void cal_each_step( size_t itemperature, size_t istep)
	{
		preprocess_each_step_beyondloop(itemperature, istep);
		// for each (atom_type,L,N), we need a set of C4
		// to describe the radial wave functions.
		// after each radial wave function is update(ne parameters),
		// we call it 'a step under a particular temperature is done'			
		for (int ic = 0; ic < metropolis->nchi; ic++)
		{
			metropolis->info.set_ic(ic).unpack_ic(Level);
			const int T = Level.wayc[ic].type;		// get the type index of this radial wave function.
			if(Level.info[T].state=="skip") continue;	// if "skip", the C4 will not changed in this 'Big step' // mohan add 2009-01-27
			for (int ie = 0; ie < Level.ne; ie++)
			{
				metropolis->info.set_ie(ie);
				try
				{
					this->Simulated_Annealing::cal_each_step(itemperature, istep);
				}
				catch(Small_Jlq &e) { continue; }
			}	
		}				
		reprocess_each_step_beyondloop(itemperature, istep);
	}
	virtual void preprocess_each_step_beyondloop(size_t itemperature, size_t istep){}
	virtual void reprocess_each_step_beyondloop(size_t itemperature, size_t istep)
	{
		metropolis->ofs_1(istep, this->Level);
		metropolis->norm_C4(istep, this->Level, this->step_num);			
	}	
	virtual void move_variable(size_t itemperature, size_t istep)
	{
		metropolis->trial_c4();
	}
	virtual void accept_process( double new_function, size_t itemperature, size_t istep)
	{			
		metropolis->accept_process(this->Level);
	}
	virtual void reject_process( double new_function, size_t itemperature, size_t istep)
	{
		metropolis->reject_process();
	}		
	Metropolis *metropolis;
	SpillageStep &Level;
};

class Metropolis::Simulated_Annealing_Orbital_Spillage: public Simulated_Annealing_Orbital
{
	public:
	Simulated_Annealing_Orbital_Spillage( Metropolis *metropolis_in, SpillageStep &Level_in)
		:Simulated_Annealing_Orbital(metropolis_in, Level_in){}
	private:
	virtual void preprocess_each_temperature(size_t itemperature)
	{
		metropolis->init_spi_states(itemperature);
		Simulated_Annealing_Orbital::preprocess_each_temperature(itemperature);
	}
	virtual void reprocess_each_temperature(size_t itemperature)
	{
		metropolis->cout_1(itemperature);
	}
	virtual void preprocess_each_step_beyondloop(size_t itemperature, size_t istep)
	{
		metropolis->file_finish_metropolis(istep);
		metropolis->change_accept_rate(istep,Level);
		metropolis->info.set_update_number(0);	// each step we renew every Jlq one by one.			
	}
	virtual void move_variable(size_t itemperature, size_t istep)
	{
		Simulated_Annealing_Orbital::move_variable(itemperature, istep);
		set_temperature( metropolis->temperature );
	}
	virtual double cal_new_function(size_t itemperature, size_t istep)
	{
		set_old_function(input.SV.cal_defined_value(0));
		return metropolis->cal_spillage(this->Level);
	}
};

class Metropolis::Simulated_Annealing_Orbital_Kinetic: public Simulated_Annealing_Orbital
{
	public:
	Simulated_Annealing_Orbital_Kinetic( Metropolis *metropolis_in, SpillageStep &Level_in)
		:Simulated_Annealing_Orbital(metropolis_in, Level_in){}
	private:
	virtual void preprocess_each_temperature(size_t itemperature)
	{
		metropolis->init_kin_states();			
		Simulated_Annealing_Orbital::preprocess_each_temperature(itemperature);
		metropolis->init_temperature_kinetic(itemperature, this->Level);
	}
	virtual void reprocess_each_temperature(size_t itemperature)
	{
		metropolis->cout_1(itemperature);
		metropolis->cout_2(itemperature, this->Level);
	}
	virtual void preprocess_each_step_beyondloop(size_t itemperature, size_t istep)
	{
		metropolis->file_finish_metropolis(istep);
		metropolis->reset_high_kinetic_energy(itemperature, istep);
		metropolis->change_accept_rate(istep, this->Level);	
		metropolis->info.update_number = 0;	// each step we renew every Jlq one by one.			
	}
	virtual void move_variable(size_t itemperature, size_t istep)
	{
		Simulated_Annealing_Orbital::move_variable(itemperature, istep);
		set_temperature( metropolis->kin_temp[metropolis->info.ic] );
	}		
	virtual bool must_reject_before(double new_function, size_t itemperature, size_t istep)
	{
		// make sure the new spillage value is in the range.
		// otherwise, we don't even calculate the kinetic energy.
		return metropolis->cal_spillage(this->Level) >= metropolis->min_spillage * metropolis->kappa;
	}
	virtual double cal_new_function(size_t itemperature, size_t istep)
	{
		set_old_function(metropolis->kin_last[metropolis->info.ic]);
		return metropolis->cal_kinetic_energy(this->Level);
	}	
	virtual void accept_process( double new_function, size_t itemperature, size_t istep)
	{
		Simulated_Annealing_Orbital::accept_process( new_function, itemperature, istep);
		metropolis->kin_last[metropolis->info.ic] = new_function;
	}			
};

#endif
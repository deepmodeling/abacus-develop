#include "Metropolis.h"
#include "../src_parallel/parallel_common.h"
#include "Calculate_C4.h"

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

	// firstly, make the initia kappa value.
	// mohan add 2009-11-12
	// this 0.001 must be added to the second step
	// to avoid the oscillation of first orbital. 
	this->kappa = 1.0 + this->delta_kappa; 
	
	assert( this->spi_ntemp>=0 );

	this->nchi = Level.nchi;

	delete[] kin_last;
	this->kin_last = new double[nchi];	
	ZEROS(kin_last, nchi);

	delete[] kin_temp;
	this->kin_temp = new double[nchi];
	ZEROS(kin_temp, nchi);

	// how many temperature we try.
	// there are two parts of this algorithm,
	// first part is about spillage, 
	// the second part is about kinetic energy.
	for (int itemp = 0; itemp < this->spi_ntemp + this->kin_ntemp; itemp++)
	{ 
		// if the kinetic energy is not small enough, we try new kappa again.
		//try_new_kappa_again:
		
		if(itemp < this->spi_ntemp)
		{
			ofs_running << "\n" << " ---> SpillageTemp " << setw(3) << itemp+1 << endl; 
			ofs_running  << "      Temperature  " << setw(6) << temperature << endl;
			this->states = "Spi";
			ofs_running  << setw(5) << "STEPS" 
				<< setw(20) << "SPILLAGE"
				<< setw(10) << "UPDATES" << endl;

			if(itemp==0)
			{
				cout << setw(5) << "STEP" << setw(10) << "TEMP"
				<< setw(20) << "SPILLAGE" << endl;
			}
		}
		else
		{
			// Please DIY, mohan 2010-04-14
			if(OPTIMIZE_METHOD==1)
			{
				this->states="Kin";
			}
			else if(OPTIMIZE_METHOD==2)
			{
				this->states = "Ecut";
			}
			else
			{
				WARNING_QUIT("move_various_temperature","Check OPTIMIZE_METHOD");			
			}	
		}
		
		this->Psi2.ofso << "spillage_ntemp = " << itemp << " kappa = " << kappa << endl;	
		
		// set the initial temperature and kinetic energy value.
		if(itemp== this->spi_ntemp)
		{
			this->min_spillage = input.SV.cal_defined_value(0); // save the spillage value.
			ofs_running << "\n\n Spillage at final temperature " << min_spillage << endl;
			
			for(int ic=0; ic<this->nchi; ic++)
			{
				const int T = Level.wayc[ic].type;
				const int L = Level.wayc[ic].L;
				const int N = Level.wayc[ic].N;
				
				double *c4_last = new double[Level.ne];						
				for(int ie2=0; ie2<Level.ne; ie2++) 
				{
					c4_last[ie2] = input.Coef.C4_old(T,L,N,ie2);
				}
						
				this->kin_last[ic] = 0.0;

				if(states=="Kin")
				{	
					for(int ie2=0; ie2<Level.ne; ie2++)
					{
						this->kin_last[ic] += pow( c4_last[ie2]*(ie2+1)*PI/RCUT,2 )
						*this->Psi2.jjnorm(L,ie2);
					}
				}
				else if(states=="Ecut")
				{
					this->kin_last[ic] = this->Psi2.get_ecut( c4_last, L );
				}

				// set the start temperature of kinetic part of Metropolis algorithms.
				this->kin_temp[ic] = this->kin_last[ic] * this->kin_ini_temp; 

				ofs_running << " Orbital " << ic+1 << " Initial E_kin " << setw(15) << this->kin_last[ic] << endl;
						
				delete[] c4_last;
			}
		}
		
		
		// mohan add 2009-09-11
		// how many steps moving in this temperature.
		int nsteps_now = 0;
		if( states == "Spi" ) nsteps_now = this->spi_nsteps;
		else nsteps_now = this->kin_nsteps;
			 

		for (int istep = 0; istep < nsteps_now; istep++)
		{
			// P.S.: a method to quit this whole program.
			// you only need to create a file name : FINISH
			// another function: if there existep a file 'FINISH'
			// in the current dir, the program quit at the first
			// iteration.
			// Here I choose a proper parameter: 500 steps,
			// you can change on your own.
			if(istep%500==0)
			{
				// mohan add 2009-08-28	
				ifstream ifs("FINISH");
				if( ifs ) 
				{
					return;
				}
			}

			// because the initia temperature may be too high
			// which caused by the high kinetic energy, so
			// we need to reset the temperature.
			if(itemp == this->spi_ntemp && istep == 50)
			{
				double high = 0.0;
				for(int ic=0; ic<this->nchi; ic++)
				{
					this->kin_temp[ic] = this->kin_last[ic] * this->kin_ini_temp; 
					high = std::max( kin_temp[ic], high );
				}

				for(int ic=0; ic<this->nchi; ic++)
				{
					this->kin_temp[ic] = high; 
					ofs_running << " Notice: Change temperature of orbital " << ic+1 << " to " << kin_temp[ic] << endl;
				}
			} 

			// Here is the part about changing Metropolis algorithm's accept rate.
			// we want to calculate the accept rate every 'ACCEPTANCE_STEPS' steps.
			// but we don't want to calculate at istep=0;
			if(istep%ACCEPTANCE_STEPS == 0)
			{
				if(istep>0)
				{
					for(int ic=0; ic<this->nchi; ic++)
					{
						const int T = Level.wayc[ic].type;
						const int L = Level.wayc[ic].L;
						const int N = Level.wayc[ic].N;
						for(int ie2=0; ie2<Level.ne; ie2++)
						{
							assert( input.Coef.accept_number(T,L,N,ie2) <= ACCEPTANCE_STEPS);
							
							double rate = (double)input.Coef.accept_number(T,L,N,ie2)/(double)ACCEPTANCE_STEPS;
							
							if(rate > ACCEPTANCE_HIGH)
							{
								input.Coef.accept_rate(T,L,N,ie2) *= 2.0; // You can DIY 0.5
							}
							else if(rate < ACCEPTANCE_LOW)
							{
								 input.Coef.accept_rate(T,L,N,ie2) *= 0.5; // You can DIY 2.0
							}

							// output information
							//if(ie2%5==0) cout << "\n";
							//cout << setw(8) << rate*100 << "%";
							//cout << setw(8) << input.Coef.accept_rate(T,L,N,ie2);
						}
					}
				}
				input.Coef.accept_number.zero_out();
			} // end accept_rate


			//========================================================
			// if we want to move step, we need the set of C4 from CS
			// in single zeta, the C4 can be either read in or random,
			// in double zeta, the C4 should start from random, in
			// most cases.
			//========================================================
			this->move_one_step(itemp, istep, Level, nsteps_now);		

		}// end istep 
			
		cout << setw(5) << itemp+1
		<< setiosflags(ios::fixed)
		<< setprecision(3)
		<< setiosflags(ios::scientific) 
		<< setiosflags(ios::showpoint)
		<< setw(10) 
		<< temperature
		<< resetiosflags(ios::fixed)
		<< setiosflags(ios::fixed)
		<< setiosflags(ios::showpoint)
		<< setprecision(10)
		<< setw(19) 
		<< input.SV.cal_defined_value(0)*100 << "%" << endl;

	
		//-------------------------------------------
		// it's about minimizing kinetic energy.
		//-------------------------------------------
		if( itemp >= this->spi_ntemp)
		{
			for(int ic=0; ic<this->nchi; ic++)
			{
				cout << setw(5) << Level.wayc[ic].spd;
				cout << setiosflags(ios::scientific)
					<< setiosflags(ios::showpoint)
					<< setprecision(3)
					<< setw(10)<< this->kin_temp[ic];
				cout << setprecision(10)
					<< resetiosflags(ios::fixed) 
					<< setiosflags(ios::fixed) 
					<< setiosflags(ios::showpoint)
					<< setw(20)<< this->kin_last[ic] << endl;
			}
	
			for(int ic=0; ic<this->nchi; ic++)
			{
				if(this->kin_last[ic] > KINETIC_MAX)
				{
					cout << "\n kin_last[" << ic << "]=" << kin_last[ic] << " KINETIC_MAX = " << KINETIC_MAX << endl;
					cout << " There maybe 2 reasons, first: The starting temperature for kinetical energy is too small." << endl;
					cout << " Second: The kappa is too small." << endl;
					//WARNING_QUIT("Minimizing Kinetic Energy","Kinetical Energy Too Large!");
                                        cout << "Kinetical Energy Too Large!"<<endl;
				}
				//this->kappa += this->delta_kappa;
				//goto try_new_kappa_again;
			}
		}

		// change temperature
		this->update_t();
	}// end item
	timer::tick("Metropolis","move_various_t");
	return;
}


void Metropolis::move_one_step(
	const int &itemp, 
	const int &istep, 
	SpillageStep &Level,
	const int &nsteps_now)
{	
	timer::tick("Metropolis", "move_one_step");

	//ofs_running << "\n move_one_step";

	// each step we renew every Jlq one by one.
	int update_number = 0;
	
	static int *min_os = new int[this->nchi];
	
	// for each (atom_type,L,N), we need a set of C4
	// to describe the radial wave functions.
	// after each radial wave function is update(ne parameters),
	// we call it 'a step under a particular temperature is done'
	for (int ic = 0; ic < this->nchi; ic++)
	{
		// get the type index of this radial wave function.
		const int T = Level.wayc[ic].type;

		// mohan add 2009-01-27
		// if "skip", the C4 will not changed in this 'Big step'
		if(Level.info[T].state=="skip") continue;

		const int L = Level.wayc[ic].L;
		const int N = Level.wayc[ic].N;
		
		for (int ie = 0; ie < Level.ne; ie++)
		{
			//cout << "ic=" << ic << " ie=" << ie << endl;
			// mohan add 2009-08-26
			// this judgement means, 
			// If users don't want to use Jlq(ie)
			// which is too oscillation(large), they
			// can use only small part of Jlq.
			// mohan add BLOCK_NE_MIN 2009-08-27
			if( (ie<BLOCK_NE_MIN) || (ie>=BLOCK_NE)) continue;

			// trial a new value of c4(T,L,N).
			input.Coef.trial_c4( T, L, N, ie);

			// calculate the spillage according to the new c4.
			for (int is = 0; is < STRNUM; is++)
			{
				input.SV.value[is] = Level.get_spillage(is, ic, ie);
			}

			// flag1 is about spillage can be accepted.
			// flag2 is about kinetic energy can be accepted.
			static bool accept_flag1 = false;
			static bool accept_flag2 = false;
					
			//================================================================
			// STEP 1 :
			//================================================================
			if( states == "Spi" )
			{
				// 0 means get old value, 1 means get new value.
				accept_flag1 = this->accept( input.SV.cal_defined_value(0), input.SV.cal_defined_value(1), temperature );
				accept_flag2 = true;
			}
			//================================================================
			// STEP 2
			//================================================================
			else
			{
				// make sure the new spillage value is in the range.
				// otherwise, we don't even calculate the kinetic energy.
				if( input.SV.cal_defined_value(1) < this->min_spillage * this->kappa)
				{
					// spillage is in the range, condition one ok!
					accept_flag1 = true;
									
					double *c4tmp = new double[Level.ne];
					for(int ie2=0; ie2<Level.ne; ie2++) 
					{
						c4tmp[ie2] = input.Coef.C4(T,L,N,ie2);
					}

					double value2 = 0.0;
					if(states=="Kin")
					{
						// get the normalized parameters: c4
						this->Psi2.norm_c4( c4tmp, L );
					
						// and then calculate the new set of kinetic energy.
						for(int ie2=0; ie2<Level.ne; ie2++)
						{
							value2 += pow( c4tmp[ie2]*(ie2+1)*PI/RCUT,2 )
							*this->Psi2.jjnorm(L,ie2);//In Rydberg Unit.
						}
					}
					else if(states=="Ecut")
					{
						value2 = this->Psi2.get_ecut( c4tmp, L );
					}

					// judge here
					accept_flag2 = this->accept( this->kin_last[ic], value2, this->kin_temp[ic] );
					if(accept_flag2)
					{
						kin_last[ic] = value2;
					}
					
					delete[] c4tmp;	
				}
				else
				{
					accept_flag1 = false;
				}
			}

			if (accept_flag1 && accept_flag2)
			{
				++update_number;
				for (int is = 0; is< STRNUM; is++)
				{
					Level.updateQS( is ); // (1) update Soverlap and Qoverlap for each structure,
					// S matrix and Q matrix has changed due to the change
					// of c4
				}
				
				input.SV.update_value();// (2) update spillage value
			
				for(int is=0; is<STRNUM; is++)
				{
					Level.data[is].Mkb = Level.data[is].Mkb_used; // (2.5) mohan add 2010-05-02	
					
					/*
					for(int ik=0; ik<NKS; ik++)
					{
						for(int ib=0; ib<NBANDS; ib++)
						{
							cout << "\n Mkb(" << ik << "," << ib << ")=" << Level.data[is].Mkb(ik,ib); 
						}
					}
					cout << endl;
					int ok;
					cin >> ok;
					*/
				}
								
				input.Coef.update_c4(T, L, N, ie); // (3) update c4 

				/*
				// mohan add 2009-09-25
				// another magic number!!
				// after 500 steps in this temperature,
				// we begin to accumulate coefficients.				
				if(istep>500)
				{
					// give the first value.
					// just in  case some value will never change.
					static bool copyflag = false;
					if(!copyflag)
					{
						for(int k=0; k<input.Coef.C4.getSize(); k++)
						{
							input.Coef.C4_accumulate.ptr[k] = input.Coef.C4_old.ptr[k];
						}
						copyflag = true;
					}
					
					// call the subroutine, accumulating the
					// new accept C4(T,L,N,ie)
					input.Coef.accumulating_C4( T, L, N, ie);
				}
				*/
				
			}
			else
			{
				input.Coef.go_back_c4( T, L, N, ie );
			}
		}	
	}


	
	// the last part, easy part, output information!
	// selectly output spillage information
	//	if(update_number>0)
//	{
		if((istep)%output_each_istep==0)
		{
			if(states=="Spi")
			{
				ofs_running << setw(5) << istep+1
				<< setprecision(10)
				<< setiosflags(ios::fixed)
				<< setw(20) << input.SV.cal_defined_value(0)
				<< setw(10) << update_number << endl;

				// output each structure's spillage value.
				//	for(int is=0; is< input.str_num; is++)
				//	{
				//		cout << setw(10) << input.SV.value_old[is]*100 << "%";
				//	}
			}
			else if(states=="Kin"||states=="Ecut")
			{
				//<< " ----> " << states << "Temp " << setw(15)
				//<< istep+1 << "\n" ; 
//				<< "       SpiVal"
//				<< setprecision(6)
//				<< setiosflags(ios::fixed)
//				<< setw(15) << input.SV.cal_defined_value(0);
				for(int ic=0; ic< this->nchi; ic++)
				{ 
					ofs_running << setw(4) << istep+1 << Level.wayc[ic].spd;

					ofs_running <<  setprecision(3) 
					<< setiosflags(ios::scientific)
					<< setw(10)<< this->kin_temp[ic];

					ofs_running << setprecision(10) 
					<< setiosflags(ios::fixed) 
					<< setw(20)<< this->kin_last[ic] << endl;
				}
			}
		}
//	}

	// because we choose c4 randomly between -0.1~0.1
	// so we need to normalize the set of c4 after a few
	// steps, here I choose a proper number: 50
	// you can change on you own, but please notice it can't
	// be too small, because we need to update all Q and S
	// matrix after this normalization, which may affect the 
	// speed of the program. You can't choose too large,
	// because it makes the update not efficient. 
	if ( (istep % 50 == 0 && istep !=0 ) || istep == nsteps_now-1) 
	{
		// use to describe the real space mesh.
		// don't need to choose too accuray, because
		// it may affect the speed of running.
		// you can change on your own, too.
		double dr = 0.01; // a.u.
	
		// norm ic for all radial wave functions.
		// after this call, all c4 will change/update.
		// but the spillage will not change.
		// (because we just normalize the wavefunctions).
		Calculate_C4::norm_ic(
	    	input.Coef.C4,
	    	input.Coef.ntype,
	    	input.Coef.lmax,
	    	input.Coef.nmax,
	    	input.Coef.enumber,
	    	input.Coef.tolerence,
		  	input.Coef.rcut,
	  	  	dr);
		
		// get new set of C4, so we copy the set of C4
		// to C4_old. This will also not affect the 
		// spillage value.
		input.Coef.copy_c4_to_old();

		// for each structure, we need to update their
		// Q matrix and S matrix, which may be time
		// consuming.
		for (int is = 0; is < STRNUM; is++)
		{
			Level.init_QS_matrix(is);
		}
	}		
	timer::tick("Metropolis", "move_one_step");
	return;
}

double Metropolis::Boltzmann_dist(const double &value, const double &t)
{
	static const double small = 1.0e-20;
	assert(t > small);
	assert(value > 0.0);
//	cout << "\n value = " << value;
//	cout << "\n Boltzmann_dist = " << exp( -value/t );
	return exp(-value / t);
}

void Metropolis::update_t(void)
{
	if(this->states == "Spi")
	{
		this->temperature *= this->spi_rate;
//		cout << " New Temperature for spillage = " << temperature << endl;
	}
	else if(this->states == "Kin"||states=="Ecut")
	{
		// mohan add high 2009-11-05
		double high = 0.0;
		for(int ic=0; ic< this->nchi; ic++)
		{
			this->kin_temp[ic] *= this->kin_rate;
			high = std::max( this->kin_temp[ic], high );
		}

		for(int ic=0; ic< this->nchi; ic++)
		{
			this->kin_temp[ic] = high;
			ofs_running << setprecision(6) << setiosflags(ios::scientific) 
			<< " Temperature for orbital " << ic+1 << " = " << kin_temp[ic] << endl;
		}
		
		ofs_running << " new kappa=" << kappa << endl;
	}
	return;
}

bool Metropolis::accept( const double &v_old, const double &v_new, const double &temperature_in)
{
	assert(temperature_in > 0.0);
	this->prob_old = v_old;
	this->prob_new = v_new;

//	cout << "\n v_new = " << v_new;
//	cout << "\n v_old = " << v_old;
//	cout << "\n t = " << temperature_in;

	const double deltav = v_new - v_old;

	if ( deltav < 0.0 )
	{
		return 1; // accept
	}
	else if ( deltav > 0.0 )
	{
		const double p = Metropolis::Boltzmann_dist(v_new - v_old, temperature_in);
		const double r = Random::between0and1();

		if (r < p) return 1; // accept
		if (r >= p) return 0; // reject
	}

	return 0;
}

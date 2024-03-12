#include "Metropolis.h"
#include "tools.h"
#include "Calculate_C4.h"
#include <algorithm>
#include <cstddef>

void Metropolis::preprocess_all_temperature(SpillageStep &Level)
{
	this->nchi = Level.nchi;			

	// firstly, make the initia kappa value.
	// mohan add 2009-11-12
	// this 0.001 must be added to the second step
	// to avoid the oscillation of first orbital. 
	this->kappa = 1.0 + this->delta_kappa; 
	
	delete[] this->kin_last;
	this->kin_last = new double[nchi];	
	ZEROS(this->kin_last, nchi);

	delete[] this->kin_temp;
	this->kin_temp = new double[nchi];
	ZEROS(this->kin_temp, nchi);
}

void Metropolis::file_finish_metropolis(size_t istep)
{
	// P.S.: a method to quit this whole program.
	// you only need to create a file name : FINISH
	// another function: if there existep a file 'FINISH'
	// in the current dir, the program quit at the first
	// iteration.
	// Here I choose a proper parameter: 500 steps,
	// you can change on your own.
	const int check_finish_file(500);
	if( 0 == (istep%check_finish_file) )
	{
		// mohan add 2009-08-28	
		ifstream ifs("FINISH");
		if( ifs ) 
		{
			throw File_Finish_Metropolis();
		}
	}			
}

void Metropolis::update_t(void)
{
	if(this->info.states == "Spi")
	{
		this->temperature *= this->spi_rate;
//		cout << " New Temperature for spillage = " << temperature << endl;
	}
	else if(this->info.states == "Kin"||this->info.states=="Ecut")
	{
		// mohan add high 2009-11-05
		for(int ic=0; ic< this->nchi; ic++)
		{
			this->kin_temp[ic] *= this->kin_rate;
		}
		const double high = *max_element(this->kin_temp, this->kin_temp+this->nchi);
		for(int ic=0; ic< this->nchi; ic++)
		{
			this->kin_temp[ic] = high;
			ofs_running << setprecision(6) //<< setiosflags(ios::scientific) 
			<< " Temperature for orbital " << ic+1 << " = " << kin_temp[ic] << endl;
		}
		
		ofs_running << " new kappa=" << kappa << endl;
	}
	return;
}

void Metropolis::change_accept_rate(size_t istep, SpillageStep &Level)
{
	// Here is the part about changing Metropolis algorithm's accept rate.
	// we want to calculate the accept rate every 'ACCEPTANCE_STEPS' steps.
	// but we don't want to calculate at istep=0;
	if(istep%ACCEPTANCE_STEPS == 0)
	{
		if(istep>0)
		{
			for(int ic = 0; ic < this->nchi; ic++)
			{
				const int T = Level.wayc[ic].type;
				const int L = Level.wayc[ic].L;
				const int N = Level.wayc[ic].N;
				for(int ie2=0; ie2<Level.ne; ie2++)
				{
					assert( input.Coef.accept_number(T,L,N,ie2) <= ACCEPTANCE_STEPS);

					const double rate = (double)input.Coef.accept_number(T,L,N,ie2)/(double)ACCEPTANCE_STEPS;
					
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
}		

void Metropolis::reset_high_kinetic_energy(size_t itemperature,size_t istep)
{
	// because the initia temperature may be too high
	// which caused by the high kinetic energy, so
	// we need to reset the temperature.
	if( 0 == itemperature && 50 == istep )
	{
		for(int ic=0; ic<this->nchi; ic++)
		{
			this->kin_temp[ic] = this->kin_last[ic] * this->kin_ini_temp;
		}
		const double high = *max_element(this->kin_temp, this->kin_temp+this->nchi);
		for(int ic=0; ic<this->nchi; ic++)
		{
			this->kin_temp[ic] = high; 
			ofs_running << " Notice: Change temperature of orbital " << ic+1 << " to " << kin_temp[ic] << endl;
		}
	}	
}

void Metropolis::init_spi_states(size_t itemperature)
{
	ofs_running << "\n" << " ---> SpillageTemp " << setw(3) << itemperature+1 << endl; 
	ofs_running  << "      Temperature  " << setw(6) << this->temperature << endl;
	this->info.states = "Spi";
	ofs_running  << setw(5) << "STEPS" 
		<< setw(20) << "SPILLAGE"
		<< setw(10) << "UPDATES" << endl;

	if(0==itemperature)
	{
		cout << setw(5) << "STEP" << setw(10) << "TEMP"
		<< setw(20) << "SPILLAGE" << endl;
	}	
}

void Metropolis::init_kin_states()
{
	// if the kinetic energy is not small enough, we try new kappa again.
	//try_new_kappa_again:
	// Please DIY, mohan 2010-04-14
	if(OPTIMIZE_METHOD==1)
	{
		this->info.states="Kin";
	}
	else if(OPTIMIZE_METHOD==2)
	{
		this->info.states = "Ecut";
	}
	else
	{
		WARNING_QUIT("move_various_temperature","Check OPTIMIZE_METHOD");			
	}	
}

void Metropolis::init_temperature_kinetic(size_t itemperature, SpillageStep &Level)
{
	// set the initial temperature and kinetic energy value.			
	if(0==itemperature)
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

			if(this->info.states=="Kin")
			{	
				for(int ie2=0; ie2<Level.ne; ie2++)
				{
					this->kin_last[ic] += pow( c4_last[ie2]*(ie2+1)*PI/RCUT,2 )
					*this->Psi2.jjnorm(L,ie2);
				}
			}
			else if(this->info.states=="Ecut")
			{
				this->kin_last[ic] = this->Psi2.get_ecut( c4_last, L );
			}

			// set the start temperature of kinetic part of Metropolis algorithms.
			this->kin_temp[ic] = this->kin_last[ic] * this->kin_ini_temp; 

			ofs_running << " Orbital " << ic+1 << " Initial E_kin " << setw(15) << this->kin_last[ic] << endl;
					
			delete[] c4_last;
		}				
	}
}

void Metropolis::norm_C4(size_t istep, SpillageStep &Level, size_t nsteps_now)
{
	// because we choose c4 randomly between -0.1~0.1
	// so we need to normalize the set of c4 after a few
	// steps, here I choose a proper number: 50
	// you can change on you own, but please notice it can't
	// be too small, because we need to update all Q and S
	// matrix after this normalization, which may affect the 
	// speed of the program. You can't choose too large,
	// because it makes the update not efficient. 
	if ( (istep % 50 == 0 && istep !=0 ) || istep == nsteps_now-1) 		// Peize Lin test
//	if ( (istep % 1 == 0 && istep !=0 ) || istep == nsteps_now-1) 		// Peize Lin test
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
			input.Coef.accept_rate,			// Peize Lin update 2015-12-24
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
		for (int istr(0); istr < STRNUM; ++istr)
		{
			Level.init_QS_matrix(istr);
		}
	}	
}

void Metropolis::small_jql()
{
	// mohan add 2009-08-26
	// this judgement means, 
	// If users don't want to use Jlq(ie)
	// which is too oscillation(large), they
	// can use only small part of Jlq.
	// mohan add BLOCK_NE_MIN 2009-08-27
	if( (this->info.ie<BLOCK_NE_MIN) || (this->info.ie>=BLOCK_NE)) 
		throw Small_Jlq();			
}

void Metropolis::trial_c4()
{
	//========================================================
	// if we want to move step, we need the set of C4 from CS
	// in single zeta, the C4 can be either read in or random,
	// in double zeta, the C4 should start from random, in
	// most cases.
	//========================================================
	const int ie = this->info.ie;
	const int T = this->info.T;
	const int L = this->info.L;
	const int N = this->info.N;	
	// trial a new value of c4(T,L,N).
	input.Coef.trial_c4( T, L, N, ie);	
}

double Metropolis::cal_spillage(SpillageStep &Level)
{
	const int ie = this->info.ie;
	const int ic = this->info.ic;
	for ( size_t istr = 0; istr < STRNUM; ++istr)
	{
		input.SV.value[istr] = Level.get_spillage(istr, ic, ie);
	}
	return input.SV.cal_defined_value(1);
}

double Metropolis::cal_kinetic_energy(SpillageStep &Level)
{
	double *c4tmp = new double[Level.ne];
	
	const int T = this->info.T;
	const int L = this->info.L;
	const int N = this->info.N;		
	for(int ie2=0; ie2<Level.ne; ie2++) 
	{
		c4tmp[ie2] = input.Coef.C4(T,L,N,ie2);
	}

	double value2 = 0.0;
	if(this->info.states=="Kin")
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
	else if(this->info.states=="Ecut")
	{
		value2 = this->Psi2.get_ecut( c4tmp, L );
	}	
	delete[] c4tmp;	
	return value2;
}

void Metropolis::accept_process(SpillageStep &Level)
{
	++this->info.update_number;
	
	for (int istr = 0; istr< STRNUM; ++istr)
	{
		Level.updateQS( istr ); // (1) update Soverlap and Qoverlap for each structure,
		// S matrix and Q matrix has changed due to the change of c4
	}

	input.SV.update_value();// (2) update spillage value

	for(int istr=0; istr<STRNUM; istr++)
	{
		Level.data[istr].Mkb = Level.data[istr].Mkb_used; // (2.5) mohan add 2010-05-02	
		
		/*
		for(int ik=0; ik<NKS; ik++)
		{
			for(int ib=0; ib<NBANDS; ib++)
			{
				cout << "\n Mkb(" << ik << "," << ib << ")=" << Level.data[istr].Mkb(ik,ib); 
			}
		}
		cout << endl;
		int ok;
		cin >> ok;
		*/
	}
	
	const int ie = this->info.ie;
	const int T = this->info.T;
	const int L = this->info.L;
	const int N = this->info.N;	
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

void Metropolis::reject_process()
{
	const int ie = this->info.ie;
	const int T = this->info.T;
	const int L = this->info.L;
	const int N = this->info.N;
	input.Coef.go_back_c4( T, L, N, ie );	
}

void Metropolis::ofs_1(size_t istep, SpillageStep &Level)
{
	if(0==(istep%this->output_each_istep))
	{
		if(this->info.states=="Spi")
		{
			ofs_running << setw(5) << istep+1
			<< setprecision(10)
			<< setiosflags(ios::fixed)
			<< setw(20) << input.SV.cal_defined_value(0)
			<< setw(10) << this->info.update_number << endl;

			// output each structure's spillage value.
			//	for(int is=0; is< input.str_num; is++)
			//	{
			//		cout << setw(10) << input.SV.value_old[is]*100 << "%";
			//	}
		}
		else if(this->info.states=="Kin"||this->info.states=="Ecut")
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
				<< setw(10)<< this->kin_temp[ic]
				<< resetiosflags(ios::scientific);

				ofs_running << setprecision(10) 
				<< setiosflags(ios::fixed) 
				<< setw(20)<< this->kin_last[ic] << endl
				<< resetiosflags(ios::fixed);
			}
		}
	}
}

void Metropolis::cout_1(size_t itemperature)
{
	cout << setw(5) << itemperature+1
	<< setiosflags(ios::showpoint)	
	<< setiosflags(ios::scientific)
	<< setprecision(3)
	<< setw(10) 
	<< this->temperature
	<< resetiosflags(ios::scientific)
	<< setiosflags(ios::fixed)
	<< setprecision(10)
	<< setw(19) 
	<< input.SV.cal_defined_value(0)*100 << "%" 
	<< resetiosflags(ios::fixed)
	<< resetiosflags(ios::showpoint)		
	<< endl;
}

void Metropolis::cout_2(size_t itemperature, SpillageStep &Level)
{
	//-------------------------------------------
	// it's about minimizing kinetic energy.
	//-------------------------------------------
	for(int ic=0; ic<this->nchi; ic++)
	{
		cout << setw(5) << Level.wayc[ic].spd;
		cout << setiosflags(ios::showpoint)
			<< setiosflags(ios::scientific)
			<< setprecision(3)
			<< setw(10)<< this->kin_temp[ic]
			<< resetiosflags(ios::scientific);
		cout << setprecision(10)
			<< setiosflags(ios::fixed) 
			<< setw(20)<< this->kin_last[ic] 
			<< resetiosflags(ios::fixed) 
			<< resetiosflags(ios::showpoint)
			<< endl;
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


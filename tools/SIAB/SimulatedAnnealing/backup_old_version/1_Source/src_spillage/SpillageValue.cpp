#include "SpillageValue.h"
#include "tools.h"

SpillageValue::SpillageValue()
{
	dim = 0;
	value = new double[1];
	value_old = new double[1];
}

SpillageValue::~SpillageValue()
{
	delete[] value;
	delete[] value_old;
}

void SpillageValue::allocate( const int &dim_in )
{
	this->dim = dim_in;
	assert( dim < 100);
	assert( dim > 0);
	
	delete[] value;
	delete[] value_old;
	value = new double[dim];
	value_old = new double[dim];

	// assume max step < 100
	const int max_levels = 100;
	value_each_level.create(dim, max_levels);
	
	this->reset();

	return;
}

void SpillageValue::reset(void)
{
	for(int i=0; i<dim; i++)
	{
		value[i] = 1.0;
		value_old[i] = 1.0;
	}
	return;
}

double SpillageValue::cal_defined_value(bool get_new_flag)
{
	// spillage value of each structure.
	double* v;
	if( get_new_flag ) 
	{
		v = this->value;
	}
	else 
	{
		v = this->value_old;
	}
	
	double mv=0.0;
	// case 1: // get the largest spillage value.
	if( abs(SCHEME_VALUE)==1)
	{
		for(int i=0; i<dim; i++)
		{
			mv = std::max( mv, v[i] );
		}
		return mv;
	}
	// case 2: // get average spillage value.
	else if( abs(SCHEME_VALUE)==2)
	{
		for(int i=0; i<dim; i++)
		{
			mv += v[i]/dim;
		}
		return mv;
	}
	// case 3: // mixture of case 1 and 2
	else if( abs(SCHEME_VALUE)==3)
	{
		double average=0.0;
		for(int i=0; i<dim; i++)
		{
			average += v[i]/dim;
		}

		double max_diverse = 0.0;
		double avg_diverse = 0.0;
		for(int i=0; i<dim; i++)
		{
//			max_diverse = std::max(max_diverse, abs(average-v[i]));
			avg_diverse += abs(average-v[i])/dim;
		}

		if(average < 0.1)
		{
			mv = average * (1.0 + avg_diverse / (average*0.1));
		}
		else
		{
			mv = average;
		}
		return mv;
	}
	else
	{
		cout << "\n Please Check SCHEME_VALUE." << endl;
		WARNING_QUIT("SpillageValue","get_value");
	}
}

void SpillageValue::update_value(void)
{
    for(int i=0; i<dim; i++)
    {
        value_old[i] = value[i];
    }
	return;
}

void SpillageValue::out(void)
{
	for(int i=0; i<dim; i++)
	{
		cout << setprecision(6);
		cout << " Structure" << setw(3) << i+1 << setw(10) << value_old[i]*100 << "%" << endl;
	}
}

void SpillageValue::save_level( const int &ilevel)
{
	// dim stands for number of structures.
	for(int i=0; i<dim; i++)
	{
		this->value_each_level(i ,ilevel) = value_old[i];	
	}
	return;
}

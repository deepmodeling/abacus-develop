// mohan mddify 2010-05-02
#ifndef SPILLAGE_VALUE_H
#define SPILLAGE_VALUE_H

#include "common.h"

class SpillageValue
{
	friend class Read_INPUT;
	
	friend class MultiZeta;
	// in multizeta, 'value_old' are calculated.

	friend class Metropolis;
	// in metropolis, 'value' are set.

	friend class Out_Orbital;
	// in OurOrbital, 'value each level' are used.

	public:

	// 1: get new spillage value.
	// 0: get old spillage value.
	double cal_defined_value( bool get_new_flag);
	
	void update_value();
	
	void out();
	
	void save_level( const int &ilevel);
	 
	private:

	// spillage value for each structure.
	double *value;
	
	// the previous spillage value for each structure.
	double *value_old;

	// spillage for all structures.
	matrix value_each_level; 

	// dimension is the number of structures.
	// eg. 5 different dimers.
	int dim;

	void allocate( const int &dim_in );
	
	void reset();	
	
	SpillageValue();
	~SpillageValue();
};

#endif

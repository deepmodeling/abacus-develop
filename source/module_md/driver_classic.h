#ifndef DRIVER_CLASSIC_H
#define DRIVER_CLASSIC_H

#include "../module_cell/unitcell_pseudo.h"

class Driver_classic
{
	public:
	
	Driver_classic();
	~Driver_classic();

	static void init();

	private:

	// reading the parameters
	static void reading();

    // convert INPUT parameters for classic MD
    static void convert(UnitCell_pseudo &ucell_c);

	// classic MD
	static void classic_world();

};

#endif
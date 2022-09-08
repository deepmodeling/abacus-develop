#ifndef DRIVER_H
#define DRIVER_H

class Driver
{
	public:
	
	Driver();
	~Driver();

	void init();

	private:

	// reading the parameters
	void reading();

	// do stuff, have fun!
	void atomic_world();

	// interface to the actual calculations
	void driver_run();

};

#endif

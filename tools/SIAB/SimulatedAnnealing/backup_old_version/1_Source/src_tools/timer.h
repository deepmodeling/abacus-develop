//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06
//==========================================================
#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;

//==========================================================
// CLASS :
// NAME : timer(calculate time)
//==========================================================
class timer
{
	timer();
	~timer();

	public:
//==========================================================
// MEMBER FUNCTIONS :
// NAME : tick(use twice at a time)
// NAME : start
// NAME : finish
// NAME : enable
// NAME : disable
// NAME : print
// NAME : print_all
//==========================================================
	static void tick(const string &class_name_in,const string &name_in);

	static void start(void);
	static void finish(void);

	static void enable(void);
	static void disable(void);

	static void print_all(void);
	static long double print_until_now(void);
	static double print(const string &name_in);

	private:

//==========================================================
// MEMBER VARIABLES : 
// NAME : disabled(if disabled , timer can't work)
// NAME : n_clock(uplimit number of clock numbers)
// NAME : n_now(the index of clocks) 
// NAME : start_flag
// NAME : cpu_start
// NAME : cput_second(record cpu time)
// NAME : name(record clock name)
// NAME : class_name(record another name for clock)
//==========================================================
	static bool disabled;
	static int n_clock;
	static int n_now;
	static int start_flag;
	static double *cpu_start;
	static double *cpu_second;
	static string *name;
	static string *class_name;
	static int *calls;

//==========================================================
// MEMBER FUNCTIONS :
// NAME : cpu_time(calculate time)
//==========================================================
	static double cpu_time(void);
	static bool delete_flag;

};
#endif


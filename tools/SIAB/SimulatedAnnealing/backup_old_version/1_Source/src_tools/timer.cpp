//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06
//==========================================================
#include "timer.h"

using namespace std;

//----------------------------------------------------------
// EXPLAIN :   
//----------------------------------------------------------
bool timer::disabled = false;
int timer::n_clock = 100;
int timer::n_now = 0;
int timer::start_flag = -1;
double* timer::cpu_second;
double* timer::cpu_start;
string* timer::name;
string* timer::class_name;
int* timer::calls;
bool timer::delete_flag=false;

timer::timer()
{
};

timer::~timer()
{
};

void timer::finish()
{
	timer::tick("","total");
	print_all();
	if(delete_flag)
	{
		delete[] cpu_second;
		delete[] cpu_start;
		delete[] name;
		delete[] class_name;
		delete[] calls;
	}
}

//----------------------------------------------------------
//
//----------------------------------------------------------
void timer::start(void)
{
	cpu_second = new double[n_clock]();
	cpu_start = new double[n_clock]();
	name = new string[n_clock]();
	class_name = new string[n_clock]();
	calls = new int[n_clock]();
	delete_flag = true;

	for (int i = 0; i < n_clock; i++)
	{
		cpu_start[i] = (double)start_flag;
		name[i]= "\0";
		class_name[i]= "\0";
		calls[i] = 0;
	}

	// first init ,then we can use tick
	timer::tick("","total");
	return;
}

double timer::cpu_time(void)
{
	clock_t t1 = 0;
//----------------------------------------------------------
// EXPLAIN : here static is important !!
// only first call can let t0 = 0,clock begin
// when enter this function second time , t0 > 0
//----------------------------------------------------------
	static clock_t t0 = 0;
	if (t0 == 0) 
	{
		t0 = clock();
	}
	t1 = clock() - t0;

	if(t1 < 0) return 0; // mohan add, abandon the cross point time 2^32 ~ -2^32 .
	else return (double)t1/CLOCKS_PER_SEC;
}

void timer::tick(const string &tagc,const string &tag)
{
//----------------------------------------------------------
// EXPLAIN : if timer is disabled , return
//----------------------------------------------------------
	if (disabled) 
	{
		return;
	}

	int find_clock=0;//counter
//----------------------------------------------------------
// EXPLAIN :  find if the tag has been used 
//----------------------------------------------------------
	for(find_clock=0; find_clock<n_now; find_clock++)
	{
		if (tag == name[find_clock] && tagc == class_name[find_clock]) 
		{
			break;
		}
	}

//----------------------------------------------------------
// EXPLAIN : if it's a new tag; add a new tag to list;
// add this new tag to name list;
//----------------------------------------------------------
	if (find_clock == n_now)
	{
		n_now++;
		name[find_clock] = tag;
		class_name[find_clock] = tagc;
	}

//----------------------------------------------------------
// EXPLAIN : if exceed the the uplimits of list 
//----------------------------------------------------------
	if (n_now >= n_clock) 
	{
		cout << "\nError! Too many timer!";
		return;
	}

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : cpu_time
//
// EXPLAIN : start_flag is minus 1, 
// so if cpu_start == start_flag,means a new clock counting
// begin, hence we record the start time of this clock 
// counting , if cpu_start != start_flag, means it's
// the end of this counting, so we add the time during
// this two 'time point'  to the clock time storage.
//----------------------------------------------------------
	if (cpu_start[find_clock] == start_flag )
	{
		cpu_start[find_clock] = cpu_time();
		calls[find_clock]++;
	}
	else
	{
		cpu_second[find_clock] += cpu_time() - cpu_start[find_clock];
		cpu_start[find_clock] = start_flag;
	}
	return;
}



void timer::enable(void)
{
	disabled = false;
	return;
}

void timer::disable(void)
{
	disabled = true;
	return;
}

double timer::print(const string &tag)
{
	int index = 0;
	for(int i=0; i<n_now; i++)
	{
		if (tag == name[i]) 
		{
			index = i;
			break;
		}
	}
//	cout << "\n " << name[i] 
//		 << " time : " 
//		 << cpu_second[i] << " (sec)" << endl;
	return cpu_second[index];
}

long double timer::print_until_now(void)
{
	// stop the clock
	timer::tick("","total");
	// start again
	timer::tick("","total");
	return print("total");
}

void timer::print_all()
{
//	cout<<"\n timer::print_all()"<<endl;
	const double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	cout<<"\n--------- CLASS NAME--------------- NAME ----- TIME(sec) -- CALLS ----- AVG ----- PER%";
	bool *print_flag = new bool[n_clock];
	for(int i=0; i<n_clock; i++) print_flag[i] = false;
	
	for (int i=0; i<n_clock; i++)
	{
		int k = 0;
		double tmp = -1.0;
		for(int j=0; j<n_clock; j++)
		{
			if(print_flag[j]) 
			{
				continue;
			}
			if(tmp < cpu_second[j])
			{
				k = j;
		//		cout << "\n k = " << k;
				tmp = cpu_second[j];
			}
		}
		print_flag[k]=true;
	
		if ((cpu_second[k] >= 0 && cpu_second[k] < small) ||
		        (cpu_second[k] <= 0 && cpu_second[k] > -small))
		{
			continue;
		}
		const long double spend_time = cpu_second[k];
		const double average_spend_time = spend_time/calls[k];

		cout  << "\n" 
			 << setw(20) << class_name[k]
			 << setw(20) << name[k]
			 << setw(15) << spend_time
			 << setw(10) << calls[k]
			 << setw(10) << setprecision(2) << average_spend_time
			 << setw(10) << spend_time / cpu_second[0] * 100 << "%";
			
	}
	cout<<"\n--------------------------------------------------------------------------------------"<<endl;
	delete[] print_flag;
	return;
}

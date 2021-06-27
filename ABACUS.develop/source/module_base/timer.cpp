//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06,
// UPDATE : Peize Lin at 2019-11-21
//==========================================================
#include "timer.h"
#include<vector>

#ifdef __MPI
#include "mpi.h"
#endif

using namespace std;

//----------------------------------------------------------
// EXPLAIN :   
//----------------------------------------------------------
bool timer::disabled = false;
size_t timer::n_now = 0;
map<string,map<string,timer::Timer_One>> timer::timer_pool;

void timer::finish(ofstream &ofs,const bool print_flag)
{
	timer::tick("","total");
	if(print_flag)
		print_all( ofs );
}

//----------------------------------------------------------
//
//----------------------------------------------------------
void timer::start(void)
{
	// first init ,then we can use tick
	timer::tick("","total");
}

double timer::cpu_time(void)
{
//----------------------------------------------------------
// EXPLAIN : here static is important !!
// only first call can let t0 = 0,clock begin
// when enter this function second time , t0 > 0
//----------------------------------------------------------
	static clock_t t0 = clock();
	const clock_t t1 = clock() - t0;
	return (t1<0) ? 0 : (double)t1/CLOCKS_PER_SEC;
		// mohan add, abandon the cross point time 2^32 ~ -2^32 .
}

void timer::tick(const string &class_name,const string &name,char level_in)
{
//----------------------------------------------------------
// EXPLAIN : if timer is disabled , return
//----------------------------------------------------------
	if (disabled)
		return;
	
	Timer_One &timer_one = timer_pool[class_name][name];

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : cpu_time
//
// EXPLAIN :
// if start_flag == true,means a new clock counting begin,
// hence we record the start time of this clock counting.
// if start_flag == false, means it's the end of this counting,
// so we add the time during this two 'time point'  to the clock time storage.
//----------------------------------------------------------
	if(timer_one.start_flag)
	{
#ifdef __MPI
		timer_one.cpu_start = MPI_Wtime();
#else
		timer_one.cpu_start = cpu_time();
#endif
		++timer_one.calls;
		timer_one.level = level_in;
		timer_one.start_flag = false;
	}
	else
	{
#ifdef __MPI
		timer_one.cpu_second += MPI_Wtime() - timer_one.cpu_start;
#else
		timer_one.cpu_second += cpu_time() - timer_one.cpu_start;
#endif
		timer_one.start_flag = true;
	}
}

long double timer::print_until_now(void)
{
	// stop the clock
	timer::tick("","total");
	// start again
	timer::tick("","total");
	return timer_pool[""]["total"].cpu_second;
}

void timer::print_all(ofstream &ofs)
{
	constexpr double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6
	
	vector<pair<pair<string,string>,Timer_One>> timer_pool_order;
	for(auto &timer_pool_A : timer_pool)
	{
		const string class_name = timer_pool_A.first;
		for(auto &timer_pool_B : timer_pool_A.second)
		{
			const string name = timer_pool_B.first;
			const Timer_One timer_one = timer_pool_B.second;
			if(timer_pool_order.size() < timer_one.order+1)
				timer_pool_order.resize(timer_one.order+1);
			timer_pool_order[timer_one.order] = {{class_name, name}, timer_one};
		}
	}
	
	cout << setprecision(2);
	ofs << setprecision(3);
	cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << endl;
	ofs <<"\n\n\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << endl;
	for(auto &timer_pool_order_A : timer_pool_order)
	{
		const string &class_name = timer_pool_order_A.first.first;
		const string &name = timer_pool_order_A.first.second;
		const Timer_One &timer_one = timer_pool_order_A.second;
		
		if(timer_one.cpu_second < small)
			continue;
		if(timer_one.level > 'X')
			continue;
		
		ofs  << " " 
			 << setw(2)  << timer_one.level
			 << setw(20) << class_name
			 << setw(20) << name
			 << setw(15) << timer_one.cpu_second
			 << setw(10) << timer_one.calls
			 << setw(10) << setprecision(2) << timer_one.cpu_second/timer_one.calls
			 << setw(10) << timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << endl;

		cout << resetiosflags(ios::scientific);
		
		cout << " " 
			 << setw(2)  << timer_one.level
			 << setw(20) << class_name
			 << setw(20) << name
			 << setw(15) << setprecision(5) << timer_one.cpu_second
			 << setw(10) << timer_one.calls
			 << setw(10) << setprecision(2) << timer_one.cpu_second/timer_one.calls
			 << setw(10) << timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << endl;		
	}
	cout<<" ----------------------------------------------------------------------------------------"<<endl;
	ofs <<" ----------------------------------------------------------------------------------------"<<endl;
}

/*
void timer::print_all(ofstream &ofs)
{
//	cout<<"\n timer::print_all()"<<endl;
	const double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	cout << setprecision(2);

	// prepare
	bool *print_flag = new bool[n_clock];
	for(int i=0; i<n_clock; i++) 
	{
		print_flag[i] = false;
	}

	int type = 1; // 2:calls 1:total_time
	bool non_reorder = 1;
	
	cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << endl;
	ofs <<"\n\n\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << endl;
	ofs << setprecision(3);
	for (int i=0; i<n_clock; i++)
	{
		int k = 0;
		double tmp = -1.0;

		if(non_reorder)
		{
			k = i;
		}
		else
		{
			// search in all clocks
			for(int j=0; j<n_clock; j++)
			{
				if(print_flag[j]) 
				{
					continue;
				}
				if(type==1)
				{
					if(tmp < cpu_second[j])
					{
						k = j;
						tmp = cpu_second[j];
					}
				}
				else if(type==2)
				{
					if(tmp < calls[j])
					{
						k = j;
						tmp = calls[j];
					}
				}
			}
		}
		print_flag[k]=true;
	
		if ((cpu_second[k] >= 0 && cpu_second[k] < small) ||
		        (cpu_second[k] <= 0 && cpu_second[k] > -small))
		{
			continue;
		}

		if( level[k] > 'X' ) continue;


		const long double spend_time = cpu_second[k];
		const double average_spend_time = spend_time/calls[k];


		ofs  << " " 
			 << setw(2) << level[k]
			 << setw(20) << class_name[k]
			 << setw(20) << name[k]
			 << setw(15) << spend_time
			 << setw(10) << calls[k]
			 << setw(10) << setprecision(2) << average_spend_time
			 << setw(10) << spend_time / cpu_second[0] * 100 << "%" << endl;


		cout << resetiosflags(ios::scientific);

		cout  << " " 
		     << setw(2) << level[k]
			 << setw(20) << class_name[k]
			 << setw(20) << name[k]
			 << setw(15) << spend_time
			 << setw(10) << calls[k]
			 << setw(10) << setprecision(2) << average_spend_time
			 << setw(10) << spend_time / cpu_second[0] * 100 << "%" << endl;
			
	}
	cout<<" ----------------------------------------------------------------------------------------"<<endl;
	ofs <<" ----------------------------------------------------------------------------------------"<<endl;
	delete[] print_flag;
	return;
}
*/
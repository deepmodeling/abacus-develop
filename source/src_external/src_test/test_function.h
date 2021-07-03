#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H

#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<string>
#include<mpi.h>
#include<sys/time.h>

// Peize Lin add 2015-11-11
static void MPI_RANK_OFSTREAM( const std::string& file_name, std::stringstream &content )
{
	std::stringstream file;
	int my_rank;	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	file<<file_name<<"_"<<my_rank;

	std::ofstream ofs(file.str().c_str(), std::ofstream::app);
	ofs<<content.str();

	content.str("");
	ofs.close();
}

// Peize Lin add 2016-06-06
template<typename T1,typename T2>
static std::ostream & operator<<( std::ostream & os, const std::pair<T1,T2> &p )
{
	os<<"<"<<p.first<<","<<p.second<<">";
	return os;
}

// Peize Lin add 2016-06-06
template<typename T>
static std::ostream & operator<<( std::ostream & os, const std::vector<T> &v )
{
	os<<"[";
	for( const T &i : v )
		os<<i<<endl;
	os<<"]";
	return os;
}

// Peize Lin add 2016-06-06
template<typename T>
static std::ostream & operator<<( std::ostream & os, const std::set<T> &v )
{
	os<<"(";
	for( const T &i : v )
		os<<i<<"\t";
	os<<")";
	return os;
}

// Peize Lin add 2016-06-06
template<typename T1, typename T2>
static std::ostream & operator<<( std::ostream & os, const std::map<T1,T2> &v )
{
	for( const auto &i : v )
//		os<<"{"<<i.first<<":"<<i.second<<"}"<<"\t";
		os<<i.first<<endl<<i.second<<endl;
	return os;
}

// Peize Lin add 2016-10-10
static double cal_time( const timeval &t_begin )
{
	timeval t_end;
	gettimeofday( &t_end, NULL);
	return (double)(t_end.tv_sec-t_begin.tv_sec) + (double)(t_end.tv_usec-t_begin.tv_usec)/1000000.0;
}

// Peize Lin add 2016-10-10
static double cut_time( timeval &t )
{
	const double time = cal_time(t);
	gettimeofday( &t, NULL);
	return time;
}

// Peize Lin add 2019-12-12
static vector<double> get_memory(const int N)
{
	vector<double> m;
	string s;
	ifstream ifs("/proc/meminfo");
	ifs>>s>>s;
	m.push_back(stoi(s)/1024);
	for(int i=1; i<N; ++i)
	{
		ifs>>s>>s>>s;
		m.push_back(stoi(s)/1024);
	}
	ifs.close();
	return m;
}

#endif		// TEST_FUNCTION_H

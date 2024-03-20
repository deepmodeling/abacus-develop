#ifndef METROPOLIS_INFO_H
#define METROPOLIS_INFO_H

struct Metropolis_Info
{
	int ie;
	int ic;
	int T;
	int L;
	int N;
	int update_number;
	string states; // mohan add 2009-10-15	

	Metropolis_Info &set_ie(int ie_in){ ie = ie_in; return *this;}
	Metropolis_Info &set_ic(int ic_in){ ic = ic_in; return *this;}
	Metropolis_Info &set_update_number(int update_number_in){ update_number = update_number_in; return *this;}

	void unpack_ic(SpillageStep &Level)
	{
		T = Level.wayc[ic].type;
		L = Level.wayc[ic].L;
		N = Level.wayc[ic].N;		
	}
};

class File_Finish_Metropolis: public exception{};

class Small_Jlq: public exception{};

#endif
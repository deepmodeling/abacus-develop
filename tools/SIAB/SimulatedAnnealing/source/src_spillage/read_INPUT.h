#ifndef READ_INPUT_H
#define READ_INPUT_H

#include "common.h"
#include "Coefficients.h"
#include "ReadData.h"
#include "SpillageValue.h"

// read information from file : INPUT
class Read_INPUT
{
	public:
	Read_INPUT();
	~Read_INPUT();

	ReadData* QS_data;
	Coefficients Coef;
	SpillageValue SV;

	ifstream ifs;
	ifstream ifs2;
	bool cal_sp;
	bool output;

	// get this value from Multi_Zeta.cpp
	int nlevel;

	void init(void);

	void read_test(void);
	void read_plot_c4(void);
//	void read_cal_c4(void);
	void read_start_c4(const int &step);

	private:
	int test;

	void bcast(void);
	void read_PW_QS_1(void);
	void read_PW_QS_2(void);
	void read_WFC(void);
	void read_WFC_2(void);

	string *qsfile; 

};

#endif

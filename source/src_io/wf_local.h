#ifndef WF_LOCAL_H
#define WF_LOCAL_H

#include "../src_pw/tools.h"

// mohan add 2010-09-09
namespace WF_Local
{
	void write_lowf(const string &name, double** ctot);
	void write_lowf_complex(const string &name, complex<double>** ctot, const int &ik);

	void distri_lowf(double** ctot, double **c);
	void distri_lowf_complex(complex<double>** ctot, complex<double> **cc);

	void distri_lowf_new(double** ctot, const int &is);
	void distri_lowf_complex_new(complex<double>** ctot, const int &ik);

	void distri_lowf_aug(double** ctot, double **c_aug);
	void distri_lowf_aug_complex(complex<double>** ctot, complex<double> **c_aug);

	int read_lowf(double **c, const int &is);

	int read_lowf_complex(complex<double> **c, const int &ik, const bool &newdm);
}

#endif

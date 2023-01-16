#ifndef CAL_TEST_H
#define CAL_TEST_H

namespace Cal_Test
{
	void test_memory(void);
	int cal_np(const double &ggcut, const int &n1, const int &n2, const int &n3);
	void print_mem(const int &nproc);

	extern double mporter;

	// about charge density.
	extern double mrho;
	extern double mrho_save;
	extern double mrho_core;

	// about pulay mixing.
	extern double mRrho;
	extern double mdRrho;
	extern double mdrho;
	extern double mrho_save2;

	// about potential on FFT grid.
	extern double mvltot;
	extern double mvr;
	extern double mvrs;
	extern double mvrs1;
	extern double mvnew;

	// about charge in g space.
	extern double mrhog;
	extern double mrhog_save;
	extern double mrhog_core;	
	
	// about others.
	extern double mhs;
	extern double mwf;
	extern double mnonzero;
	extern double mspar_hsrho;
	
	extern double mgvec;
	extern double mig2fftw;
	extern double mig2fftc;
	extern double mgg;
	extern double mig123;
	extern double mstrucFac;
	extern double meigts123;

	extern double mtot;

}
#endif

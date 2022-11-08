#ifndef EPSILON0_VASP_H
#define EPSILON0_VASP_H

#include "../src_pw/wavefunc.h"
#include "../src_parallel/parallel_global.h"

class Epsilon0_vasp 
{

public:
		
		Epsilon0_vasp();
		~Epsilon0_vasp();
		
		bool epsilon;
		double domega;
		int nomega;
		double eta;
		
		int oband;
		int uband;

		void Cal_psi_nu(int ik);
		
		void cal_epsilon0(const elecstate::ElecState* pelec);
		void Init();
		void Delete();
		void Cal_b(int ik);
		void Cal_psi(int ik);
		void Cal_psi_nabla(int ik);
		void Cal_epsilon0s(const elecstate::ElecState* pelec);
		void Cal_T();
		void Cal_epsilon0();
		void CMatrixMul(int m, int n, int l, std::complex<double>** A, std::complex<double>** B, std::complex<double>** C);
		
		std::complex<double> **eps0;

private:
		
		bool init_finish;
		std::complex<double> ***b;
		std::complex<double> **psi;
		std::complex<double> ***psi_nabla;
		std::complex<double> ***psi_nu;
		std::complex<double> **eps0s;
		std::complex<double> **T;
};

namespace GlobalC
{
extern Epsilon0_vasp epsilon0_vasp;
}

#endif	

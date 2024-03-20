#include "memory_calculation.h"

#include "../src_tools/complexmatrix_inline.h"

Memory::Memory()
{}

Memory::~Memory()
{}

void Memory::calculation(void)
{
	TITLE("Memory","calculation");
	
	const double dk = 0.01;
	int kmesh = static_cast<int>(sqrt(ECUT) / dk) +1 + 4;
	if (kmesh % 2 == 0)++kmesh;

	const int n1 = NTYPE;
	const int n2 = LMAXUSED + 1;
	const int n3 = NMAXUSED;
	const int n4 = kmesh; 
	const int n1234 = n1 * n2 * n3 * n4;

	cout << "\n psi1d dimension(NTYPE, LMAXUSED+1, NMAXUSED, kmesh) : " << n1234 
	<< " = " << n1 << " * " << n2 << " * "
	<< n3 << " * " << n4 << " = " << n1234 * 8 / (double)1024 / (double)1024 << " MB";

	for(int i=0; i<mz.nlevel; i++)
	{
		int nwfc2 = 0;
		for(int it=0; it<NTYPE; it++)
		{
			for(int l=0; l<mz.lmax_type[it]+1; l++)
			{
				nwfc2 += NA[it] * (2*l+1);
				cout << "\n na=" << NA[it] << " 2*l+1=" << 2*l+1 << " n=1";
			}
		}
		cout << "\n newfc2=" << nwfc2;
		cout << "\n psi3d dimension(NKSTOT=" << NKSTOT <<", nwfc2=" << nwfc2 << ", npwx=" << PW.npwx << ")="
		<< NKSTOT * nwfc2 * PW.npwx * 16 / (double)1024 / (double)1024 << " MB";

		cout << "\n S matrix dim = " << nwfc2 * nwfc2;
		
//		complex<double>* haha = new complex<double>[PW.npwx];
//		ZEROS(haha, PW.npwx);

		ComplexMatrix haha(3,PW.npwx);
	
		complex<double> sum = complex<double>(0,0);
		for(int t=0; t<nwfc2*nwfc2; t++)
		{
			timer::tick("Memory","test");
			for(int j=0; j<PW.npwx; j++)
			{
//				sum += conj(haha[j]) * haha[j];
				sum += conj(haha(2,j)) * haha(2,j);
			}
			timer::tick("Memory","test");
		}

	}

	return;
}

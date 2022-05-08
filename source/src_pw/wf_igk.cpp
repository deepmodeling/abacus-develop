#include "wf_igk.h"

#include "global.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

WF_igk::WF_igk()
{
	g2kin = new double[1];
}

WF_igk::~WF_igk()
{
	if (GlobalV::test_deconstructor)
	{
		std::cout << " ~WF_igk()" << std::endl;
	}
	delete[] g2kin;
#ifdef __CUDA
	if (this->d_g2kin)
	{
		cudaFree(this->d_g2kin);
	}
#endif
}


//--------------------------------------------------------
// Compute kinetic energy for each k-point
//--------------------------------------------------------
void WF_igk::ekin(const int ik)
{
	ModuleBase::timer::tick("WF_igk", "ekin");
	ModuleBase::GlobalFunc::ZEROS(g2kin, this->npwx);

	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		//--------------------------------------------------------
		// EXPLAIN : Put the correct units on the kinetic energy
		//--------------------------------------------------------
		this->g2kin[ig] = GlobalC::pw.get_GPlusK_cartesian(ik, this->igk(ik, ig)).norm2() * GlobalC::ucell.tpiba2;
	}
#ifdef __CUDA
	cudaMemcpy(this->d_g2kin, this->g2kin, GlobalC::kv.ngk[ik] * sizeof(double), cudaMemcpyHostToDevice);
#endif
	ModuleBase::timer::tick("WF_igk", "ekin");
	return;
}

ModuleBase::Vector3<double> WF_igk::get_1qvec_cartesian(const int ik, const int ig) const
{
	ModuleBase::Vector3<double> qvec = GlobalC::pw.get_GPlusK_cartesian(ik, this->igk(ik, ig));

	/*
	if(igk(ik,ig)==0)
	{
		std::cout << " g add = " << GlobalC::pw.gcar << std::endl;
		std::cout << "\n igk = " << igk(ik,ig);
		std::cout << " k=" << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " <<
	GlobalC::kv.kvec_c[ik].z; std::cout << " g=" << GlobalC::pw.gcar[ this->igk(ik, ig) ].x << " " << GlobalC::pw.gcar[
	this->igk(ik, ig) ].y << " " << GlobalC::pw.gcar[ this->igk(ik, ig) ].z;
	}
	*/

	return qvec;
}

double *WF_igk::get_qvec_cartesian(const int &ik)
{
	double *qmod = new double[GlobalC::kv.ngk[ik]];
	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		// cartesian coordinate
		// modulus, in GlobalC::ucell.tpiba unit.
		const double q2 = GlobalC::pw.get_GPlusK_cartesian(ik, this->igk(ik, ig)).norm2();
		qmod[ig] = GlobalC::ucell.tpiba * sqrt(q2); // sqrt(q2) * GlobalC::ucell.tpiba;
	}
	return qmod;
}

std::complex<double> *WF_igk::get_sk(const int ik, const int it, const int ia) const
{
	ModuleBase::timer::tick("WF_igk", "get_sk");
	const double arg = (GlobalC::kv.kvec_c[ik] * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
	const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));
	std::complex<double> *sk = new std::complex<double>[GlobalC::kv.ngk[ik]];
	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		const int iig = this->igk(ik, ig);
		const int iat = GlobalC::ucell.itia2iat(it, ia);
		sk[ig] = kphase * GlobalC::pw.eigts1(iat, GlobalC::pw.ig1[iig]) * GlobalC::pw.eigts2(iat, GlobalC::pw.ig2[iig])
				 * GlobalC::pw.eigts3(iat, GlobalC::pw.ig3[iig]);
	}
	ModuleBase::timer::tick("WF_igk", "get_sk");
	return sk;
}

std::complex<double> *WF_igk::get_skq(int ik,
									  const int it,
									  const int ia,
									  ModuleBase::Vector3<double> q) // pengfei 2016-11-23
{
	std::complex<double> *skq = new std::complex<double>[GlobalC::kv.ngk[ik]];

	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		ModuleBase::Vector3<double> qkq = GlobalC::pw.get_GPlusK_cartesian(ik, this->igk(ik, ig)) + q;
		double arg = (qkq * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
		skq[ig] = std::complex<double>(cos(arg), -sin(arg));
	}

	return skq;
}

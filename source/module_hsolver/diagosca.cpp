//=====================
// AUTHOR : Peize Lin
// DATE : 2021-11-02
// REFACTORING AUTHOR : Daye Zheng
// DATE : 2022-04-14
//=====================

#include "diagosca.h"
#include "module_hamilt/matrixblock.h"
#include "module_base/global_variable.h"
#include "module_base/scalapack_connector.h"
#include "module_base/global_function.h"
#include <cstring>
#include <cassert>

typedef ModuleHamilt::MatrixBlock<double> matd;
typedef ModuleHamilt::MatrixBlock<std::complex<double>> matcd;

namespace ModuleHSolver
{

void DiagoSca::diag(
    ModuleHamilt::Hamilt* phm_in,
    ModulePsi::Psi<double> &psi,
    double *eigenvalue_in)
{
    matd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
	assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);
    this->pdsygvx_diag(h_mat.desc, h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigenvalue_in, psi);
}

void DiagoSca::diag(
    ModuleHamilt::Hamilt* phm_in,
    ModulePsi::Psi<std::complex<double>> &psi,
    double *eigenvalue_in)
{
    matcd h_mat, s_mat;
    phm_in->matrix(h_mat, s_mat);
	assert(h_mat.col == s_mat.col && h_mat.row == s_mat.row && h_mat.desc == s_mat.desc);
    this->pzhegvx_diag(h_mat.desc, h_mat.col, h_mat.row, h_mat.p, s_mat.p, eigenvalue_in, psi);
}

std::pair<int,std::vector<int>> DiagoSca::pdsygvx_once(
	const int*const desc,
	const int ncol,
	const int nrow,
	const double*const h_mat,
	const double*const s_mat,
	double*const ekb,
	ModulePsi::Psi<double> &wfc_2d) const
{
	ModuleBase::matrix h_tmp(ncol, nrow, false);
	memcpy( h_tmp.c, h_mat, sizeof(double)*ncol*nrow );
	ModuleBase::matrix s_tmp(ncol, nrow, false);
	memcpy( s_tmp.c, s_mat, sizeof(double)*ncol*nrow );

	const char jobz='V', range='I', uplo='U';
	const int itype=1, il=1, iu=GlobalV::NBANDS, one=1;
	int M=0, NZ=0, lwork=-1, liwork=-1, info=0;
	const double abstol=0, orfac=-1;
	std::vector<double> work(1,0);
	std::vector<int> iwork(1,0);
	std::vector<int> ifail(GlobalV::NLOCAL,0);
	std::vector<int> iclustr(2*GlobalV::DSIZE);
	std::vector<double> gap(GlobalV::DSIZE);

	pdsygvx_(&itype, &jobz, &range, &uplo,
		&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
		NULL, NULL, &il, &iu, &abstol,
		&M, &NZ, ekb, &orfac, wfc_2d.get_pointer(), &one, &one, desc,
		work.data(), &lwork, iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
	if (info)
		throw std::runtime_error("info = "+ModuleBase::GlobalFunc::TO_STRING(info)+".\n"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

//	GlobalV::ofs_running<<"lwork="<<work[0]<<"\t"<<"liwork="<<iwork[0]<<std::endl;
	lwork = work[0];
	work.resize(lwork,0);
	liwork = iwork[0];
	iwork.resize(liwork,0);

	pdsygvx_(&itype, &jobz, &range, &uplo,
		&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
		NULL, NULL, &il, &iu, &abstol,
		&M, &NZ, ekb, &orfac, wfc_2d.get_pointer(), &one, &one, desc,
		work.data(), &lwork,
		iwork.data(), &liwork,
		ifail.data(), iclustr.data(), gap.data(), &info);
//	GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

	if(info==0)
		return std::make_pair(info, std::vector<int>{});
	else if(info<0)
		return std::make_pair(info, std::vector<int>{});
	else if(info%2)
		return std::make_pair(info, ifail);
	else if(info/2%2)
		return std::make_pair(info, iclustr);
	else if(info/4%2)
		return std::make_pair(info, std::vector<int>{M,NZ});
	else if(info/16%2)
		return std::make_pair(info, ifail);
	else
		throw std::runtime_error("info = "+ModuleBase::GlobalFunc::TO_STRING(info)+".\n"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
}



std::pair<int,std::vector<int>> DiagoSca::pzhegvx_once(
	const int*const desc,
	const int ncol,
	const int nrow,
	const std::complex<double>*const h_mat,
	const std::complex<double>*const s_mat,
	double*const ekb,
	ModulePsi::Psi<std::complex<double>> &wfc_2d) const
{
	ModuleBase::ComplexMatrix h_tmp(ncol, nrow, false);
	memcpy( h_tmp.c, h_mat, sizeof(std::complex<double>)*ncol*nrow );
	ModuleBase::ComplexMatrix s_tmp(ncol, nrow, false);
	memcpy( s_tmp.c, s_mat, sizeof(std::complex<double>)*ncol*nrow );

	const char jobz='V', range='I', uplo='U';
	const int itype=1, il=1, iu=GlobalV::NBANDS, one=1;
	int M=0, NZ=0, lwork=-1, lrwork=-1, liwork=-1, info=0;
	const double abstol=0, orfac=-1;
	std::vector<std::complex<double>> work(1,0);
	std::vector<double> rwork(1,0);
	std::vector<int> iwork(1,0);
	std::vector<int> ifail(GlobalV::NLOCAL,0);
	std::vector<int> iclustr(2*GlobalV::DSIZE);
	std::vector<double> gap(GlobalV::DSIZE);

	pzhegvx_(&itype, &jobz, &range, &uplo,
		&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
		NULL, NULL, &il, &iu, &abstol,
		&M, &NZ, ekb, &orfac, wfc_2d.get_pointer(), &one, &one, desc,
		work.data(), &lwork, rwork.data(), &lrwork,
		iwork.data(), &liwork, ifail.data(), iclustr.data(), gap.data(), &info);
	if (info)
		throw std::runtime_error("info="+ModuleBase::GlobalFunc::TO_STRING(info)+". "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));

//	GlobalV::ofs_running<<"lwork="<<work[0]<<"\t"<<"lrwork="<<rwork[0]<<"\t"<<"liwork="<<iwork[0]<<std::endl;
	lwork = work[0].real();
	work.resize(lwork,0);
	lrwork = rwork[0] + this->degeneracy_max * GlobalV::NLOCAL;
	rwork.resize(lrwork,0);
	liwork = iwork[0];
	iwork.resize(liwork,0);

	pzhegvx_(&itype, &jobz, &range, &uplo,
		&GlobalV::NLOCAL, h_tmp.c, &one, &one, desc, s_tmp.c, &one, &one, desc,
		NULL, NULL, &il, &iu, &abstol,
		&M, &NZ, ekb, &orfac, wfc_2d.get_pointer(), &one, &one, desc,
		work.data(), &lwork,
		rwork.data(), &lrwork,
		iwork.data(), &liwork,
		ifail.data(), iclustr.data(), gap.data(), &info);
//	GlobalV::ofs_running<<"M="<<M<<"\t"<<"NZ="<<NZ<<std::endl;

	if(info==0)
		return std::make_pair(info, std::vector<int>{});
	else if(info<0)
		return std::make_pair(info, std::vector<int>{});
	else if(info%2)
		return std::make_pair(info, ifail);
	else if(info/2%2)
		return std::make_pair(info, iclustr);
	else if(info/4%2)
		return std::make_pair(info, std::vector<int>{M,NZ});
	else if(info/16%2)
		return std::make_pair(info, ifail);
	else
		throw std::runtime_error("info = "+ModuleBase::GlobalFunc::TO_STRING(info)+".\n"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
}

void DiagoSca::pdsygvx_diag(
	const int*const desc,
	const int ncol,
	const int nrow,
	const double*const h_mat,
	const double*const s_mat,
	double*const ekb,
	ModulePsi::Psi<double> &wfc_2d)
{
	while(true)
	{
		const std::pair<int, std::vector<int>> info_vec = pdsygvx_once(desc, ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
		post_processing(info_vec.first, info_vec.second);
		if(info_vec.first==0)
			break;
	}
}



void DiagoSca::pzhegvx_diag(
	const int*const desc,
	const int ncol,
	const int nrow,
	const std::complex<double>*const h_mat,
	const std::complex<double>*const s_mat,
	double*const ekb,
	ModulePsi::Psi<std::complex<double>> &wfc_2d)
{
	while(true)
	{
		const std::pair<int, std::vector<int>> info_vec = pzhegvx_once(desc, ncol, nrow, h_mat, s_mat, ekb, wfc_2d);
		post_processing(info_vec.first, info_vec.second);
		if(info_vec.first==0)
			break;
	}
}



void DiagoSca::post_processing(const int info, const std::vector<int> &vec)
{
	const std::string str_info = "info = "+ModuleBase::GlobalFunc::TO_STRING(info)+".\n";
	const std::string str_FILE = ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)+".\n";
	const std::string str_info_FILE = str_info + str_FILE;

	if(info==0)
	{
		return;
	}
	else if(info<0)
	{
		const int info_negative = -info;
		const std::string str_index = (info_negative>100)
			? ModuleBase::GlobalFunc::TO_STRING(info_negative/100)+"-th argument "+ModuleBase::GlobalFunc::TO_STRING(info_negative%100)+"-entry is illegal.\n"
			: ModuleBase::GlobalFunc::TO_STRING(info_negative)+"-th argument is illegal.\n";
		throw std::runtime_error(str_info_FILE + str_index);
	}
	else if(info%2)
	{
		std::string str_ifail = "ifail = ";
		for(const int i : vec)
			str_ifail += ModuleBase::GlobalFunc::TO_STRING(i) + " ";
		throw std::runtime_error(str_info_FILE + str_ifail);
	}
	else if(info/2%2)
	{
		int degeneracy_need = 0;
		for(int irank=0; irank<GlobalV::DSIZE; ++irank)
			degeneracy_need = std::max(degeneracy_need, vec[2*irank+1]-vec[2*irank]);
		const std::string str_need = "degeneracy_need = "+ModuleBase::GlobalFunc::TO_STRING(degeneracy_need)+".\n";
		const std::string str_saved = "degeneracy_saved = "+ModuleBase::GlobalFunc::TO_STRING(this->degeneracy_max)+".\n";
		if(degeneracy_need <= this->degeneracy_max)
		{
			throw std::runtime_error(str_info_FILE + str_need + str_saved);
		}
		else
		{
			GlobalV::ofs_running<<str_need<<str_saved;
			this->degeneracy_max = degeneracy_need;
			return;
		}
	}
	else if(info/4%2)
	{
		const std::string str_M = "M = "+ModuleBase::GlobalFunc::TO_STRING(vec[0])+".\n";
		const std::string str_NZ = "NZ = "+ModuleBase::GlobalFunc::TO_STRING(vec[1])+".\n";
		const std::string str_NBANDS = "GlobalV::NBANDS = "+ModuleBase::GlobalFunc::TO_STRING(GlobalV::NBANDS)+".\n";
		throw std::runtime_error(str_info_FILE + str_M + str_NZ + str_NBANDS);
	}
	else if(info/16%2)
	{
		const std::string str_npos = "not positive definite = "+ModuleBase::GlobalFunc::TO_STRING(vec[0])+".\n";
		throw std::runtime_error(str_info_FILE + str_npos);
	}
	else
	{
		throw std::runtime_error(str_info_FILE);
	}
}

}
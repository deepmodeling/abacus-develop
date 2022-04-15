#include "hamiltlcao.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

namespace ModuleHamilt
{

template<typename T, typename T1>
void HamiltLCAO<T,T1>::constructFixedReal()
{
    //update the flag
    if(HamiltLCAO<T,T1>::isFixedDone)return;
    else HamiltLCAO<T,T1>::isFixedDone=true;

    if(GlobalV::GAMMA_ONLY_LOCAL) this->fixedRealM.resize(this->LM->ParaV->nloc, T1(0));
    else this->fixedRealM.resize(this->LM->ParaV->nnr, T1(0));

    ModuleBase::TITLE("HamiltLCAO","constructFixedReal"); 
    ModuleBase::timer::tick("HamiltLCAO","constructFixedReal"); 

	// (9) compute S, T, Vnl, Vna matrix.
    this->genH->calculate_S_no(this->fixedRealM.getS());

    if(GlobalV::VNL_IN_H)
    {
        this->genH->calculate_NL_no(this->fixedRealM.getH());
    }

    if(GlobalV::T_IN_H)
    {
        this->genH->calculate_T_no(this->fixedRealM.getH());
    }

    ModuleBase::timer::tick("HamiltLCAO","constructFixedReal"); 
	return;

    return;
}

template<typename T, typename T1>
void HamiltLCAO<T,T1>::constructUpdateReal()
{
    return;
}

template<typename T, typename T1>
void HamiltLCAO<T,T1>::hk_fixed_mock(const int ik)
{
    if(GlobalV::GAMMA_ONLY_LOCAL) this->kM.resize(this->LM->ParaV->nloc, T(0));
    else this->kM.resize(this->LM->ParaV->nnr, T(0));

	// folding_fixedH() should be refactored to there, 
    // but now deepks code in this function, should be moved to another place
	return;
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::hk_update_mock(const int ik)
{
    //update_Hk and update_Hgamma should be refactored to there
    return;
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::getMatrix(MatrixBlock<T> &hk_in, MatrixBlock<T> &sk_in)
{
    hk_in = MatrixBlock<T>{kM.hloc.data(), this->LM->ParaV->nrow, this->LM->ParaV->ncol, this->LM->ParaV->desc};
    sk_in = MatrixBlock<T>{kM.sloc.data(), this->LM->ParaV->nrow, this->LM->ParaV->ncol, this->LM->ParaV->desc};
}


// case for nspin==4
/*void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}*/
// case for nspin<4, multi-k-points 
void HamiltLCAO<std::complex<double>, double>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

// case for nspin<4, gamma_only
void HamiltLCAO<double, double>::matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)
{
    this->getMatrix(hk_in, sk_in);
}

// nspin==4 case not supported yet
/*void HamiltLCAO<std::complex<double>, std::complex<double>>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}*/
void HamiltLCAO<std::complex<double>, double>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}
void HamiltLCAO<double, double>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}

/*void HamiltLCAO<std::complex<double>, std::complex<double>>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};*/
void HamiltLCAO<double, std::complex<double>>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};
void HamiltLCAO<double, double>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};

}
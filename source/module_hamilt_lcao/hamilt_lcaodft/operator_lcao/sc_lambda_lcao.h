#ifndef SC_LAMBDA_LCAO_H
#define SC_LAMBDA_LCAO_H

#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __OPLAMBDATEMPLATE
#define __OPLAMBDATEMPLATE

template <class T>
class OperatorScLambda : public T
{
};

#endif

template<typename T>
class OperatorScLambda<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    OperatorScLambda<OperatorLCAO<T>>(LCAO_Matrix* LM_in,
                                    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                    std::vector<double>* HR_pointer_in,
                                    std::vector<T>* HK_pointer_in)
        : HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_sc_lambda;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

    void cal_h_lambda(int ik, T* h_lambda);

    void cal_weight_func(const std::vector<T>& Sloc2);

    // setters
    void set_nat(int nat_in);

    void set_nloc(int nloc_in);

    void set_lambda(const std::vector<ModuleBase::Vector3<double>>& lambda_in);

    void set_iwt2iat(const int* iwt2iat_in);

    void set_Wi(const std::vector<T>& Wi_in);

    // getters
    int get_nat();

    int get_nloc();

    const std::vector<T>& get_Wi() const;

    const std::vector<ModuleBase::Vector3<double>>& get_lambda() const;

    const std::vector<int>& get_iwt2iat() const;

  private:
    std::vector<double>* HR_pointer = nullptr;
    std::vector<T>* HK_pointer = nullptr;
    std::vector<ModuleBase::Vector3<double>> lambda_;
    std::vector<int> iwt2iat_;
    std::vector<T> Wi_;
    int nloc_;
    int nat_;
};

} // namespace hamilt

#endif // SC_LAMBDA_LCAO_H
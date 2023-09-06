#include "spin_constrain.h"
#include "module_base/timer.h"
#include "module_base/tool_title.h"
#include "module_base/global_function.h"
#include <algorithm>

void SpinConstrain::cal_h_lambda(const std::vector<std::complex<double>>& Sloc2, std::complex<double>* h_lambda)
{
    ModuleBase::TITLE("SpinConstrain","cal_h_lambda");
    ModuleBase::timer::tick("SpinConstrain", "cal_h_lambda");
    const Parallel_Orbitals* pv = this->ParaV;
    for (const auto& sc_elem1 : this->get_atomCounts())
    {
        int it1 = sc_elem1.first;
        int nat_it1 = sc_elem1.second;
        int nw_it1 = this->get_orbitalCounts().at(it1);
        for (int ia1 = 0; ia1 < nat_it1; ia1++)
        {
            int iat1 = this->get_iat(it1, ia1);
            for (int iw1 = 0; iw1 < nw_it1; iw1++)
            {
                int iwt1 = this->get_iwt(it1, ia1, iw1);
                const int mu = pv->global2local_row(iwt1);
                if (mu < 0) continue;
                for (const auto& sc_elem2 : this->get_atomCounts())
                {
                    int it2 = sc_elem2.first;
                    int nat_it2 = sc_elem2.second;
                    int nw_it2 = this->get_orbitalCounts().at(it2);
                    for (int ia2 = 0; ia2 < nat_it2; ia2++)
                    {
                        int iat2 = this->get_iat(it2, ia2);
                        for (int iw2 = 0; iw2 < nw_it2; iw2++)
                        {
                            int iwt2 = this->get_iwt(it2, ia2, iw2);
                            const int nu = pv->global2local_col(iwt2);
                            if (nu < 0) continue;
                            int icc;
                            if (ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER())
						    {
                                icc = mu + nu * pv->nrow;
                            }
                            else
                            {
                                icc = mu * pv->ncol + nu;
                            }
                            this->Wi_[icc] = (iat1 == iat2) ? Sloc2[icc] : Sloc2[icc]*0.5;
                            if (iwt1 % 2 == 0)
                            {
                                h_lambda[icc] = (iwt2 % 2 == 0) ?
                                    Wi_[icc] * lambda_[iat1][2] :
                                    Wi_[icc] * (lambda_[iat1][0] + lambda_[iat1][1] * std::complex<double>(0, -1));
                            }
                            else
                            {
                                h_lambda[icc] = (iwt2 % 2 == 0) ?
                                    Wi_[icc] * (lambda_[iat1][0] + lambda_[iat1][1] * std::complex<double>(0, 1)) :
                                    Wi_[icc] * (-lambda_[iat1][2]);
                            }
                        }
                    }

                }
            }
        }
    }
    ModuleBase::timer::tick("SpinConstrain", "cal_h_lambda");
    return;
}
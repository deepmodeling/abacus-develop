#ifndef LANMA_LOOP_H
#define LANMA_LOOP_H

#include <vector>

void lambda_loop(const std::vector<std::vector<double>>& M_CONSTR,
                 const std::vector<std::vector<int>>& CONSTRL,
                 const int NIONS,
                 const int NTYP,
                 const std::vector<int>& NITYP,
                 const double INISC,
                 const double SCDIFF,
                 const std::vector<double>& SCCONV_GRAD,
                 const int NSC,
                 const int NSCMIN,
                 const double SCCUT,
                 const int N,
                 std::vector<std::vector<double>>& MW,
                 std::vector<std::vector<double>>& OUT_LAMBDA);
//  TODO datatype
//<> CHTOT,
//<> CHTOTL,
//<> W)

#endif // LANMA_FIELD_H
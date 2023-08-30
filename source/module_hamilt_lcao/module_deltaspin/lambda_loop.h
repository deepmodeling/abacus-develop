#ifndef LAMBDA_LOOP_H
#define LAMBDA_LOOP_H

#include <iostream>
#include <vector>
#include "module_base/vector3.h"

class LambdaLoop
{
    public:
         LambdaLoop();
        ~LambdaLoop();

        void init_input_parameters();

        void run_lambda_loop(int outer_step);

        std::vector<ModuleBase::Vector3<double>> target_spin; // which is M_CONSTR from INPUT
        std::vector<ModuleBase::Vector3<int>> constrain; // which is CONSTRL from INPUT
        int nat; // NIONS changed to nat
        int ntype; // NTYP changed to ntype
        std::vector<int> iat2it; // NITYP changed to iat2it
        double alpha_trial; // which is INISC from INPUT
        double epsilon; // which is SCDIFF from INPUT
        std::vector<double> bound_gradient; // which is SCCONV_GRAD from INPUT
        int num_step; // which is NSC from INPUT
        int num_min; // which is NSCMIN from INPUT
        double restrict_current; // which is SCCUT from INPUT
        int outer_step; // which is N in the original code
        std::vector<ModuleBase::Vector3<double>> spin; // which is MW from INPUT, the initial spin
        std::vector<ModuleBase::Vector3<double>> out_lambda; // which is OUT_LAMBDA from INPUT
        bool debug = false; // which is LDESC in the original code
};

#endif // LAMBDA_LOOP_H
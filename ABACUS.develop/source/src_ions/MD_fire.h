#ifndef MD_FIRE_H
#define MD_FIRE_H

#include "src_pw/tools.h"

class MD_fire
{
    public:
        MD_fire();
        ~MD_fire(){};
        void check_FIRE(
            const int& numIon, 
            const Vector3<double>* force,
            double& deltaT, 
            Vector3<double>* vel);
    private:
        double alpha_start ;//alpha_start begin
        double alpha;//alpha begin
        double finc;//finc begin
        double fdec;//fdec begin
        double f_alpha;
        int N_min;
        double dt_max;//dt_max
        int negative_count;//Negative count
};

#endif

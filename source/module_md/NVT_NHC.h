#ifndef NVT_NHC_H
#define NVT_NHC_H

#include "mdrun.h"

class NVT_NHC : public MDrun
{
public:
    NVT_NHC(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~NVT_NHC();

    void setup(ModuleESolver::ESolver *p_ensolve);
    void first_half();
    void second_half();
    void outputMD(std::ofstream &ofs, bool cal_stress);
    void write_restart();
    void restart();
    void integrate();
    void update_mass();

    double t_target;       // target temperature
    const static int nc = 1;
    const static int nsy = 7;
    double w[nsy];         // scale evolution operator

    double *Q;             // MHC mass
    double *G;             // NHC acceleration
    double *eta;           // NHC position
    double *veta;          // NHC velocity

};

#endif
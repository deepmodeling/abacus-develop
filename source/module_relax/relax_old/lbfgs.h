#ifndef LBFGS_H
#define LBFGS_H

/**
 * @file bfgs.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2024-11-28
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <vector>
#include <tuple> 
#include<algorithm>
#include<cmath>
#include"module_base/lapack_connector.h"

#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell.h"
#include "module_esolver/esolver.h"
#include "module_esolver/esolver_ks.h"


class LBFGS
{
public:
    
    double alpha;//initialize H,diagonal element is alpha
    double maxstep;//every movement smaller than maxstep
    int size;//number of etoms
    int memory;
    double damping;
    double H0;
    int iteration;
    double energy;
    double alpha_k;
    double xtol;
    double stpmin;
    double stpmax;
    bool bracket;
    double xtrapl;
    double xtrapu;
    double old_stp;
    bool no_update;
    int task;
    ModuleESolver::ESolver* solver;

    
    std::vector<std::vector<double>> H;
    std::vector<double> force0;
    std::vector<std::vector<double>> force;
    std::vector<double> pos0;
    std::vector<std::vector<double>> pos;
    std::vector<double> pos_taud0;
    std::vector<std::vector<double>> pos_taud;
    std::vector<std::vector<double>> dpos;
    std::vector<std::vector<double>> s;
    std::vector<std::vector<double>> y;
    std::vector<double> rho;
    std::vector<int> isave;
    std::vector<double> dsave;
    std::vector<double> steplength;

    /**
     * @brief 
     * 
     * @param _size 
     */
    void allocate(const int _size);//initialize parameters
    void relax_step(const ModuleBase::matrix _force,UnitCell& ucell,const double &etot,ModuleESolver::ESolver* p_esolver);//
    void PrepareStep(std::vector<std::vector<double>>& force,std::vector<std::vector<double>>& pos,std::vector<std::vector<double>>& H,std::vector<double>& pos0,std::vector<double>& force0,std::vector<std::vector<double>>& dpos,UnitCell& ucell,const double &etot);
    void IsRestrain(std::vector<std::vector<double>>& dpos);
    

private:
    bool sign;
    
    void CalculateLargestGrad(const ModuleBase::matrix& _force,UnitCell& ucell);
    void GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos);
    void GetPostaud(UnitCell& ucell,std::vector<std::vector<double>>& pos_taud);

    void Update(std::vector<std::vector<double>>& pos_taud, 
                std::vector<double>& pos_taud0, 
                std::vector<double>& force,
                std::vector<double>& force0, 
                UnitCell& ucell,
                int iteration,
                int memory,
                std::vector<std::vector<double>>& s,
                std::vector<std::vector<double>>& y,
                std::vector<double>& rho);
    double LineSearch(UnitCell& ucell,std::vector<std::vector<double>>& r,std::vector<std::vector<double>>& g,double e);
    double Step(double stp,double f,double g,double c1,double c2,double xtol,std::vector<int>& isave,std::vector<double>& dsave);
    void UpdateLineSearch(double& stx,double& fx,double& gx,double& sty,double& fy,double& gy,double& stp,double& fp,double& gp,double& stpmin,double& stpmax);
    double DetermineStep(double stp);
    void DetermineStep(std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,double& maxstep);
    void UpdatePos(UnitCell& ucell);
    void Save(int a,double b,double c,double d,double e,double f,double g,double h,double i,double j,double k,double l,double m,double n);
    double GetEnergy(UnitCell& ucell,double stp);
    std::vector<double> GetForce(UnitCell& ucell,double stp);
    

    // matrix method
    std::vector<double> ReshapeMToV(std::vector<std::vector<double>>& matrix);
    std::vector<std::vector<double>> MAddM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b);
    std::vector<double> VSubV(std::vector<double>& a, std::vector<double>& b);
    std::vector<double> VAddV(std::vector<double>& a, std::vector<double>& b);
    std::vector<std::vector<double>> ReshapeVToM(std::vector<double>& matrix);
    std::vector<double> DotInMAndV1(std::vector<std::vector<double>>& matrix, std::vector<double>& vec);
    std::vector<double> DotInMAndV2(std::vector<std::vector<double>>& matrix, std::vector<double>& vec);
    double DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2);
    std::vector<std::vector<double>> OuterVAndV(std::vector<double>& a, std::vector<double>& b);
    std::vector<std::vector<double>> MPlus(std::vector<std::vector<double>>& a, double b);
    std::vector<std::vector<double>> MSubM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b);
    std::vector<double> DotInVAndFloat(std::vector<double>& vec, double b); 
};

#endif // BFGS_H

#ifndef BFGS_H
#define BFGS_H

/**
 * @file bfgs.h
 * @author 19hello (you@domain.com)
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



class BFGS
{
public:
    
    double alpha;//initialize H,diagonal element is alpha
    double maxstep;//every movement smaller than maxstep
    int size;//number of atoms

    
    std::vector<double> steplength;//the length of atoms displacement 
    std::vector<std::vector<double>> H;//Hessian matrix
    std::vector<double> force0;//force in previous step
    std::vector<std::vector<double>> force;
    std::vector<double> pos0;//atom pos in previous step(cartesian coordinates)
    std::vector<std::vector<double>> pos;
    std::vector<double> pos_taud0;//atom pos in previous step(relative coordinates)
    std::vector<std::vector<double>> pos_taud;
    std::vector<std::vector<double>> dpos;

    /**
     * @brief 
     * 
     * @param _size 
     */
    void allocate(const int _size);//initialize parameters
    void relax_step(const ModuleBase::matrix& _force,UnitCell& ucell);//a full iteration step
    void PrepareStep(std::vector<std::vector<double>>& force,std::vector<std::vector<double>>& pos,std::vector<std::vector<double>>& H,std::vector<double>& pos0,std::vector<double>& force0,std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,UnitCell& ucell);//calculate the atomic displacement in one iteration step
    void IsRestrain(std::vector<std::vector<double>>& dpos);//check if converged

private:
    bool sign;//check if this is the first iteration
    
    void CalculateLargestGrad(const ModuleBase::matrix& _force,UnitCell& ucell);
    void GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos);
    void GetPostaud(UnitCell& ucell,std::vector<std::vector<double>>& pos_taud);
    void Update(std::vector<double>& pos, std::vector<double>& force,std::vector<std::vector<double>>& H,UnitCell& ucell);//update hessian matrix
    void DetermineStep(std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,double& maxstep);//normalize large atomic displacements based on maxstep
    void UpdatePos(UnitCell& ucell);//update ucell with the new coordinates 
    

    // matrix method
    std::vector<double> ReshapeMToV(std::vector<std::vector<double>>& matrix);
    std::vector<std::vector<double>> MAddM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b);
    std::vector<double> VSubV(std::vector<double>& a, std::vector<double>& b);
    std::vector<std::vector<double>> ReshapeVToM(std::vector<double>& matrix);
    std::vector<double> DotInMAndV1(std::vector<std::vector<double>>& matrix, std::vector<double>& vec);
    std::vector<double> DotInMAndV2(std::vector<std::vector<double>>& matrix, std::vector<double>& vec);
    double DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2);
    std::vector<std::vector<double>> OuterVAndV(std::vector<double>& a, std::vector<double>& b);
    std::vector<std::vector<double>> MPlus(std::vector<std::vector<double>>& a, double& b);
    std::vector<std::vector<double>> MSubM(std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& b);
};

#endif // BFGS_H

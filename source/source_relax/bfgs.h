#ifndef BFGS_H
#define BFGS_H

#include <vector>
#include <tuple> 
#include<algorithm>
#include<cmath>
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_cell/unitcell.h"
#include "matrix_methods.h"

class BFGS
{
public:
    std::vector<double> steplength;//the length of atoms displacement 
    std::vector<std::vector<double>> H;//Hessian matrix
    std::vector<double> force0;//force in previous step
    std::vector<std::vector<double>> force;
    std::vector<double> pos0;//atom pos in previous step(cartesian coordinates)
    std::vector<std::vector<double>> pos;
    std::vector<double> pos_taud0;//atom pos in previous step(relative coordinates)
    std::vector<std::vector<double>> pos_taud;
    std::vector<std::vector<double>> dpos;

    void allocate(const int _size);//initialize parameters
    void relax_step(const ModuleBase::matrix& _force,UnitCell& ucell);//a full iteration step
    void PrepareStep(std::vector<std::vector<double>>& force,std::vector<std::vector<double>>& pos,std::vector<std::vector<double>>& H,std::vector<double>& pos0,std::vector<double>& force0,std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,UnitCell& ucell);//calculate the atomic displacement in one iteration step

private:
    bool sign;//check if this is the first iteration
    double alpha;//initialize H,diagonal element is alpha
    double maxstep;//every movement smaller than maxstep
    int size;//number of atoms
    
    void IsRestrain(std::vector<std::vector<double>>& dpos);//check if converged
    void CalculateLargestGrad(const ModuleBase::matrix& _force,UnitCell& ucell);
    void GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos);
    void GetPostaud(UnitCell& ucell,std::vector<std::vector<double>>& pos_taud);
    void Update(std::vector<double>& pos, std::vector<double>& force,std::vector<std::vector<double>>& H,UnitCell& ucell);//update hessian matrix
    void DetermineStep(std::vector<double>& steplength,std::vector<std::vector<double>>& dpos,double& maxstep);//normalize large atomic displacements based on maxstep
    void UpdatePos(UnitCell& ucell);//update ucell with the new coordinates
    
};

#endif // BFGS_H

#ifndef MATRIX_METHODS
#define MATRIX_METHODS

#include <vector>
#include <cassert>
#include "source_base/vector3.h"
#include "source_base/matrix.h"


std::vector<double> ReshapeMToV(std::vector<ModuleBase::Vector3<double>>& matrix);
ModuleBase::matrix MAddM(ModuleBase::matrix& a, ModuleBase::matrix& b);
std::vector<double> VSubV(std::vector<double>& a, std::vector<double>& b);
std::vector<double> VAddV(std::vector<double>& a, std::vector<double>& b);
std::vector<ModuleBase::Vector3<double>> ReshapeVToM(std::vector<double>& matrix);
std::vector<double> DotInMAndV1(ModuleBase::matrix& matrix, std::vector<double>& vec);
std::vector<double> DotInMAndV2(ModuleBase::matrix& matrix, std::vector<double>& vec);
double DotInVAndV(std::vector<double>& vec1, std::vector<double>& vec2);
ModuleBase::matrix OuterVAndV(std::vector<double>& a, std::vector<double>& b);
ModuleBase::matrix MPlus(ModuleBase::matrix& a, double b);
ModuleBase::matrix MSubM(ModuleBase::matrix& a, ModuleBase::matrix& b);
std::vector<double> DotInVAndFloat(std::vector<double>& vec, double b); 



#endif
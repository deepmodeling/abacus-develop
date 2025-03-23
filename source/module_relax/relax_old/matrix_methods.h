#ifndef MATRIX_METHODS
#define MATRIX_METHODS

#include <vector>



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



#endif
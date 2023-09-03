#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <vector>
#include <map>
#include "module_base/vector3.h"

struct ScAtomData;

class SpinConstrain
{
public:
    // Public method to access the Singleton instance
    static SpinConstrain& getInstance();

    // Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;             // Copy construct
    SpinConstrain(SpinConstrain&&) = delete;                  // Move construct

    /**
    * @brief parse json input file for non-collinear spin-constrained DFT
    */
    void parseScJsonFile(const std::string& filename, std::map<std::string, std::vector<ScAtomData>>& data);

    void set_nat(int nat);
    int get_nat();

private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain(){};                              // Destructor
    SpinConstrain& operator=(SpinConstrain const&) = delete;  // Copy assign
    SpinConstrain& operator=(SpinConstrain &&) = delete;      // Move assign
    std::vector<ModuleBase::Vector3<double>> lambda;
    std::vector<ModuleBase::Vector3<double>> sc_mag;
    std::vector<double> sc_spin_val;
    std::vector<double> sc_spin_angle1;
    std::vector<double> sc_spin_angle2;
    int nat;
};


/**
 * @brief struct for storing parameters of non-collinear spin-constrained DFT
 */
struct ScAtomData {
    int index;
    std::vector<double> lambda;
    std::vector<double> sc_mag;
    double sc_spin_val;
    double sc_spin_angle1;
    double sc_spin_angle2;
};

#endif // SPIN_CONSTRAIN_H
#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <vector>
#include <map>
#include "module_base/vector3.h"

struct ScAtomData;

class SpinConstrain
{
public:
    /// Public method to access the Singleton instance
    static SpinConstrain& getInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// parse json input file for non-collinear spin-constrained DFT
    void Set_ScData_From_Json(const std::string& filename);
    /// get sc_data
    std::map<int, std::vector<ScAtomData>>& get_ScData();
    /// clear sc_data
    void clear_ScData();
    /// set element index to atom index map
    void set_itia(const std::map<int, int>& itia);
    /// get element index to atom index map
    std::map<int, int>& get_itia();
    /// get sc_lambda
    std::vector<ModuleBase::Vector3<double>> get_sc_lambda();
    /// get sc_mag
    std::vector<ModuleBase::Vector3<double>> get_sc_mag();

private:
    SpinConstrain(){};                               // Private constructor
    ~SpinConstrain(){};                              // Destructor
    SpinConstrain& operator=(SpinConstrain const&) = delete;  // Copy assign
    SpinConstrain& operator=(SpinConstrain &&) = delete;      // Move assign
    std::map<int, std::vector<ScAtomData>> ScData;
    std::map<int, int> itia;
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
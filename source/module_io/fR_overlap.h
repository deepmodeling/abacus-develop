#ifndef FR_OVERLAP_H
#define FR_OVERLAP_H
#ifdef __LCAO
#include <complex>
#include <functional>
#include "module_basis/module_ao/parallel_orbitals.h"
#include "module_cell/module_neighbor/sltk_grid_driver.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_lcao/module_hcontainer/hcontainer.h"
#include "module_base/math_lebedev_laikov.h"


template <typename T>
class FR_overlap
{
public:
    using fr_ptr = std::function<T(ModuleBase::Vector3<double>)>;

    FR_overlap();

    void set_parameters(
        fr_ptr fr_in, 
        const UnitCell* ucell_in, 
        Grid_Driver* GridD_in, 
        const Parallel_Orbitals* paraV,
        int radial_grid_num = 140,
        int degree = 110
    );

    FR_overlap(const FR_overlap<T>& FR_in);

    FR_overlap(FR_overlap<T>&& FR_in);

    ~FR_overlap();

    void calculate_FR();

    hamilt::HContainer<T>* get_FR_pointer() const
    {
        return this->FR_container;
    }

protected:
    void initialize_FR(Grid_Driver* GridD, const Parallel_Orbitals* paraV);

    void cal_FR_IJR(const int& iat1, const int& iat2, const Parallel_Orbitals* paraV, const ModuleBase::Vector3<double>& dtau, T* data_pointer);

    std::map<std::pair<int, int>, double> psi_inter(const int &T1, const std::set<std::pair<int, int>> &LN_pair1, const double &r_norm);

    double Polynomial_Interpolation(
        const double *psi_r,
        const int &mesh_r,
        const double &dr,
        const double &x	
    );

    fr_ptr fr = nullptr;
    const UnitCell* ucell = nullptr;
    int radial_grid_num = 140;
    ModuleBase::Lebedev_laikov_grid *Leb_grid = nullptr;
    hamilt::HContainer<T> *FR_container = nullptr;
};
#endif
#endif
#include "source_cell/module_neighbor_search/unitcell_interface.h"
#include <vector>

class UnitCellPlus : public IAtomProvider
{
public:
    UnitCellPlus()=default;
    ~UnitCellPlus()=default;

    double get_lat0() const override {
        return lat0;
    }

    double get_omega() const override {
        return omega;
    }

    const ModuleBase::Matrix3& get_latvec() const override {
        return latvec;
    }

    int get_natom() const override {
        return nat;
    }

    int get_na(int i) const override {
        return na[i];
    }

    int get_ntype() const override {
        return ntype;
    }

    ModuleBase::Vector3<double> get_tauu(int i, int j) const override {
        int count=0;
        for(int k=0;k<i;k++)
        {
            count += na[k];
        }
        count+=j;
        return tau[count];
    }

    double lat0;
    double omega;
    int nat;
    std::vector<int> na;
    int ntype;
    ModuleBase::Matrix3 latvec;
    std::vector<ModuleBase::Vector3<double>> tau;

};
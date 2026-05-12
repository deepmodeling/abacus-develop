#pragma once

#include "source_cell/module_neighlist/unitcell_interface.h"
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
        assert(i >= 0 && i < ntype);
        return na[i];
    }

    int get_ntype() const override {
        return ntype;
    }

    ModuleBase::Vector3<double> get_tauu(int i, int j) const override {
        assert(i >= 0 && i < ntype);
        assert(j >= 0 && j < na[i]);
        //assert(naa.size() > 0);
        if(i==0)
        {
            return tau[j];
        }
        return tau[naa[i-1]+j];
    }

    double lat0;
    double omega;
    int nat;
    std::vector<int> na;
    // naa is the cumulative sum of na: naa[i] = na[0] + na[1] + ... + na[i]
    std::vector<int> naa;
    int ntype;
    ModuleBase::Matrix3 latvec;
    std::vector<ModuleBase::Vector3<double>> tau;

    // Compute cumulative counts in naa from current na vector.
    // Assumption: naa[i] == sum_{k=0..i} na[k]. Call this after na is set/modified.
    void compute_naa()
    {
        naa.resize(na.size());
        if(naa.size()>0)
        {
            naa[0] = na[0];
        }
        for(int i=1;i<naa.size();i++)
        {
            naa[i] = naa[i-1] + na[i];
        }
    }

};
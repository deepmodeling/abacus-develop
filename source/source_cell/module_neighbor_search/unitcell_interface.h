#ifndef UNITCELL_INTERFACE_H
#define UNITCELL_INTERFACE_H


#include "source_base/vector3.h"
#include "source_base/matrix3.h"


class IAtomProvider {
public:
    virtual ~IAtomProvider() = default;

    virtual double get_lat0() const = 0;
    virtual double get_omega() const = 0;
    virtual const ModuleBase::Matrix3& get_latvec() const = 0;


 
    virtual int get_natom() const = 0;
    virtual int get_na(int i) const = 0;
    virtual int get_ntype() const = 0;
    virtual ModuleBase::Vector3<double> get_tauu(int i,int j) const = 0;
};

#endif // UNITCELL_INTERFACE_H
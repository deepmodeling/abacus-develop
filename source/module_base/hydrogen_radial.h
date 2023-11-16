#include "associate_lagurre.h"

template<typename T>
class Hydrogen_Radial
{
    public:
        Hydrogen_Radial(){};
        ~Hydrogen_Radial(){};

        static T value(const int& Z, const int& n, const int& l, const T& r);
        
        static void generate(const int& Z, const int& n, const int& l, T* r, T* Rnl, const int& size);
    private:
        T Z = 1; // charge
        int l = 0; // angular momentum

        const T me = 9.10938356e-31; // electron mass, kg
        const T e = 1.60217662e-19; // electron charge, C
        const T hbar = 1.0545718e-34; // reduced Planck constant, J*s
        Associate_Laguerre<T> L;
};
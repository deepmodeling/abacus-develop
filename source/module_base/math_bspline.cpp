#include "math_bspline.h"

#include <assert.h>

#include <cstdlib>

namespace ModuleBase
{
    Bspline::Bspline()
    {
        bezier = nullptr;
        norder = 0;
        xi = 0;
        Dx = 1.0; 
    }
    Bspline::~Bspline()
    {
        delete[] bezier;
    }

    void Bspline::init(int norderin, double Dxin, double xiin)
    {
        this->xi = xiin;
        this->Dx = Dxin;
        this->norder = norderin;
        assert(Dx > 0);
        //norder must be a positive even number.
        assert(norder > 0);
        assert(norder % 2 == 0); 
        delete[] bezier; bezier = new double [this->norder+1];
        for(int i = 0 ; i < norder+1 ; ++i)
        {
            bezier[i] = 0;
        }
    }

    double Bspline::bezier_ele(int n)
    {
        return this->bezier[n];
    }

    void Bspline::getbspline(double x)
    {
        bezier[0] = 1.0;
        for(int k = 1 ; k <= norder ; ++k)
        {
            //for n>=1
            for(int n = k; n >= 1; --n )
            {
                this->bezier[n] = ((x + n*this->Dx - this->xi)*this->bezier[n] + 
                (this->xi + (k-n+1)*Dx - x)*this->bezier[n-1])/(k*this->Dx);
            }

            //for n = 0
            this->bezier[0] = (x - this->xi)*this->bezier[0] / (k*this->Dx);
        }
    }
}



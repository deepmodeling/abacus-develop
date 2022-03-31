#ifndef IONS_MOVE_CG_H
#define IONS_MOVE_CG_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
class Ions_Move_CG
{
public:
    Ions_Move_CG();
    ~Ions_Move_CG();

    void allocate(void);
    void start(const ModuleBase::matrix& force, const double &etot);
        
        static double RELAX_CG_THR_E;
	int sd_step;
	int cg_step;

private:
	double* pos0;
	double* grad0;
        double* cg_grad0;
	double* move0;
	double  e0;
        // setup gradients.
	void setup_cg_grad(double *grad, const double *grad0, double *cg_grad, const double *cg_grad0, const int &ncggrad, int &flag); //LiuXh fix bug of lpf, 20180515
	void setup_move( double *move, double *cg_gradn, const double &trust_radius );
	void setup_etot_cg(const double &energy_in, const bool sd, const bool trial);
        void Brent(double &fa,double &fb,double &fc,double &xa,double &xb,double &xc,double &best_x,double &xpt);
        void f_cal(const double *g0,const double *g1,const int &dim,double &f_value);
        void third_order(const double &e0,const double &e1,const double &fa, const double &fb, const double x, double &best_x);
};


#endif

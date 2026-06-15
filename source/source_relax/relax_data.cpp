#include "relax_data.h"

int Relax_Data::dim = 0;
int Relax_Data::dim_lattice = 0;

double Relax_Data::largest_grad = 0.0;

void Relax_Data::allocate(const int dim_in) {
    pos.resize(dim_in, 0.0);
    grad.resize(dim_in, 0.0);
    move.resize(dim_in, 0.0);
    pos_p.resize(dim_in, 0.0);
    grad_p.resize(dim_in, 0.0);
    move_p.resize(dim_in, 0.0);
}
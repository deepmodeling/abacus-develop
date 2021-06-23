#include "gint_gamma.h"
#include "src_pw/global.h"
#include "src_global/ylm.h"
#include "module_neighbor/sltk_atom_arrange.h"

Gint_Gamma::Gint_Gamma()
{
    ylm1 = new double[100];
    ylm2 = new double[100]; // can used for L=9
    iq = new int[1];
    x0 = new double[1];
    x1 = new double[1];
    x2 = new double[1];
    x3 = new double[1];
    x12 = new double[1];
    x03 = new double[1];
}

Gint_Gamma::~Gint_Gamma()
{
    delete[] ylm1;
    delete[] ylm2;
    delete[] iq;
    delete[] x0;
    delete[] x1;
    delete[] x2;
    delete[] x3;
    delete[] x12;
    delete[] x03;
}


void Gint_Gamma::save_atoms_on_grid(const Grid_Technique &gt)
{
    TITLE("Grid_Integral","save_atoms_on_grid");

    // mohan change.
    max_size = gt.max_atom;

	if(max_size == 0)
	{
		// mohan add return 2011-03-15
		ofs_warning << " processor " << MY_RANK << ": no atom on sub-fft-grid." << endl;
		return;
	}

    delete[] iq;
    delete[] x0;
    delete[] x1;
    delete[] x2;
    delete[] x3;
    delete[] x12;
    delete[] x03;
    this->iq = new int[max_size];
    this->x0 = new double[max_size];
    this->x1 = new double[max_size];
    this->x2 = new double[max_size];
    this->x3 = new double[max_size];
    this->x12 = new double[max_size];
    this->x03 = new double[max_size];

	ZEROS(iq, max_size);
	ZEROS(x0, max_size);
	ZEROS(x1, max_size);
	ZEROS(x2, max_size);
	ZEROS(x3, max_size);
	ZEROS(x12, max_size);
	ZEROS(x03, max_size);

	this->vfactor = std::abs(this->latvec0.Det())/gt.ncxyz;

    //OUT(ofs_running,"Max atom number on sub-FFT-grid",max_size);
    //ofs_running << "\n dense(DIY) = " << dense;
    //ofs_running << "\n count_dense = " << (double)count_dense/nxyz*100 << "%";
    //ofs_running << "\n count_sparse = " << (double)count_sparse/nxyz*100 << "%" << endl;

    return;
}

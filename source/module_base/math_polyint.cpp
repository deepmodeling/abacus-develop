#include "math_polyint.h"
#include "timer.h"
namespace ModuleBase
{

void PolyInt::Polynomial_Interpolation
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
    ModuleBase::realArray &y,
    const int &dim_y,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
    ModuleBase::timer::tick("PolyInt","Poly_Interpo_1");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
	if(iq>=table_length-4)
	{
		std::cout << "\n x = " << x;
		std::cout << "\n iq = " << iq << " table_length = " << table_length << std::endl;
	}	
    assert(iq < table_length-4);

    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    y(dim1, dim2, dim_y)=
        table(dim1, dim2, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, iq+3) * x1 * x2 * x0 / 6.0 ;

    ModuleBase::timer::tick("PolyInt","Poly_Interpo_1");
    return;
}

double PolyInt::Polynomial_Interpolation
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	ModuleBase::timer::tick("PolyInt","Poly_Interpo_2");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);

    if (iq > table_length - 4)
    {
        std::cout << "\n x = " << x;
        std::cout << "\n table_interval = " << table_interval;
        std::cout << "\n iq=" << iq << " table_length = " << table_length << std::endl;
        std::cout << "\n Not enough space allocated for radial FFT." << std::endl;
        std::cout << " It is due to the rapid change of the size of cell:" << std::endl;
        std::cout << " Try reseting a larger cell_factor parameter in INPUT" << std::endl; // LiuXh add 20180619
        std::cout << " Or try reseting a smaller relax_scale_force parameter in INPUT\n" << std::endl;
        exit(0);
    }

    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
        table(dim1, dim2, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, iq+3) * x1 * x2 * x0 / 6.0 ;

//	ModuleBase::timer::tick("PolyInt","Poly_Interpo_2");
    return y;
}

double PolyInt::Polynomial_Interpolation            // pengfei Li 2018-3-23
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
	const int &dim3,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	ModuleBase::timer::tick("PolyInt","Poly_Interpo_3");
    assert(table_interval>0.0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
    
	if(iq>table_length-4)
	{
		std::cout << "\n x = " << x;
		std::cout << "\n table_interval = " << table_interval;
		std::cout << "\n iq=" << iq << " table_length = " << table_length << std::endl;
	}
	assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;
    const double y=
        table(dim1, dim2, dim3, iq)   * x1 * x2 * x3 / 6.0 +
        table(dim1, dim2, dim3, iq+1) * x0 * x2 * x3 / 2.0 -
        table(dim1, dim2, dim3, iq+2) * x1 * x0 * x3 / 2.0 +
        table(dim1, dim2, dim3, iq+3) * x1 * x2 * x0 / 6.0 ;

//	ModuleBase::timer::tick("PolyInt","Poly_Interpo_3");
    return y;
}

double PolyInt::Polynomial_Interpolation
(
    const double *table,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	assert(table_interval>0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
//	if(iq >= table_length-4)
//		std::cout << "\n iq = " << iq << " table_length = " << table_length;
  
   assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;

    /*
    const double y=
    	table[iq]   * x1 * x2 * x3 / 6.0 +
    	table[iq+1] * x0 * x2 * x3 / 2.0 -
    	table[iq+2] * x1 * x0 * x3 / 2.0 +
    	table[iq+3] * x1 * x2 * x0 / 6.0 ;
    	*/

    return x1*x2*(table[iq]*x3+table[iq+3]*x0)/6.0
         + x0*x3*(table[iq+1]*x2-table[iq+2]*x1)/2.0;
}

double PolyInt::Polynomial_Interpolation_xy
(
    const double *xpoint,
    const double *ypoint,
    const int table_length,
    const double &x             // input value
)
{
    int position = -1;

    if (x < xpoint[0])
    {
        return ypoint[0];
    }
    // ModuleBase::timer::tick("PolyInt","Poly_Inter_xy");

    for (int ik = 0; ik < table_length; ik++)
    {
        if (x < xpoint[ik])
        {
            break;
        }
        else
        {
            position ++;
        }
    }

    assert(position >= 0);
    assert(position <= table_length-1);

    if (position + 6 < table_length)
    {
        double dx1, dx2, dx3, dx4, dx5, dx6;
        dx1 = x - xpoint[position];
        dx2 = x - xpoint[position+1];
        dx3 = x - xpoint[position+2];
        dx4 = x - xpoint[position+3];
        dx5 = x - xpoint[position+4];
        dx6 = x - xpoint[position+5];


        double x12, x13, x14, x15, x16, x23, x24, x25, x26, x34, x35, x36, x45, x46, x56;
        x12 = xpoint[position] - xpoint[position+1];
        x13 = xpoint[position] - xpoint[position+2];
        x14 = xpoint[position] - xpoint[position+3];
        x15 = xpoint[position] - xpoint[position+4];
        x16 = xpoint[position] - xpoint[position+5];


        x23 = xpoint[position+1] - xpoint[position+2];
        x24 = xpoint[position+1] - xpoint[position+3];
        x25 = xpoint[position+1] - xpoint[position+4];
        x26 = xpoint[position+1] - xpoint[position+5];

        x34 = xpoint[position+2] - xpoint[position+3];
        x35 = xpoint[position+2] - xpoint[position+4];
        x36 = xpoint[position+2] - xpoint[position+5];

        x45 = xpoint[position+3] - xpoint[position+4];
        x46 = xpoint[position+3] - xpoint[position+5];

        x56 = xpoint[position+4] - xpoint[position+5];

        double part1, part2, part3, part4, part5, part6;
        part1 = dx2 * dx3 * dx4 * dx5 * dx6 / x12 / x13 / x14 / x15 / x16 * ypoint[position];
        part2 = dx1 * dx3 * dx4 * dx5 * dx6 / (-x12) / x23 / x24 / x25 / x26 * ypoint[position+1];
        part3 = dx1 * dx2 * dx4 * dx5 * dx6 / (-x13) / (-x23) / x34 / x35 / x36 * ypoint[position+2];
        part4 = dx1 * dx2 * dx3 * dx5 * dx6 / (-x14) / (-x24) / (-x34) / x45 / x46 * ypoint[position+3];
        part5 = dx1 * dx2 * dx3 * dx4 * dx6 / (-x15) / (-x25) / (-x35) / (-x45) / x56 * ypoint[position+4];
        part6 = dx1 * dx2 * dx3 * dx4 * dx5 / (-x16) / (-x26) / (-x36) / (-x46) / (-x56) * ypoint[position+5];

        // 	ModuleBase::timer::tick("PolyInt","Poly_Inter_xy");
        return part1 + part2 + part3 + part4 + part5 + part6;
    }
    else
    {
        // 	ModuleBase::timer::tick("PolyInt","Poly_Inter_xy");
        return ypoint[position];
    }
}

}

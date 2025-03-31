#include "symmetry.h"
using namespace ModuleSymmetry;

void Symmetry::hermite_normal_form(const ModuleBase::Matrix3 &s3, 
		ModuleBase::Matrix3 &h3, 
		ModuleBase::Matrix3 &b3) const
{
    ModuleBase::TITLE("Symmetry","hermite_normal_form");
    // check the non-singularity and integer elements of s
#ifdef __DEBUG
    assert(!equal(s3.Det(), 0.0));
#endif

	auto near_equal = [this](double x, double y) 
	{
		return fabs(x - y) < 10 * epsilon;
	};

    ModuleBase::matrix s = s3.to_matrix();

	for (int i = 0; i < 3; ++i) 
	{
		for (int j = 0;j < 3;++j)
		{
            double sij_round = std::round(s(i, j));
#ifdef __DEBUG
            assert(near_equal(s(i, j), sij_round));
#endif
            s(i, j) = sij_round;
        }
    }

    // convert Matrix3 to matrix
    ModuleBase::matrix h=s, b(3, 3, true);
    b(0, 0)=1; 
    b(1, 1)=1; 
    b(2, 2)=1;

    // transform H into lower triangular form
    auto max_min_index = [&h, this](int row, int &i1_to_max, int &i2_to_min)
    {
        if(fabs(h(row, i1_to_max)) < fabs(h(row, i2_to_min)) - epsilon)
        {
            int tmp = i2_to_min;
            i2_to_min = i1_to_max;
            i1_to_max = tmp;
        }
        return;
    };

    auto max_min_index_row1 = [&max_min_index, &h, this](int &imax, int &imin)
    {
        int imid=1;
        imax=0; imin=2;
        max_min_index(0, imid, imin);
        max_min_index(0, imax, imid);
		max_min_index(0, imid, imin);
		if (equal(h(0, imin), 0)) 
		{
			imin = imid;
		} 
		else if (equal(h(0, imax), 0)) 
		{
			imax = imid;
		}
        return;
    };

    auto swap_col = [&h, &b](int c1, int c2)
    {
        double tmp=0.0;
        for(int r=0;r<3;++r)
        {
            tmp = h(r, c2);
            h(r, c2)=h(r, c1);
            h(r, c1)=tmp;
            tmp = b(r, c2);
            b(r, c2)=b(r, c1);
            b(r, c1)=tmp;
        }
        return;
    };

    // row 1 
    int imax=0;
    int imin=0;

    while(int(equal(h(0, 0), 0)) + int(equal(h(0, 1), 0)) + int(equal(h(0, 2), 0)) < 2)
    {
        max_min_index_row1(imax, imin);
        double f = floor((fabs(h(0, imax) )+ epsilon)/fabs(h(0, imin)));
		if (h(0, imax) * h(0, imin) < -epsilon) 
		{
			f *= -1;
		}
		for(int r=0;r<3;++r) 
		{
			h(r, imax) -= f*h(r, imin); 
			b(r, imax) -= f*b(r, imin); 
		}
    }

	if (equal(h(0, 0), 0)) 
	{
		equal(h(0, 1), 0) ? swap_col(0, 2) : swap_col(0, 1);
	}

	if (h(0, 0) < -epsilon) 
	{
		for (int r = 0; r < 3; ++r) 
		{
			h(r, 0) *= -1;
            b(r, 0) *= -1;
        }
    }

    //row 2
	if (equal(h(1, 1), 0)) 
	{
		swap_col(1, 2);
	}

    while(!equal(h(1, 2), 0))
    {
        imax=1, imin=2;
        max_min_index(1, imax, imin);
        double f = floor((fabs(h(1, imax) )+ epsilon)/fabs(h(1, imin)));

		if (h(1, imax) * h(1, imin) < -epsilon) 
		{
			f *= -1;
		}

		for(int r=0;r<3;++r) 
		{
			h(r, imax) -= f*h(r, imin); 
			b(r, imax) -= f*b(r, imin); 
		}

		if (equal(h(1, 1), 0)) 
		{
			swap_col(1, 2);
		}
    }
	if (h(1, 1) < -epsilon) 
	{
		for (int r = 0; r < 3; ++r) 
		{
            h(r, 1) *= -1;
            b(r, 1) *= -1;
        }
    }

    //row3
    if (h(2, 2) < -epsilon) 
	{
		for (int r = 0; r < 3; ++r) 
		{
			h(r, 2) *= -1;
            b(r, 2) *= -1;
        }
    }

    // deal with off-diagonal elements
	while (h(1, 0) > h(1, 1) - epsilon) 
	{
		for(int r=0;r<3;++r) 
		{
			h(r, 0) -= h(r, 1); 
			b(r, 0) -= b(r, 1);
		}
	}
	while (h(1, 0) < -epsilon) 
	{
		for(int r=0;r<3;++r) 
		{
			h(r, 0) += h(r, 1); 
			b(r, 0) += b(r, 1);
		}
	}
    for(int j=0;j<2;++j)
    {
		while (h(2, j) > h(2, 2) - epsilon) 
		{
			for(int r=0;r<3;++r) 
			{
				h(r, j) -= h(r, 2); 
				b(r, j) -= b(r, 2);
			}
        }
		while (h(2, j) < -epsilon) 
		{
			for(int r=0;r<3;++r) 
			{
				h(r, j) += h(r, 2); 
				b(r, j) += b(r, 2);
			}
        }
    }

    //convert matrix to Matrix3
    h3.e11=h(0, 0); h3.e12=h(0, 1); h3.e13=h(0, 2);
    h3.e21=h(1, 0); h3.e22=h(1, 1); h3.e23=h(1, 2);
    h3.e31=h(2, 0); h3.e32=h(2, 1); h3.e33=h(2, 2);
    b3.e11=b(0, 0); b3.e12=b(0, 1); b3.e13=b(0, 2);
    b3.e21=b(1, 0); b3.e22=b(1, 1); b3.e23=b(1, 2);
    b3.e31=b(2, 0); b3.e32=b(2, 1); b3.e33=b(2, 2);

    //check s*b=h
    ModuleBase::matrix check_zeros = s3.to_matrix() * b - h;
#ifdef __DEBUG
	for (int i = 0;i < 3;++i)
	{
		for(int j=0;j<3;++j)
		{
			assert(near_equal(check_zeros(i, j), 0));
		}
	}
#endif
    return;
}

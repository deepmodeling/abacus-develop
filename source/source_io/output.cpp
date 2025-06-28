/* output.cpp */

#include "output.h"

void output::printrm(std::ofstream &ofs,const std::string &s, const ModuleBase::matrix &m, const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n " << s << " nr=" << b1 << " nc=" << b2 ;
    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
//        ofs<<"\n row = "<<i;
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) ofs << "\n ";

			if (std::abs(m(i,j)) > limit)
			{
				ofs << std::setprecision(6) << std::setw(12) << m(i,j);
            }
			else
			{
				ofs << std::setw(12) << "0";
			}
        }
    }
	ofs << std::endl;
    return;
}

void output::printrm(const std::string &s, const ModuleBase::matrix &m, const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n " << s << " nr=" << b1 << " nc=" << b2 ;
    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
		//std::cout << "\n row=" << i;
        for (int j = 0;j < b2;j++)
        {
//            if (j % 4 == 0) std::cout << "\n ";
			if (j % 8 == 0) std::cout << "\n ";
//            if (std::abs(m(i,j)) > limit) std::cout << std::setprecision(15) << std::setw(20) << m(i,j);
            if (std::abs(m(i,j)) > limit) std::cout << std::setprecision(6) << std::setw(12) << m(i,j);
//            else std::cout<<std::setw(20)<<"0";
            else std::cout<<std::setw(12)<<"0";
        }
//		std::cout << "\n";
    }

	std::cout << std::endl;
    return;
}

/*
void output::printcm_norm(std::ofstream &ofs, const std::string &s, const ModuleBase::ComplexMatrix &m, const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n" << s << " nr=" <<  b1 << " nc=" << b2 ;

    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) ofs << "\n ";
            double norm = ( conj( m(i, j) ) * m(i, j) ).real();
            if (std::abs(norm) > limit) ofs<< std::setw(12) << sqrt(norm);
            else ofs<<std::setw(12)<<"0";
        }
    }
    return;
}
*/

/*
void output::printcm_norm(const std::string &s, const ModuleBase::ComplexMatrix &m, const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n" << s << " nr=" <<  b1 << " nc=" << b2 ;

    if (b1*b2 == 0) return;
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) std::cout << "\n ";
            double norm = ( conj( m(i, j) ) * m(i, j) ).real();
            if (std::abs(norm) > limit) std::cout<< std::setw(12) << sqrt(norm);
            else std::cout<<std::setw(12)<<"0";
        }
    }
    return;
}
*/


/*
void output::printcm(std::ofstream &ofs, const std::string &s, const ModuleBase::ComplexMatrix &m)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    ofs << "\n" << s << " " <<  b1 << " " << b2 ;

    if (b1*b2 == 0) return;

    ofs.setf(std::ios::scientific, std::ios::floatfield);
    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j%2 == 0) ofs << "\n ";
            ofs << std::setw(20) << m(i, j).real() << std::setw(20) << m(i, j).imag();
        }
    }
    ofs.unsetf(std::ios::scientific);
    return;
}//end print cm
*/

/*
void output::printcm(const std::string &s, const ModuleBase::ComplexMatrix &m)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n" << s << " nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 4 == 0) std::cout << "\n ";
            std::cout << std::setw(18) << m(i, j);
        }
    }
    return;
}
*/

/*
void output::printcm_real_limit_hermit(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n" << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) std::cout << "\n ";
            const double x = m(i,j).real();
            if (std::abs(x) < limit || i < j) std::cout<<std::setw(12)<<"0";
            else std::cout << std::setw(12) << x;
        }
    }
    return;
}
*/

/*
void output::printcm_real(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n " << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) std::cout << "\n ";
            const double x = m(i,j).real();
            if (std::abs(x) < limit ) std::cout<<std::setw(12)<<"0";
            else std::cout << std::setw(12) << x;
        }
    }
    return;
}
*/


/*
void output::printcm_imag(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit)
{
    const int b1 = m.nr;
    const int b2 = m.nc;
    std::cout << "\n " << s << "  nr = " <<  b1 << " nc = " << b2 ;

    if (b1*b2 == 0) return;

    for (int i = 0;i < b1;i++)
    {
        for (int j = 0;j < b2;j++)
        {
            if (j % 8 == 0) std::cout << "\n ";
            const double x = m(i,j).imag();
            if (std::abs(x) < limit ) std::cout<<std::setw(12)<<"0";
            else std::cout << std::setw(12) << x;
        }
    }
    return;
}
*/





void output::printr3_d(std::ofstream &ofs, const std::string &s,const ModuleBase::realArray &u)
{
    const int b1 = u.getBound1();
    const int b2 = u.getBound2();
    const int b3 = u.getBound3();
    ofs << "\n\n " << s << "  b1 = " << b1
    << " b2 = " << b2 << " b3 = " << b3;

    if (u.getSize() == 0) return;

    for (int i = 0;i < b1;i++)
        for (int j = 0;j < b2;j++)
            for (int k = 0;k < b3;k++)
            {
                if (k % 4 == 0)
                {
                    ofs << "\n";
                }
                ofs << std::setw(18) << u(i, j, k);
            }
    return;
}//end printr3_d

/*
void output::printr4_d(std::ofstream &ofs, const std::string &s,const ModuleBase::realArray &u)
{
    const int b1 = u.getBound1();
    const int b2 = u.getBound2();
    const int b3 = u.getBound3();
    const int b4 = u.getBound4();
    ofs << "\n\n " << s << "  b1 = " << b1 << " b2 = " << b2
    << " b3 = " << b3 << " b4 = " << b4;

    if (u.getSize() == 0) return;

    for (int i = 0;i < b1;i++)
        for (int j = 0;j < b2;j++)
            for (int k = 0;k < b3;k++)
                for (int m = 0;m < b4;m++)
                {
                    if (m % 4 == 0)
                    {
                        ofs << "\n";
                    }
                    ofs << std::setw(15) << u(i, j, k, m);
                }
}//end print4_d
*/

void output::printM3(std::ofstream &ofs,const std::string &description, const ModuleBase::Matrix3 &m)
{
    ofs << " " << description << std::endl;
	ofs << std::setiosflags(std::ios::showpos);
    ofs << " " << std::setw(20) << m.e11 << std::setw(20) << m.e12 << std::setw(20) << m.e13
	<< "\n " << std::setw(20) << m.e21 << std::setw(20) << m.e22 << std::setw(20) << m.e23
	<< "\n " << std::setw(20) << m.e31 << std::setw(20) << m.e32 << std::setw(20) << m.e33 << std::endl;
	ofs << std::resetiosflags(std::ios::showpos);
    return;
}

void output::printM3(const std::string &description, const ModuleBase::Matrix3 &m)
{
    std::cout << "\n " << description << std::endl;
    std::cout << std::setw(20) << m.e11 << std::setw(20) << m.e12 << std::setw(20) << m.e13 << "\n"
         << std::setw(20) << m.e21 << std::setw(20) << m.e22 << std::setw(20) << m.e23 << "\n"
         << std::setw(20) << m.e31 << std::setw(20) << m.e32 << std::setw(20) << m.e33 << std::endl;
    return;
}

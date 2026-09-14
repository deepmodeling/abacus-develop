#include "symmetry.h"
#include "source_base/mymath.h"

bool ModuleSymmetry::test_brav = 0;

namespace ModuleSymmetry
{
// Find the type of bravais lattice.
std::string Symmetry_Basic::get_brav_name(const int ibrav) const
{
    switch(ibrav)
    {
        case  1: return "01. Cubic P (simple)";
        case  2: return "02. Cubic I (body-centered)";
        case  3: return "03. Cubic F (face-centered)";
        case  4: return "04. Hexagonal cell";
        case  5: return "05. Tetrogonal P (simple)";
        case  6: return "06. Tetrogonal I (body-centered)";
        case  7: return "07. Rhombohedral (Trigonal) cell";
        case  8: return "08. Orthorhombic P(simple)";
        case  9: return "09. Orthorhombic I (body-centered)";
        case 10: return "10. Orthorhombic F (face-centered)";
        case 11: return "11. Orthorhombic C (base-centered)";
        case 12: return "12. Monoclinic P (simple)";
        case 13: return "13. Monoclinic A (base-center)";
        case 14: return "14. Triclinic cell";
        case 15: return "wrong !! ";
    }
    // return "Congratulations! You have found a bravais lattice that never existed!";
    return "Unknown Bravais lattice";
}

// Control the accuracy
bool Symmetry_Basic::equal(const double &m, const double &n) const
{
    //if( fabs(m-n) < 1.0e-5 )
    if (fabs(m-n) < epsilon) //LiuXh add 2021-08-12, use accuracy for symmetry
    {
        return true;
    }
    return false;
}

// check the boundary condition of atom positions.
void Symmetry_Basic::check_boundary(double &x)const
{
    if(equal(x,-0.5) || equal(x,0.5)) x=-0.5;
}

double Symmetry_Basic::get_translation_vector(const double& x1, const double& x2) const
{
    double t=0.0; // "t"ranslation
    t = x2 - x1;
    t = fmod(t+100.0, 1.0);
    if( fabs(t-1) < epsilon * 0.5) { t = 0.0; }
    return t;
}

void Symmetry_Basic::check_translation(double &x, const double &t) const
{
    x += t;
    //impose the periodic boundary condition
    x = fmod(x + 100.5,1) - 0.5;
    return;
}

double Symmetry_Basic::check_diff(const double& x1, const double& x2)const
{
    double diff = x1 - x2;
    diff = fmod(diff + 100,1);
    //for reasons of safety
    if(fabs(diff - 1.0) < epsilon)
    {
        diff = 0;
    }
    return diff;
}


void Symmetry_Basic::order_atoms(double* pos, const int& nat, const int* index) const
{
    double** tmp = new double*[nat];
    for(int ia=0; ia<nat; ia++)
    {
        tmp[ia] = new double[3];
    }

    for(int ia=0; ia<nat; ia++)
    {
        tmp[ia][0] = pos[index[ia]*3+0];
        tmp[ia][1] = pos[index[ia]*3+1];
        tmp[ia][2] = pos[index[ia]*3+2];
    }

    for(int ia=0; ia<nat; ia++)
    {
        pos[ia*3+0] = tmp[ia][0];
        pos[ia*3+1] = tmp[ia][1];
        pos[ia*3+2] = tmp[ia][2];
    }

    for(int ia=0; ia<nat; ia++)
    {
        delete[] tmp[ia];
    }
    delete[] tmp;

    return;
}

// convert a set of vectors (va) represented in the basis vectors old1, old2, old3
// to a set of vectors (vb) represented in the basis vectors new1, new2, new3
void Symmetry_Basic::veccon(
        double *carpos,
        double *rotpos,
        const int num,
        const ModuleBase::Vector3<double> &old1,
        const ModuleBase::Vector3<double> &old2,
        const ModuleBase::Vector3<double> &old3,
        const ModuleBase::Vector3<double> &new1,
        const ModuleBase::Vector3<double> &new2,
        const ModuleBase::Vector3<double> &new3
        )
{

    GlobalV::ofs_running << "\n old1:" << old1.x << " " << old1.y << " " << old1.z;
    GlobalV::ofs_running << "\n old2:" << old2.x << " " << old2.y << " " << old2.z;
    GlobalV::ofs_running << "\n old3:" << old3.x << " " << old3.y << " " << old3.z;

    GlobalV::ofs_running << "\n new1:" << new1.x << " " << new1.y << " " << new1.z;
    GlobalV::ofs_running << "\n new2:" << new2.x << " " << new2.y << " " << new2.z;
    GlobalV::ofs_running << "\n new3:" << new3.x << " " << new3.y << " " << new3.z;

    ModuleBase::Matrix3 oldlat;
    oldlat.e11 = old1.x;
    oldlat.e12 = old1.y;
    oldlat.e13 = old1.z;
    oldlat.e21 = old2.x;
    oldlat.e22 = old2.y;
    oldlat.e23 = old2.z;
    oldlat.e31 = old3.x;
    oldlat.e32 = old3.y;
    oldlat.e33 = old3.z;

    ModuleBase::Matrix3 newlat;
    newlat.e11 = new1.x;
    newlat.e12 = new1.y;
    newlat.e13 = new1.z;
    newlat.e21 = new2.x;
    newlat.e22 = new2.y;
    newlat.e23 = new2.z;
    newlat.e31 = new3.x;
    newlat.e32 = new3.y;
    newlat.e33 = new3.z;

    ModuleBase::Matrix3 GT = newlat.Inverse();

    ModuleBase::Vector3<double> car;
    ModuleBase::Vector3<double> direct_old;
    ModuleBase::Vector3<double> direct_new;

    //calculate the reciprocal vectors rb1, rb2, rb3 for the vectors new1, new2, new3
    //this->recip(1.0, new1, new2, new3, rb1, rb2, rb3);

    for(int i = 0; i < num; ++i)
    {
        direct_old.x = carpos[i * 3 + 0];
        direct_old.y = carpos[i * 3 + 1];
        direct_old.z = carpos[i * 3 + 2];

        car = direct_old * oldlat;
        direct_new = car * GT;

        rotpos[i * 3 + 0] = direct_new.x;
        rotpos[i * 3 + 1] = direct_new.y;
        rotpos[i * 3 + 2] = direct_new.z;
    }
    return;
}

void Symmetry_Basic::rotate( ModuleBase::Matrix3 &gmatrix, ModuleBase::Vector3<double> &gtrans,
        int i, int j, int k, // FFT grid index.
        const int nr1, const int nr2, const int nr3, // dimension of FFT grid.
        int &ri, int &rj, int &rk)
{
    static ModuleBase::Matrix3 g;
    g.e11 = gmatrix.e11;
    g.e21 = gmatrix.e21 * (double)nr1 / (double)nr2;
    g.e31 = gmatrix.e31 * (double)nr1 / (double)nr3;
    g.e12 = gmatrix.e12 * (double)nr2 / (double)nr1;
    g.e22 = gmatrix.e22;
    g.e32 = gmatrix.e32 * (double)nr2 / (double)nr3;
    g.e13 = gmatrix.e13 * (double)nr3 / (double)nr1;
    g.e23 = gmatrix.e23 * (double)nr3 / (double)nr2;
    g.e33 = gmatrix.e33;

    ri = int(g.e11 * i + g.e21 * j + g.e31 * k) + (int)(gtrans.x *  nr1);
    if (ri < 0)
    {
        ri += 10 * nr1;
    }
    ri = ri%nr1;
    rj = static_cast<int>(g.e12 * i + g.e22 * j + g.e32 * k) + static_cast<int>(gtrans.y  * nr2);
    if (rj < 0)
    {
        rj += 10 * nr2;
    }
    rj = rj%nr2;
    rk = static_cast<int>(g.e13 * i + g.e23 * j + g.e33 * k) + static_cast<int>(gtrans.z  * nr3);
    if (rk < 0)
    {
        rk += 10 * nr3;
    }
    rk = rk%nr3;
    return;
}
}

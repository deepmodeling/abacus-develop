#include <memory>
#include <array>
#include "symmetry.h"

namespace ModuleSymmetry
{
int Symmetry::symm_flag = 0;
bool Symmetry::symm_autoclose = false;
bool Symmetry::pricell_loop = true;

void Symmetry::set_atom_map(const Atom* atoms)
{
    ModuleBase::TITLE("Symmetry", "set_atom_map");
    if (this->isym_rotiat_.size() == this->nrotk) {
        return;
    }
    this->isym_rotiat_.resize(this->nrotk);
    for (int i = 0; i < this->nrotk; ++i) {
        this->isym_rotiat_[i].resize(this->nat, -1);
    }

    double* pos = this->newpos;
    double* rotpos = this->rotpos;
    ModuleBase::GlobalFunc::ZEROS(pos, this->nat * 3);
    int iat = 0;
    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = 0; ia < this->na[it]; ia++)
        {
            pos[3 * iat] = atoms[it].taud[ia].x;
            pos[3 * iat + 1] = atoms[it].taud[ia].y;
            pos[3 * iat + 2] = atoms[it].taud[ia].z;
            for (int k = 0; k < 3; ++k)
            {
                this->check_translation(pos[iat * 3 + k], -floor(pos[iat * 3 + k]));
                this->check_boundary(pos[iat * 3 + k]);
            }
            iat++;
        }
    }
    for (int it = 0; it < this->ntype; it++)
    {
        for (int ia = istart[it]; ia < istart[it] + na[it]; ++ia)
        {
			const int xx = ia * 3; 
			const int yy = ia * 3 + 1; 
			const int zz = ia * 3 + 2;

			for (int k = 0;k < this->nrotk;++k)
            {
				rotpos[xx] = pos[xx] * gmatrix[k].e11 
					+ pos[yy] * gmatrix[k].e21 
					+ pos[zz] * gmatrix[k].e31 + gtrans[k].x;
				rotpos[yy] = pos[xx] * gmatrix[k].e12 
					+ pos[yy] * gmatrix[k].e22 
					+ pos[zz] * gmatrix[k].e32 + gtrans[k].y;
				rotpos[zz] = pos[xx] * gmatrix[k].e13 
					+ pos[yy] * gmatrix[k].e23 
					+ pos[zz] * gmatrix[k].e33 + gtrans[k].z;

                check_translation(rotpos[xx], -floor(rotpos[xx]));
                check_boundary(rotpos[xx]);
                check_translation(rotpos[yy], -floor(rotpos[yy]));
                check_boundary(rotpos[yy]);
                check_translation(rotpos[zz], -floor(rotpos[zz]));
                check_boundary(rotpos[zz]);

                for (int ja = istart[it]; ja < istart[it] + na[it]; ++ja)
                {
                    double diff1 = check_diff(pos[ja * 3], rotpos[xx]);
                    double diff2 = check_diff(pos[ja * 3 + 1], rotpos[yy]);
                    double diff3 = check_diff(pos[ja * 3 + 2], rotpos[zz]);
                    if (equal(diff1, 0.0) && equal(diff2, 0.0) && equal(diff3, 0.0))
                    {
                        this->isym_rotiat_[k][ia] = ja;

                        break;
                    }
                }
            }
        }
    }
}

void Symmetry::symmetrize_vec3_nat(double* v)const   // pengfei 2016-12-20
{
    ModuleBase::TITLE("Symmetry", "symmetrize_vec3_nat");
    double* vtot;
    int* n;
    vtot = new double[nat * 3]; ModuleBase::GlobalFunc::ZEROS(vtot, nat * 3);
    n = new int[nat]; ModuleBase::GlobalFunc::ZEROS(n, nat);

    for (int j = 0;j < nat; ++j)
    {
        const int jx = j * 3; const int jy = j * 3 + 1; const int jz = j * 3 + 2;
        for (int k = 0; k < nrotk; ++k)
        {
            int l = this->isym_rotiat_[k][j];
            if (l < 0) {
                continue;
            }
            vtot[l*3] = vtot[l*3] + v[jx] * gmatrix[k].e11 + v[jy] * gmatrix[k].e21 + v[jz] * gmatrix[k].e31;
            vtot[l*3+1] = vtot[l*3+1] + v[jx] * gmatrix[k].e12 + v[jy] * gmatrix[k].e22 + v[jz] * gmatrix[k].e32;
            vtot[l*3+2] = vtot[l*3+2] + v[jx] * gmatrix[k].e13 + v[jy] * gmatrix[k].e23 + v[jz] * gmatrix[k].e33;
            n[l]++;
        }
	}
    for (int j = 0;j < nat; ++j)
    {
        v[j * 3] = vtot[j * 3] / n[j];
        v[j * 3 + 1] = vtot[j * 3 + 1] / n[j];
        v[j * 3 + 2] = vtot[j * 3 + 2] / n[j];
    }
    delete[] vtot;
    delete[] n;
	return;
}

void Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat)const   //zhengdy added 2017
{
    ModuleBase::matrix A = lat.latvec.to_matrix();
    ModuleBase::matrix AT = lat.latvec.Transpose().to_matrix();
    ModuleBase::matrix invA = lat.GT.to_matrix();
    ModuleBase::matrix invAT = lat.G.to_matrix();
    ModuleBase::matrix tot_sigma(3, 3, true);
    sigma = A * sigma * AT;
    for (int k = 0; k < nrotk; ++k) {
        tot_sigma += invA * gmatrix[k].to_matrix() * sigma
                     * gmatrix[k].Transpose().to_matrix() * invAT;
    }
    sigma = tot_sigma * static_cast<double>(1.0 / nrotk);
	return;
}

void Symmetry::gmatrix_convert_int(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b) const
{
    auto round = [](double x){return (x>0.0)?floor(x+0.5):ceil(x-0.5);};
    ModuleBase::Matrix3 ai = a.Inverse();
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          sb[i]=b*ai*sa[i]*a*bi;
          //to int 
          sb[i].e11=round(sb[i].e11);
          sb[i].e12=round(sb[i].e12);
          sb[i].e13=round(sb[i].e13);
          sb[i].e21=round(sb[i].e21);
          sb[i].e22=round(sb[i].e22);
          sb[i].e23=round(sb[i].e23);
          sb[i].e31=round(sb[i].e31);
          sb[i].e32=round(sb[i].e32);
          sb[i].e33=round(sb[i].e33);
    }
}

void Symmetry::gmatrix_convert(const ModuleBase::Matrix3* sa, ModuleBase::Matrix3* sb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const
{
    ModuleBase::Matrix3 ai = a.Inverse();
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          sb[i]=b*ai*sa[i]*a*bi;
    }
}

void Symmetry::gtrans_convert(const ModuleBase::Vector3<double>* va, ModuleBase::Vector3<double>* vb, 
        const int n, const ModuleBase::Matrix3 &a, const ModuleBase::Matrix3 &b)const
{
    ModuleBase::Matrix3 bi = b.Inverse();
    for (int i=0;i<n;++i)
    {
          vb[i]=va[i]*a*bi;
    }
}

void Symmetry::gmatrix_invmap(const ModuleBase::Matrix3* s, const int n, int* invmap) const
{
    ModuleBase::Matrix3 eig(1, 0, 0, 0, 1, 0, 0, 0, 1);
    ModuleBase::Matrix3 tmp;
    for (int i=0;i<n;++i)
    {
        for (int j=i;j<n;++j)
        {
            tmp=s[i]*s[j];
            if(equal(tmp.e11, 1) && equal(tmp.e22, 1) && equal(tmp.e33, 1) &&
                equal(tmp.e12, 0) && equal(tmp.e21, 0) && equal(tmp.e13, 0) &&
                equal(tmp.e31, 0) && equal(tmp.e23, 0) && equal(tmp.e32, 0))
            {
                invmap[i]=j;
                invmap[j]=i;
                break;
            }
        }
    }
}

void Symmetry::get_shortest_latvec(ModuleBase::Vector3<double> &a1, 
        ModuleBase::Vector3<double> &a2, ModuleBase::Vector3<double> &a3) const
{
    double len1=a1.norm();
    double len2=a2.norm();
    double len3=a3.norm();
    bool flag=true; //at least one iter
    auto loop = [this, &flag](ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double>&v2, double &len)
    {
        bool fa=false, fb=false;
        // loop a
        double tmp_len=(v1-v2).norm();
        while (tmp_len < len-epsilon)
        {
            v1=v1-v2;
            len=v1.norm();
            tmp_len=(v1-v2).norm();
            fa=true;
        }
        // loop b
        tmp_len=(v1+v2).norm();
        while(tmp_len < len-epsilon)
        {
            assert(!fa);
            v1=v1+v2;
            len=v1.norm();
            tmp_len=(v1+v2).norm();
            fb=true;
        }
        if (fa || fb) {
            flag = true;
        }
        return;
    };
    while(flag) //iter
    {
        flag=false;
        // if any of a1, a2, a3 is updated, flag will become true.
        // which means a further search is needed.
        loop(a1, a2, len1);
        loop(a1, a3, len1);
        loop(a2, a1, len2);
        loop(a2, a3, len2);
        loop(a3, a1, len3);
        loop(a3, a2, len3);
    }
    return;
}

void Symmetry::get_optlat(ModuleBase::Vector3<double> &v1, ModuleBase::Vector3<double> &v2, 
        ModuleBase::Vector3<double> &v3, ModuleBase::Vector3<double> &w1, 
        ModuleBase::Vector3<double> &w2, ModuleBase::Vector3<double> &w3, 
        int& real_brav, double* cel_const, double* tmp_const) const
{
    ModuleBase::Vector3<double> r1, r2, r3;
    double cos1 = 1;
    double cos2 = 1;
    double cos3 = 1;
    int nif = 0;
    int ibrav = 0;
    for (int n33 = -2; n33 < 3; ++n33)
    {
        for (int n32 = -2; n32 < 3; ++n32)
        {
            for (int n31 = -2; n31 < 3; ++n31)
            {
                for (int n23 = -2; n23 < 3; ++n23)
                {
                    for (int n22 = -2; n22 < 3; ++n22)
                    {
                        for (int n21 = -2; n21 < 3; ++n21)
                        {
                            for (int n13 = -2; n13 < 3; ++n13)
                            {
                                for (int n12 = -2; n12 < 3; ++n12)
                                {
                                    for (int n11 = -2; n11 < 3; ++n11)
                                    {
                                        ModuleBase::Matrix3 mat(n11, n12, n13, n21, n22, n23, n31, n32, n33);

                                        if (equal(mat.Det(),1.0))
                                        {
                                            r1.x = n11 * v1.x + n12 * v2.x + n13 * v3.x;
                                            r1.y = n11 * v1.y + n12 * v2.y + n13 * v3.y;
                                            r1.z = n11 * v1.z + n12 * v2.z + n13 * v3.z;
                                     
									        r2.x = n21 * v1.x + n22 * v2.x + n23 * v3.x;
                                            r2.y = n21 * v1.y + n22 * v2.y + n23 * v3.y;
                                            r2.z = n21 * v1.z + n22 * v2.z + n23 * v3.z;
                                     
									        r3.x = n31 * v1.x + n32 * v2.x + n33 * v3.x;
                                            r3.y = n31 * v1.y + n32 * v2.y + n33 * v3.y;
                                            r3.z = n31 * v1.z + n32 * v2.z + n33 * v3.z;
											
                                            ibrav = standard_lat(r1, r2, r3, cel_const);

                                            if ( ibrav < real_brav || ( ibrav == real_brav
                                                    && ( fabs(cel_const[3]) < (cos1-1.0e-9) )
                                                    && ( fabs(cel_const[4]) < (cos2-1.0e-9) )
                                                    && ( fabs(cel_const[5]) < (cos3-1.0e-9) )) //mohan fix bug 2012-01-15, not <=
                                               )
                                            {
                                                real_brav = ibrav;
												
                                                cos1 = fabs(cel_const[3]);
                                                cos2 = fabs(cel_const[4]);
                                                cos3 = fabs(cel_const[5]);

                                                for (int i = 0; i < 6; ++i)
                                                {
                                                    tmp_const[i] = cel_const[i];
                                                }
                                                w1 = r1;
                                                w2 = r2;
                                                w3 = r3;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

bool Symmetry::is_all_movable(const Atom* atoms, const Statistics& st)const
{
    bool all_mbl = true;
    for (int iat = 0;iat < st.nat;++iat)
    {
        int it = st.iat2it[iat];
        int ia = st.iat2ia[iat];
        if (!atoms[it].mbl[ia].x || !atoms[it].mbl[ia].y || !atoms[it].mbl[ia].z)
        {
            all_mbl = false;
            break;
        }
    }
    return all_mbl;
}

}

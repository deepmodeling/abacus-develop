//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// DATA : 2006-11 ~ 2008-11
//==========================================================
#ifndef OUTPUT_H
#define OUTPUT_H

//#include "../src_pw/tools.h"
#include "../module_base/realarray.h"
#include "../module_base/matrix3.h"
#include "../module_base/complexmatrix.h"
class output
{
public:
    //============================
    // Print realArray (3D or 4D)
    //============================
    void printr3_d(ofstream &ofs,const string &s,const realArray &u) const;
    void printr4_d(ofstream &ofs,const string &s,const realArray &u) const;

    //===========================
    // print matrix3
    //===========================
    void printM3(ofstream &ofs,const string& description, const Matrix3 &m)const;
    void printM3(const string &description, const Matrix3 &m)const;

    //===============================
    // print matrix
    //===============================
    void printrm(ofstream &ofs,const string &s, const matrix &m, const double &limit = 1.0e-15) const;
    void printrm(const string &s, const matrix &m, const double &limit = 1.0e-15) const;


    //===============================
    // print ComplexMatrix
    //===============================
    void printcm(ofstream &ofs,const string &s, const ComplexMatrix &m) const;

    void printcm(const string &s, const ComplexMatrix &m) const;

    void printcm_real(const string &s, const ComplexMatrix &m,const double &limit = 1.0e-15) const;

    void printcm_real_limit_hermit(const string &s, const ComplexMatrix &m,const double &limit) const;

    void printcm_imag(const string &s, const ComplexMatrix &m,const double &limit = 1.0e-15) const;
    void printcm_norm(const string &s, const ComplexMatrix &m, const double &limit)const;
    void printcm_norm(ofstream &ofs, const string &s, const ComplexMatrix &m, const double &limit)const;

    //***************
    // Template
    //***************
public:

    template <class T>
    void printr1_d(ofstream &ofs, const string &s,T *u, int n1) const
    {
        ofs<<"\n\n "<<s<< "  n1 = "<< n1;
        if (n1<=0)return;
        if (u==0) return;
        //	ofs.setf(ios::scientific,ios::floatfield);
        for (int i=0;i<n1;i++)
        {
            if (i%8==0)ofs<<"\n";
            ofs<< setw(12)<<u[i];
        }
        //	ofs.unsetf(ios::scientific);
    }

    //===================================================
    // print one dimension array (double,int,string ...)
    //===================================================
    template <class T>
    void printr1_d(const string &s, T *u,const int n1) const
    {
        cout << "\n " << s << "  Dimension = " << n1;
        if (n1 <= 0)return;
        if (u == 0) return;
        //cout.setf(ios::scientific, ios::floatfield);
        for (int i = 0;i < n1;i++)
        {
            if (i % 8 == 0) cout << "\n";
            cout << setw(12) << u[i];
        }
        cout<<endl;
        //cout.unsetf(ios::scientific);
    }

    //================================
    // output Vector3 : printV3
    // output Vector3[] : printv31_d
    //================================
    template <class T>
    void printV3(ofstream &ofs, const Vector3 <T> v)const
    {
        ofs << " ";
        ofs << setw(18) << v.x << setw(18) << v.y << setw(18) << v.z << endl;
    }

    template <class T>
    void printV3(const Vector3 <T> v)const
    {
        cout << " ";
        cout << setw(18) << v.x << setw(18) << v.y << setw(18) << v.z << endl;
    }

    template <class T>
    void printv31_d(ofstream &ofs, const string &s, Vector3<T> *u, int n1) const
    {
        ofs << " " << s << " Dimension = " << n1 << endl;
        if (n1 <= 0)return;
        if (u == 0) return;
        for (int i = 0;i < n1;i++)
        {
            printV3(ofs, u[i]);
        }
    }

    template <class T>
    void printv31_d(const string &s, Vector3<T> *u, int n1) const
    {
        cout << "\n " << s << "  dimension = " << n1;
        if (n1 <= 0)return;
        if (u == 0) return;
        for (int i = 0;i < n1;i++)
        {
            printV3(u[i]);
        }
    }

};

#endif // OUTPUT_H

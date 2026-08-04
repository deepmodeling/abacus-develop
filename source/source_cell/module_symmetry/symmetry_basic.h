/**
 * @file symmetry_basic.h
 * @author Zhengpan, mohan, spshu
 * @date 2007-9
 * @brief Basic symmetry operations class.
 */
#ifndef SYMMETRY_BASIC_H
#define SYMMETRY_BASIC_H
#include "symm_other.h"
#include "source_base/mymath.h"
#include "source_base/ylm.h"
#include "source_base/matrix3.h"
namespace ModuleSymmetry
{
/**
 * @brief Basic symmetry operations class.
 */
class Symmetry_Basic
{
    public:

        Symmetry_Basic() {};
        ~Symmetry_Basic() {};

        double epsilon;        ///< the precision of symmetry operation
        double epsilon_input;   ///< the input value of symmetry_prec, should not be changed

        /**
         * @brief Control accuracy - check if two doubles are equal.
         *
         * @param m first value
         * @param n second value
         * @return true if equal within epsilon
         */
        bool equal(const double &m, const double &n)const;

        /**
         * @brief Check boundary and wrap value to [0,1).
         *
         * @param x value to check
         */
        void check_boundary(double &x)const;

        /**
         * @brief Get translation vector.
         *
         * @param x1 first coordinate
         * @param x2 second coordinate
         * @return translation vector
         */
        double get_translation_vector(const double& x1, const double& x2)const;

        /**
         * @brief Check translation.
         *
         * @param x coordinate
         * @param t translation
         */
        void check_translation(double &x, const double &t) const;

        /**
         * @brief Check difference between two values.
         *
         * @param x1 first value
         * @param x2 second value
         * @return difference
         */
        double check_diff(const double& x1, const double& x2) const;
        
        /**
         * @brief Convert vectors from one basis to another.
         *
         * @param va output vectors
         * @param vb input vectors
         * @param num number of vectors
         * @param aa1 first vector of source basis
         * @param aa2 second vector of source basis
         * @param aa3 third vector of source basis
         * @param bb1 first vector of target basis
         * @param bb2 second vector of target basis
         * @param bb3 third vector of target basis
         */
        void veccon(
                double *va,
                double *vb,
                const int num,
                const ModuleBase::Vector3<double> &aa1,
                const ModuleBase::Vector3<double> &aa2,
                const ModuleBase::Vector3<double> &aa3,
                const ModuleBase::Vector3<double> &bb1,
                const ModuleBase::Vector3<double> &bb2,
                const ModuleBase::Vector3<double> &bb3
                );

        /**
         * @brief Generate symmetry operations from generators.
         *
         * @param symgen generators
         * @param ngen number of generators
         * @param symop output symmetry operations
         * @param nop number of operations
         */
        void matrigen(ModuleBase::Matrix3 *symgen, const int ngen, ModuleBase::Matrix3* symop, int &nop) const;

        /**
         * @brief Set up symmetry group.
         *
         * @param symop symmetry operations
         * @param nop number of operations
         * @param ibrav Bravais lattice type
         * @param cal_symm_repr control for symmetry representation output
         */
        void setgroup(ModuleBase::Matrix3 *symop, int &nop, const int &ibrav, 
                      const int* cal_symm_repr) const;

        /**
         * @brief Rotate symmetry operation.
         *
         * @param gmatrix rotation matrix
         * @param gtrans translation vector
         * @param i index in x
         * @param j index in y
         * @param k index in z
         */
        void rotate(
                ModuleBase::Matrix3 &gmatrix, ModuleBase::Vector3<double> &gtrans, 
                int i, int j, int k, const int, const int, const int, int&, int&, int&);

        /**
         * @brief Test atom ordering.
         *
         * @param posi atom positions
         * @param natom number of atoms
         * @param subindex subindex array
         */
        void test_atom_ordering(double *posi, const int natom, int *subindex) const;

        /**
         * @brief Find out the greatest subgroup according to the number of operations of certain type.
         *
         * Used to deal with incomplete group due to a subtle symmetry_prec.
         *
         * @param nrot number of rotations
         * @param ninv number of inversions
         * @param nc2 number of C2 operations
         * @param nc3 number of C3 operations
         * @param nc4 number of C4 operations
         * @param nc6 number of C6 operations
         * @param ns1 number of S1 operations
         * @param ns3 number of S3 operations
         * @param ns4 number of S4 operations
         * @param ns6 number of S6 operations
         * @return subgroup index
         */
        int subgroup(const int& nrot, const int& ninv, const int& nc2, const int& nc3, const int& nc4, const int& nc6,
            const int& ns1, const int& ns3, const int& ns4, const int& ns6)const;

        /**
         * @brief Determine point group.
         *
         * @param nrot number of rotations
         * @param pgnumber point group number
         * @param pgname point group name
         * @param gmatrix rotation matrices
         * @param ofs_running output file stream
         * @param cal_symm_repr control for symmetry representation output
         * @return true if successful
         */
        bool pointgroup(const int& nrot, int& pgnumber, std::string& pgname, const ModuleBase::Matrix3* gmatrix, std::ofstream& ofs_running,
                        const int* cal_symm_repr)const;

protected:
    /**
     * @brief Get Bravais lattice name.
     *
     * @param ibrav Bravais lattice type
     * @return lattice name
     */
    std::string get_brav_name(const int ibrav) const;

    /**
     * @brief Order atoms.
     *
     * @param posi atom positions
     * @param natom number of atoms
     * @param subindex subindex array
     */
    void atom_ordering(double *posi, const int natom, int *subindex);

    /**
     * @brief Order atoms (new version).
     *
     * @param posi atom positions
     * @param natom number of atoms
     * @param subindex subindex array
     */
    void atom_ordering_new(double *posi, const int natom, int *subindex) const;

private:
    /**
     * @brief Order atoms according to index.
     *
     * @param pos atom positions
     * @param nat number of atoms
     * @param index index array
     */
    void order_atoms(double* pos, const int &nat, const int *index) const;

    /**
     * @brief Order atoms along y direction.
     *
     * @param pos atom positions
     * @param oldpos old position
     * @param newpos new position
     */
    void order_y(double *pos, const int &oldpos, const int &newpos);

    /**
     * @brief Order atoms along z direction.
     *
     * @param pos atom positions
     * @param oldpos old position
     * @param newpos new position
     */
    void order_z(double *pos, const int &oldpos, const int &newpos);
};

/// @brief for test only
extern bool test_brav;

}//end of define namespace

#endif

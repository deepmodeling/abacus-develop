/**
 * @file sltk_atom.h
 * @brief FAtom class for storing atom information in neighbor search.
 */
#ifndef INCLUDE_FATOM
#define INCLUDE_FATOM

#include <memory>
#include <vector>

/**
 * @brief A class containing atom position, type and index.
 */
class FAtom
{
public:
    double x;           ///< x coordinate
    double y;           ///< y coordinate
    double z;           ///< z coordinate

    int type;           ///< atom type
    int natom;          ///< atom index

    int cell_x;         ///< cell index in x direction
    int cell_y;         ///< cell index in y direction
    int cell_z;         ///< cell index in z direction

    /**
     * @brief Default constructor.
     */
    FAtom();

    /**
     * @brief Constructor with parameters.
     *
     * @param x_in x coordinate
     * @param y_in y coordinate
     * @param z_in z coordinate
     * @param type_in atom type
     * @param natom_in atom index
     * @param cell_x_in cell index in x direction
     * @param cell_y_in cell index in y direction
     * @param cell_z_in cell index in z direction
     */
    FAtom(const double& x_in, const double& y_in, const double& z_in, 
            const int& type_in, const int& natom_in, 
            const int& cell_x_in, const int& cell_y_in, const int& cell_z_in)
    {
        x = x_in;
        y = y_in;
        z = z_in;
        type = type_in;
        natom = natom_in;
        cell_x = cell_x_in;
        cell_y = cell_y_in;
        cell_z = cell_z_in;
    }

    /**
     * @brief Destructor.
     */
    ~FAtom()
    {
    }
};

#endif

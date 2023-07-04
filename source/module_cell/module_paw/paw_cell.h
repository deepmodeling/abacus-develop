// The Paw_Cell class stores PAW-related information in unitcell
// including number of atoms, number of element types, and so on
// it basically serves as an interface between the unitcell class
// and the libpaw library from ABINIT

#ifndef PAW_CELL
#define PAW_CELL

class Paw_Cell
{
    public:

    Paw_Cell();
    ~Paw_Cell();

};

#endif
// The Paw_Atom class stores PAW-related on-site terms
// It also contains subroutines for connecting with libpaw from ABINIT

#ifndef PAW_ATOM
#define PAW_ATOM

#include <vector>

class Paw_Atom
{
    public:

    Paw_Atom();
    ~Paw_Atom();

    void calculate_rhoij(); //calculate <psi|ptilde><ptilde|psi>

    void convert_rhoij(); //convert to format in libpaw

    private:

    std::vector<std::vector<double>> ca(:,:) //coefficients <psi|ptilde>
    std::vector<std::vector<double>> rhoij(:,:) //on-site density matrix
    
    //for libpaw
    integer :: nrhoijsel; //#. nonzero elements
    std::vector<int>    rhoijselect; //index of nonzero elements
    std::vector<double> rhoij_paw;   //stores nonzero elements
};

#endif
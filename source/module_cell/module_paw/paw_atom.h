// The Paw_Atom class stores PAW-related on-site terms
// It also contains subroutines for connecting with libpaw from ABINIT

#ifndef PAW_ATOM
#define PAW_ATOM

#include <vector>
#include <complex>

class Paw_Atom
{
    public:

    Paw_Atom(){};
    ~Paw_Atom(){};

    //set the sizes of data structures
    void init_paw_atom(const int nproj_in);

    //pass <psi|ptilde> from outside and saves it
    void set_ca(std::vector<std::complex<double>> & ca_in, const double weight_in);

    void reset_rhoij(); //set rhoij = 0
    void accumulate_rhoij(); //calculate and accumulate <psi|ptilde><ptilde|psi> from <psi|ptilde>

    // not sure this is gonna be used, but it is nice to have
    // an interface that returns rhoij and rhoijp I suppose
    std::vector<double> get_rhoij(){return rhoij;}
    std::vector<double> get_rhoijp(){return rhoijp;}
    std::vector<int> get_rhoijselect(){return rhoijselect;}
    int get_nrhoijsel(){return nrhoijsel;}

    void convert_rhoij(); //convert to format in libpaw

    private:

    int nproj;

    std::vector<std::complex<double>> ca; //coefficients <psi|ptilde> for a given psi
    std::vector<double> rhoij; //on-site density matrix, upper triangular

    double weight; //weight of current band
    
    //for libpaw
    int nrhoijsel; //#. nonzero elements
    std::vector<int>    rhoijselect; //index of nonzero elements
    std::vector<double> rhoijp; //rhoij in 'packed' format, only non-zero elements stored;
};

#endif
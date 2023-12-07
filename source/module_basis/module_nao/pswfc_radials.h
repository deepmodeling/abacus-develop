#ifndef PSWFC_RADIALS_H_
#define PSWFC_RADIALS_H_

#include "module_basis/module_nao/radial_set.h"
#include <map>
#include <vector>

class PswfcRadials : public RadialSet {
    public:
        PswfcRadials() {};
        PswfcRadials& operator=(const PswfcRadials& rhs);
        PswfcRadials* clone() const { return new PswfcRadials(*this); }
        ~PswfcRadials() {};

        void build(const std::string& file = "", 
                   const int itype = 0,
                   const double screening_coeff = 0.1,
                   std::ofstream* ptr_log = nullptr, 
                   const int rank = 0);

        void read_upf_pswfc(std::ifstream& ifs,               //!< input file stream from orbital file
                            const double screening_coeff,     //!< screening coefficient
                            std::ofstream* ptr_log = nullptr, //!< output file stream for logging
                            const int rank = 0                //!< MPI rank
        );

        bool startswith(std::string word, std::string pattern);
        std::string read_keyword_value(std::ifstream& ifs, std::string word);
        std::string steal_from_quotes(std::string word);
        std::string steal_from_quotes(std::ifstream& ifs, std::string word);
    private:

};
#endif // PSWFC_RADIALS_H_
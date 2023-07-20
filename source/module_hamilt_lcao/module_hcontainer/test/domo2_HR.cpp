
#include "../hcontainer.h"

void add_HR(hamilt::HContainer<double> HR_in)
{
    if(this->HR == nullptr) 
    {
        this->HR = new hamilt::HContainer<double>(HR_in);
    }
    else 
    {
        for(int iap=0;iap<HR_in.size_atom_pairs();iap++)
        {
            auto tmp = HR_in.get_atom_pair(iap);
            this->HR->insert_pair(tmp);
        }
    }
}

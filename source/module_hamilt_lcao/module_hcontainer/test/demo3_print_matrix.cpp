#include "../hcontainer.h"

void print_HR(const HContainer& HR_in, const UnitCell& ucell, const int* R_index, ofstream& ofs)
{
    HR_in.fix_R(R_index[0], R_index[1], R_index[2]);
    ofs<<"R_index: "<<R_index[0]<<" "<<R_index[1]<<" "<<R_index[2]<<"\n";
    for(int iap=0;iap<HR_in.size_atom_pairs();iap++)
    {
        auto tmp = HR_in.get_atom_pair(iap);
        int iat1 = tmp.get_atom_i();
        int iat2 = tmp.get_atom_j();
        int T1 = ucell.get_element(iat1);
        int T2 = ucell.get_element(iat2);
        int I1 = ucell.get_element_atom(iat1);
        int I2 = ucell.get_element_atom(iat2);
        const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
        const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
        int size1 = GlobalC::ucell.get_nw(T1);
        int size2 = GlobalC::ucell.get_nw(T2);
        for(int iw1=0;iw1<size1;++iw1)
        {
            const int global_index1 = start1 + iw1;
            for(int iw2=0;iw2<size2;++iw2)
            {
                const int global_index2 = start2 + iw2;
                ofs << global_index1<<" "<<global_index2<<" "<< tmp.get_matrix_value(start1 + iw1, start2 + iw2) << '\n';
            }
        }
    }
    HR_in.unfix_R();
}
#ifdef __EXX
#include "op_exx_lcao.h"
#include "source_base/blacs_connector.h"

namespace hamilt
{
    RI::Cell_Nearest<int, int, 3, double, 3> init_cell_nearest(const UnitCell& ucell, const std::array<int, 3>& Rs_period)
    {
        RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
        std::map<int, std::array<double, 3>> atoms_pos;
        for (int iat = 0; iat < ucell.nat; ++iat) {
            atoms_pos[iat] = RI_Util::Vector3_to_array3(
                ucell.atoms[ucell.iat2it[iat]]
                .tau[ucell.iat2ia[iat]]);
        }
        const std::array<std::array<double, 3>, 3> latvec
            = { RI_Util::Vector3_to_array3(ucell.a1),
               RI_Util::Vector3_to_array3(ucell.a2),
               RI_Util::Vector3_to_array3(ucell.a3) };
        cell_nearest.init(atoms_pos, latvec, Rs_period);
        return cell_nearest;
    }

template<>
void OperatorEXX<OperatorLCAO<double, double>>::add_loaded_Hexx(const int ik)
{
    BlasConnector::axpy(this->hR->get_paraV()->get_local_size(), 1.0, this->Hexxd_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, double>>::add_loaded_Hexx(const int ik)
{
    BlasConnector::axpy(this->hR->get_paraV()->get_local_size(), 1.0, this->Hexxc_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}
template<>
void OperatorEXX<OperatorLCAO<std::complex<double>, std::complex<double>>>::add_loaded_Hexx(const int ik)
{
    BlasConnector::axpy(this->hR->get_paraV()->get_local_size(), 1.0, this->Hexxc_k_load[ik].data(), 1, this->hsk->get_hk(), 1);
}

} // namespace hamilt
#endif
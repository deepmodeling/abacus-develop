#include "./dm_band.h"
#include "source_lcao/module_ri/ri_util.h"
namespace LR
{
    // using TC = std::array<int, 3>;
    // using TAC = std::pair<int, TC>;
    // template <typename T>
    // using TDM = std::map<int, std::map<TAC, RI::Tensor<T>>>;
    template<>
    void DMBand<double>::cal_dm_band(const int iband1, const int iband2, const int ik,
        DMBand<double>::TDM& dm_band, const double fac,
        const std::vector<int> nws1, const std::vector<int> nws2) const
    {
        ModuleBase::TITLE("DMBand", "cal_dm_band");
        // NOTICE: DM_onebase will be passed into `cal_energy` interface and conjugated by "zdotc". 
        // So the formula should be the same as RHS. instead of LHS of the A-matrix, 
        // i.e. c1v · conj(c2o) ·  e^{-ik(R2-R1)}
        assert(ik == 0);
        const bool use_nws1 = !nws1.empty();
        const bool use_nws2 = !nws2.empty();
        for (auto cell : bvk_cells_)
        {
            for (int it1 = 0;it1 < ucell_.ntype;++it1)
                for (int ia1 = 0; ia1 < ucell_.atoms[it1].na; ++ia1)
                    for (int it2 = 0;it2 < ucell_.ntype;++it2)
                        for (int ia2 = 0;ia2 < ucell_.atoms[it2].na;++ia2)
                        {
                            const int iat1 = ucell_.itia2iat(it1, ia1);
                            const int iat2 = ucell_.itia2iat(it2, ia2);
                            const std::size_t nw1 = ucell_.atoms[it1].nw;
                            const std::size_t nw2 = ucell_.atoms[it2].nw;
                            RI::Tensor<double> dm_tmp({ nw1, nw2 });
                            for (int iw1 = 0;iw1 < nw1;++iw1)
                                for (int iw2 = 0;iw2 < nw2;++iw2)
                                {
                                    const int iwt1 = use_nws1 ? nws1[it1] : ucell_.itiaiw2iwt(it1, ia1, iw1);
                                    const int iwt2 = use_nws2 ? nws2[it2] : ucell_.itiaiw2iwt(it2, ia2, iw2);
                                    if (pmat_.in_this_processor(iwt1, iwt2))
                                        dm_tmp(iw1, iw2) = fac * c1_(ik, iband1, iwt1) * c2_(ik, iband2, iwt2);
                                }
                            dm_band[iat1][std::make_pair(iat2, cell)] = dm_tmp;
                        }
        }
    }

    template<>
    void DMBand<std::complex<double>>::cal_dm_band(const int iband1, const int iband2, const int ik,
        DMBand<std::complex<double>>::TDM& dm_band, const std::complex<double> fac,
        const std::vector<int> nws1, const std::vector<int> nws2) const
    {
        ModuleBase::TITLE("DMBand", "cal_dm_band");
        // NOTICE: dm_band will be passed into `cal_energy` interface and conjugated by "zdotc". 
        // So the formula should be the same as RHS. instead of LHS of the A-matrix, 
        // i.e. c1v · conj(c2o) ·  e^{-ik(R2-R1)}
        const bool use_nws1 = !nws1.empty();
        const bool use_nws2 = !nws2.empty();
        for (auto cell : bvk_cells_)
        {
            std::complex<double> fac_phase = RI::Global_Func::convert<std::complex<double>>(std::exp(
                -ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT * (kvec_c_.at(ik) * (RI_Util::array3_to_Vector3(cell) * ucell_.latvec)))) * fac;
            for (int it1 = 0;it1 < ucell_.ntype;++it1)
                for (int ia1 = 0; ia1 < ucell_.atoms[it1].na; ++ia1)
                    for (int it2 = 0;it2 < ucell_.ntype;++it2)
                        for (int ia2 = 0;ia2 < ucell_.atoms[it2].na;++ia2)
                        {
                            const int iat1 = ucell_.itia2iat(it1, ia1);
                            const int iat2 = ucell_.itia2iat(it2, ia2);
                            const std::size_t nw1 = ucell_.atoms[it1].nw;
                            const std::size_t nw2 = ucell_.atoms[it2].nw;
                            RI::Tensor<std::complex<double>> dm_tmp({ nw1, nw2 });
                            for (int iw1 = 0;iw1 < nw1;++iw1)
                                for (int iw2 = 0;iw2 < nw2;++iw2)
                                {
                                    const int iwt1 = use_nws1 ? nws1[it1] : ucell_.itiaiw2iwt(it1, ia1, iw1);
                                    const int iwt2 = use_nws2 ? nws2[it2] : ucell_.itiaiw2iwt(it2, ia2, iw2);
                                    if (pmat_.in_this_processor(iwt1, iwt2))
                                        dm_tmp(iw1, iw2) = fac_phase * c1_(ik, iband1, iwt1) * std::conj(c2_(ik, iband2, iwt2));
                                }
                            dm_band[iat1][std::make_pair(iat2, cell)] = dm_tmp;
                        }
        }
    }
    
}
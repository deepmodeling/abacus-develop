
#include "source_cell/unitcell.h"
#include "source_base/parallel_2d.h"
#include "source_psi/psi.h"
#include <RI/global/Tensor.h>
namespace LR
{
    template<typename T>
    class DMBand
    {
        using TC = std::array<int, 3>;
        using TAC = std::pair<int, TC>;
        using TDM = std::map<int, std::map<TAC, RI::Tensor<T>>>;
    public:
        DMBand(const UnitCell& ucell,
            const Parallel_2D& pmat,
            const std::vector<ModuleBase::Vector3<double>>& kvec_c,
            const std::vector<TC>& bvk_cells,
            const psi::Psi<T>& c1, const psi::Psi<T>& c2)
            : ucell_(ucell), pmat_(pmat), 
             kvec_c_(kvec_c), bvk_cells_(bvk_cells), c1_(c1), c2_(c2) {};
        DMBand() = delete;
        ~DMBand() = default;

        void cal_dm_band(const int iband1, const int iband2, const int ik, TDM& dm_band, const T fac = 1.0,
            const std::vector<int> nws1 = {}, const std::vector<int> nws2 = {}) const;
        void eval(const int iband1, const int iband2, const int ik, const T fac = 1.0)
        {
            cal_dm_band(iband1, iband2, ik, data_, fac);
        }

        DMBand<T> operator+(const DMBand<T>& rhs) const
        {
            DMBand<T> res = *this;
            for (auto& ia1_map1 : rhs.data_)
            {
                int iat1 = ia1_map1.first;
                for (auto& ia2_cell_tensor : ia1_map1.second)
                {
                    const auto& iat2_cell = ia2_cell_tensor.first;
                        const RI::Tensor<T>&tensor = ia2_cell_tensor.second;
                    res.data_[iat1][iat2_cell] += tensor;
                }
            }
            return res;
        }
        DMBand<T>& operator+=(const DMBand<T>& rhs)
        {
            for (auto ia1_map1 = rhs.data_.begin(); ia1_map1 != rhs.data_.end(); ++ia1_map1)
            {
                int iat1 = ia1_map1->first;
                for (auto ia2_cell_tensor = ia1_map1->second.begin(); ia2_cell_tensor != ia1_map1->second.end(); ++ia2_cell_tensor)
                {
                    const auto& iat2_cell = ia2_cell_tensor->first;
                    const RI::Tensor<T>& tensor = ia2_cell_tensor->second;
                    this->data_[iat1][iat2_cell] += tensor;
                }
            }
            return *this;
        }
        DMBand<T> operator-(const DMBand<T>& rhs) const
        {
            DMBand<T> res = *this;
            for (auto ia1_map1 = rhs.data_.begin(); ia1_map1 != rhs.data_.end(); ++ia1_map1)
            {
                int iat1 = ia1_map1->first;
                for (auto ia2_cell_tensor = ia1_map1->second.begin(); ia2_cell_tensor != ia1_map1->second.end(); ++ia2_cell_tensor)
                {
                    const auto& iat2_cell = ia2_cell_tensor->first;
                    const RI::Tensor<T>& tensor = ia2_cell_tensor->second;
                    res.data_[iat1][iat2_cell] -= tensor;
                }
            }
            return res;
        }
        DMBand<T>& operator-=(const DMBand<T>& rhs)
        {
            for (auto ia1_map1 = rhs.data_.begin(); ia1_map1 != rhs.data_.end(); ++ia1_map1)
            {
                int iat1 = ia1_map1->first;
                for (auto ia2_cell_tensor = ia1_map1->second.begin(); ia2_cell_tensor != ia1_map1->second.end(); ++ia2_cell_tensor)
                {
                    const auto& iat2_cell = ia2_cell_tensor->first;
                    const RI::Tensor<T>& tensor = ia2_cell_tensor->second;
                    this->data_[iat1][iat2_cell] -= tensor;
                }
            }
            return *this;
        }

        const TDM& get_data() const { return data_; }
        
    private:
        const UnitCell& ucell_;
        const std::vector<TC> bvk_cells_;
        const std::vector<ModuleBase::Vector3<double>>& kvec_c_;
        const Parallel_2D& pmat_;
        const psi::Psi<T>& c1_; // band 1 (global)
        const psi::Psi<T>& c2_; // band 2 (global)
        TDM data_;
    };
}
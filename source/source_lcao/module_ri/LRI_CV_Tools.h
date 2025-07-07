//=======================
// AUTHOR : Peize Lin
// DATE :   2022-10-24
//=======================

#ifndef LRI_CV_TOOLS_H
#define LRI_CV_TOOLS_H

#include "source_base/abfs-vector3_order.h"
#include "source_lcao/module_ri/abfs.h"

#include <RI/global/Array_Operator.h>
#include <RI/global/Tensor.h>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace LRI_CV_Tools
{
template <typename Tdata>
extern RI::Tensor<Tdata> cal_I(const RI::Tensor<Tdata>& m);
template <typename Tdata>
extern std::vector<std::vector<RI::Tensor<Tdata>>> cal_I(const std::vector<std::vector<RI::Tensor<Tdata>>>& ms);

template <typename Tdata>
inline RI::Tensor<Tdata> transform_Rm(const RI::Tensor<Tdata>& V);
template <typename Tdata>
inline std::array<RI::Tensor<Tdata>, 3> transform_Rm(const std::array<RI::Tensor<Tdata>, 3>& dV);

// template<typename T> inline bool exist(const T &V);

// template<typename T1, typename T2, typename Treturn>
// extern Treturn mul1(const T1 &t1, const T2 &t2);
// template<typename T1, typename T2, typename Treturn>
// extern Treturn mul2(const T1 &mat, const T2 &vec);

template <typename Tdata>
inline bool exist(const RI::Tensor<Tdata>& V);
template <typename T, std::size_t N>
inline bool exist(const std::array<T, N>& dV);

template <typename Tdata>
extern RI::Tensor<Tdata> mul1(const RI::Tensor<Tdata>& t1, const RI::Tensor<Tdata>& t2);
template <typename T>
extern std::array<T, 3> mul1(const std::array<T, 3>& t1, const T& t2);

template <typename Tdata>
extern std::vector<RI::Tensor<Tdata>> mul2(const std::vector<std::vector<RI::Tensor<Tdata>>>& mat,
                                           const std::vector<RI::Tensor<Tdata>>& vec);
template <typename T1, typename T2>
extern std::array<T2, 3> mul2(const T1& t1, const std::array<T2, 3>& t2);
template <typename T>
extern RI::Tensor<T> mul2(const T& t1, const RI::Tensor<T>& t2);
template <typename T, typename TkeyA, typename TkeyB, typename Tvalue>
extern std::map<TkeyA, std::map<TkeyB, Tvalue>> mul2(const T& t1, const std::map<TkeyA, std::map<TkeyB, Tvalue>>& t2);

// template<typename T, std::size_t N>
// std::array<T,N> operator-(const std::array<T,N> &v1, const std::array<T,N>
// &v2); template<typename T> std::vector<T> operator-(const std::vector<T> &v1,
// const std::vector<T> &v2);
template <typename T, std::size_t N>
extern std::vector<std::array<T, N>> minus(const std::vector<std::array<T, N>>& v1,
                                           const std::vector<std::array<T, N>>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
extern std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> minus(
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v1,
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
inline std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> minus(
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v1,
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue>
extern std::map<TkeyA, std::map<TkeyB, Tvalue>> minus(std::map<TkeyA, std::map<TkeyB, Tvalue>>& v1,
                                                      std::map<TkeyA, std::map<TkeyB, Tvalue>>& v2);

template <typename T, std::size_t N>
extern std::vector<std::array<T, N>> add(const std::vector<std::array<T, N>>& v1,
                                         const std::vector<std::array<T, N>>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
extern std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> add(
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v1,
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
inline std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> add(
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v1,
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v2);
template <typename TkeyA, typename TkeyB, typename Tvalue>
extern std::map<TkeyA, std::map<TkeyB, Tvalue>> add(std::map<TkeyA, std::map<TkeyB, Tvalue>>& v1,
                                                    std::map<TkeyA, std::map<TkeyB, Tvalue>>& v2);

template <typename T, std::size_t N>
extern std::array<T, N> negative(const std::array<T, N>& v_in);

// template<typename T> T transpose12(const T &c_in);
template <typename Tdata>
RI::Tensor<Tdata> transpose12(const RI::Tensor<Tdata>& c_in);
template <typename T, std::size_t N>
std::array<T, N> transpose12(const std::array<T, N>& c_in);

template <typename T, std::size_t N>
extern std::array<std::vector<T>, N> change_order(std::vector<std::array<T, N>>&& ds_in);
template <typename T, std::size_t N>
std::vector<std::array<T, N>> change_order(std::array<std::vector<T>, N>&& ds_in);
template <typename T, std::size_t N>
extern std::array<std::vector<std::vector<T>>, N> change_order(std::vector<std::vector<std::array<T, N>>>&& ds_in);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
extern std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> change_order(
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>&& ds_in);
template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
extern std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> change_order(
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>&& ds_in);

template <typename Tcell>
extern std::array<Tcell, 3> cal_latvec_range(const double& rcut_times,
                                             const UnitCell& ucell,
                                             const std::vector<double>& orb_cutoff);

template <typename TA, typename Tcell, typename Tdata>
extern std::map<int, std::map<int, std::map<Abfs::Vector3_Order<double>, RI::Tensor<Tdata>>>> get_CVws(
    const UnitCell& ucell,
    const std::map<TA, std::map<std::pair<TA, std::array<Tcell, 3>>, RI::Tensor<Tdata>>>& CVs);
template <typename TA, typename Tcell, typename Tdata>
extern std::map<int, std::map<int, std::map<Abfs::Vector3_Order<double>, std::array<RI::Tensor<Tdata>, 3>>>> get_dCVws(
    const UnitCell& ucell,
    const std::map<TA, std::map<std::pair<TA, std::array<Tcell, 3>>, std::array<RI::Tensor<Tdata>, 3>>>& dCVs);
template <typename TA, typename TC, typename Tdata>
extern std::array<std::array<std::map<TA, std::map<std::pair<TA, TC>, RI::Tensor<Tdata>>>, 3>, 3> cal_dMRs(
    const UnitCell& ucell,
    const std::array<std::map<TA, std::map<std::pair<TA, TC>, RI::Tensor<Tdata>>>, 3>& dMs);

using TC = std::array<int, 3>;
using TAC = std::pair<int, TC>;
template <typename T>
using TLRI = std::map<int, std::map<TAC, RI::Tensor<T>>>;
template <typename T>
TLRI<T> read_Cs_ao(const std::string& file_path, const double& threshold = 1e-10);
template <typename T>
void write_Cs_ao(const TLRI<T>& Vs, const std::string& file_path);
template <typename T>
TLRI<T> read_Vs_abf(const std::string& file_path, const double& threshold = 1e-10);
template <typename T>
void write_Vs_abf(const TLRI<T>& Vs, const std::string& file_path);

template <typename T>
struct is_std_array : std::false_type
{
};
template <typename T, std::size_t N>
struct is_std_array<std::array<T, N>> : std::true_type
{
};
template <typename T>
struct is_tensor : std::false_type
{
};
template <typename T>
struct is_tensor<RI::Tensor<T>> : std::true_type
{
};

template <typename Tout>
struct TinType;

template <typename T>
struct TinType<RI::Tensor<T>>
{
    using type = T;
};

template <typename T, std::size_t N>
struct TinType<std::array<RI::Tensor<T>, N>>
{
    using type = T;
};

template <typename Tdata, typename = std::enable_if_t<!is_std_array<Tdata>::value>>
inline void init_elem(Tdata& data, const size_t ndim0, const size_t ndim1)
{
    data = Tdata({ndim0, ndim1});
};
template <typename T, std::size_t N>
extern void init_elem(std::array<RI::Tensor<T>, N>& data, const size_t ndim0, const size_t ndim1);

template <typename Tdata, typename = std::enable_if_t<!is_std_array<Tdata>::value && !is_tensor<Tdata>::value>>
inline void add_elem(Tdata& data, const Tdata& val, const Tdata& frac)
{
    data += frac * val;
};
template <typename T, std::size_t N>
extern void add_elem(std::array<T, N>& data, const T& val, const T& frac);
template <typename Tdata, typename = std::enable_if_t<is_tensor<Tdata>::value>>
inline void add_elem(const Tdata& data,
                     const int lmp,
                     const int lmq,
                     const typename TinType<Tdata>::type& val,
                     const typename TinType<Tdata>::type& frac)
{
    data(lmp, lmq) += frac * val;
};
template <typename T, std::size_t N>
extern void add_elem(std::array<RI::Tensor<T>, N>& data,
                     const int lmp,
                     const int lmq,
                     const std::array<T, N>& val,
                     const T& frac);
template <typename Tdata, typename = std::enable_if_t<is_tensor<Tdata>::value>>
inline void add_elem(Tdata& data,
                     const int lmp0,
                     const int lmq0,
                     const Tdata& val,
                     const int lmp1,
                     const int lmq1,
                     const typename TinType<Tdata>::type& frac)
{
    data(lmp0, lmq0) += frac * val(lmp1, lmq1);
};
template <typename T, std::size_t N>
extern void add_elem(std::array<RI::Tensor<T>, N>& data,
                     const int lmp0,
                     const int lmq0,
                     const std::array<RI::Tensor<T>, N>& val,
                     const int lmp1,
                     const int lmq1,
                     const T& frac);

template <typename Tout, typename Tin>
inline RI::Tensor<Tout> convert(RI::Tensor<Tin>&& data);
template <typename Tout, typename Tin, std::size_t N>
extern std::array<RI::Tensor<Tout>, N> convert(std::array<RI::Tensor<Tin>, N>&& data);

// template <typename T>
// typename std::enable_if<!RI::Global_Func::is_complex<T>::value, T>::type
// inline check_zero(T value) {
//     return (std::abs(value) < 1e-8) ? static_cast<T>(0) : value;
// }

// template <typename T>
// typename std::enable_if<RI::Global_Func::is_complex<T>::value, T>::type
// inline check_zero(const T& value) {
//     using RealType = typename T::value_type;
//     RealType real_part = std::real(value);
//     RealType imag_part = std::imag(value);

//     real_part = (std::abs(real_part) < 1e-8) ? 0 : real_part;
//     imag_part = (std::abs(imag_part) < 1e-8) ? 0 : imag_part;

//     return std::complex<RealType>(real_part, imag_part);
// }

// template <typename T>
// extern RI::Tensor<T> check_zero(RI::Tensor<T>&& data);
// template <typename T, std::size_t N>
// extern std::array<RI::Tensor<T>, N> check_zero(std::array<RI::Tensor<T>, N>&& data);

template <typename T>
struct plus
{
    T operator()(const T& lhs, const T& rhs) const
    {
        using namespace RI::Array_Operator;
        return lhs + rhs;
    }
};
} // namespace LRI_CV_Tools

#include "LRI_CV_Tools.hpp"
#include "write_ri_cv.hpp"
#endif

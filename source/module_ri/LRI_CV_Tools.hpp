//=======================
// AUTHOR : Peize Lin
// DATE :   2022-10-24
//=======================

#ifndef LRI_CV_TOOLS_HPP
#define LRI_CV_TOOLS_HPP


#include "../source_base/mathzone.h"
#include "Inverse_Matrix.h"
#include "LRI_CV_Tools.h"
#include "RI_Util.h"

#include <RI/global/Global_Func-1.h>
#include <RI/global/Map_Operator.h>
#include "../module_hamilt_pw/hamilt_pwdft/global.h"

template <typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::cal_I(const RI::Tensor<Tdata>& m) {
    Inverse_Matrix<Tdata> I;
    I.input(m);
    I.cal_inverse(Inverse_Matrix<Tdata>::Method::potrf);
    return I.output();
}

template <typename Tdata>
std::vector<std::vector<RI::Tensor<Tdata>>>
    LRI_CV_Tools::cal_I(const std::vector<std::vector<RI::Tensor<Tdata>>>& ms) {
    Inverse_Matrix<Tdata> I;
    I.input(ms);
    I.cal_inverse(Inverse_Matrix<Tdata>::Method::potrf);
    return I.output({ms[0][0].shape[0], ms[1][0].shape[0]},
                    {ms[0][0].shape[1], ms[0][1].shape[1]});
}

template <typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::transform_Rm(const RI::Tensor<Tdata>& V) {
    return V.transpose();
}

template <typename Tdata>
std::array<RI::Tensor<Tdata>, 3>
    LRI_CV_Tools::transform_Rm(const std::array<RI::Tensor<Tdata>, 3>& dV) {
    return std::array<RI::Tensor<Tdata>, 3>{-dV[0].transpose(),
                                            -dV[1].transpose(),
                                            -dV[2].transpose()};
}

template <typename Tdata>
bool LRI_CV_Tools::exist(const RI::Tensor<Tdata>& V) {
    return !V.empty();
}

template <typename T, std::size_t N>
bool LRI_CV_Tools::exist(const std::array<T, N>& dV) {
    for (size_t i = 0; i < 3; ++i)
        if (!dV[i].empty())
            return true;
    return false;
}

template <typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::mul1(const RI::Tensor<Tdata>& t1,
                                     const RI::Tensor<Tdata>& t2) {
    const size_t sa0 = t1.shape[0], sa1 = t2.shape[0], sl0 = t2.shape[1],
                 sl1 = t2.shape[2];
    return (t1 * t2.reshape({sa1, sl0 * sl1})).reshape({sa0, sl0, sl1});
}
template <typename T>
std::array<T, 3> LRI_CV_Tools::mul1(const std::array<T, 3>& t1, const T& t2) {
    return std::array<T, 3>{mul1(t1[0], t2), mul1(t1[1], t2), mul1(t1[2], t2)};
}
/*
template<typename T>
std::array<T,3> LRI_CV_Tools::mul1(
    const T &t1,
    const std::array<T,3> &t2)
{
    return std::array<T,3>{
        mul1(t1,t2[0]), mul1(t1,t2[1]), mul1(t1,t2[2]) };
}
*/

template <typename Tdata>
std::vector<RI::Tensor<Tdata>>
    LRI_CV_Tools::mul2(const std::vector<std::vector<RI::Tensor<Tdata>>>& mat,
                       const std::vector<RI::Tensor<Tdata>>& vec) {
    const size_t sa0 = vec[0].shape[0], sa1 = vec[1].shape[0],
                 sl0 = vec[0].shape[1], sl1 = vec[0].shape[2];
    const RI::Tensor<Tdata> vec0 = vec[0].reshape({sa0, sl0 * sl1}),
                            vec1 = vec[1].reshape({sa1, sl0 * sl1});
    return std::vector<RI::Tensor<Tdata>>{
        (mat[0][0] * vec0 + mat[0][1] * vec1).reshape({sa0, sl0, sl1}),
        (mat[1][0] * vec0 + mat[1][1] * vec1).reshape({sa1, sl0, sl1})};
}
/*
template<typename T1, typename T2>
std::array<T2,3> LRI_CV_Tools::mul2(
    const std::array<T1,3> &t1,
    const T2 &t2)
{
    return std::array<T2,3>{
        mul2(t1[0],t2), mul2(t1[1],t2), mul2(t1[2],t2) };
}
*/
template <typename T1, typename T2>
std::array<T2, 3> LRI_CV_Tools::mul2(const T1& t1,
                                     const std::array<T2, 3>& t2) {
    return std::array<T2, 3>{mul2(t1, t2[0]), mul2(t1, t2[1]), mul2(t1, t2[2])};
}

template <typename T>
RI::Tensor<T> LRI_CV_Tools::mul2(const T& t1, const RI::Tensor<T>& t2) {
    return t1 * t2;
}

template <typename T, typename TkeyA, typename TkeyB, typename Tvalue>
std::map<TkeyA, std::map<TkeyB, Tvalue>>
    LRI_CV_Tools::mul2(const T& t1,
                       const std::map<TkeyA, std::map<TkeyB, Tvalue>>& t2) {
    std::map<TkeyA, std::map<TkeyB, Tvalue>> res;
    for (const auto& outerPair: t2) {
        const TkeyA keyA = outerPair.first;
        const std::map<TkeyB, Tvalue>& innerMap = outerPair.second;
        std::map<TkeyB, Tvalue> newInnerMap;

        for (const auto& innerPair: innerMap) {
            const TkeyB keyB = innerPair.first;
            const Tvalue value = innerPair.second;
            newInnerMap[keyB] = mul2(t1, value);
        }

        res[keyA] = newInnerMap;
    }

    return res;
}

/*
template<typename T, std::size_t N>
std::array<T,N> LRI_CV_Tools::operator-(const std::array<T,N> &v1, const
std::array<T,N> &v2)
{
    std::array<T,N> v;
    for(std::size_t i=0; i<N; ++i)
        v[i] = v1[i] - v2[i];
    return v;
}
template<typename T>
std::vector<T> LRI_CV_Tools::operator-(const std::vector<T> &v1, const
std::vector<T> &v2)
{
    assert(v1.size()==v2.size());
    std::vector<T> v(v1.size());
    for(std::size_t i=0; i<v.size(); ++i)
        v[i] = v1[i] - v2[i];
    return v;
}
*/
template <typename T, std::size_t N>
std::vector<std::array<T, N>>
    LRI_CV_Tools::minus(const std::vector<std::array<T, N>>& v1,
                        const std::vector<std::array<T, N>>& v2) {
    assert(v1.size() == v2.size());
    std::vector<std::array<T, N>> v(v1.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < N; ++j)
            v[i][j] = v1[i][j] - v2[i][j];
    return v;
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> LRI_CV_Tools::minus(
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v1,
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v2) {
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> v1_order
        = change_order(std::move(v1));
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> v2_order
        = change_order(std::move(v2));
    auto dv = minus(v1_order, v2_order);
    return change_order(std::move(dv));
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> LRI_CV_Tools::minus(
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v1,
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v2) {
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> dv;
    for (size_t i = 0; i != N; ++i)
        dv[i] = minus(v1[i], v2[i]);
    return dv;
}

template <typename TkeyA, typename TkeyB, typename Tvalue>
std::map<TkeyA, std::map<TkeyB, Tvalue>>
    LRI_CV_Tools::minus(std::map<TkeyA, std::map<TkeyB, Tvalue>>& v1,
                        std::map<TkeyA, std::map<TkeyB, Tvalue>>& v2) {
    assert(v1.size() == v2.size());
    using namespace RI::Map_Operator;
    using namespace RI::Array_Operator;

    std::map<TkeyA, std::map<TkeyB, Tvalue>> dv;
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it1 != v1.end() && it2 != v2.end()) {
        assert(it1->first == it2->first);
        const TkeyA& keyA = it1->first;
        const std::map<TkeyB, Tvalue>& map1 = it1->second;
        const std::map<TkeyB, Tvalue>& map2 = it2->second;
        dv[keyA] = map1 - map2;
        ++it1;
        ++it2;
    }
    return dv;
}

template <typename T, std::size_t N>
std::vector<std::array<T, N>>
    LRI_CV_Tools::add(const std::vector<std::array<T, N>>& v1,
                        const std::vector<std::array<T, N>>& v2) {
    assert(v1.size() == v2.size());
    std::vector<std::array<T, N>> v(v1.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        for (std::size_t j = 0; j < N; ++j)
            v[i][j] = v1[i][j] + v2[i][j];
    return v;
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> LRI_CV_Tools::add(
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v1,
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>& v2) {
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> v1_order
        = change_order(std::move(v1));
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> v2_order
        = change_order(std::move(v2));
    auto dv = add(v1_order, v2_order);
    return change_order(std::move(dv));
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> LRI_CV_Tools::add(
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v1,
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>& v2) {
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> dv;
    for (size_t i = 0; i != N; ++i)
        dv[i] = add(v1[i], v2[i]);
    return dv;
}

template <typename TkeyA, typename TkeyB, typename Tvalue>
std::map<TkeyA, std::map<TkeyB, Tvalue>>
    LRI_CV_Tools::add(std::map<TkeyA, std::map<TkeyB, Tvalue>>& v1,
                        std::map<TkeyA, std::map<TkeyB, Tvalue>>& v2) {
    assert(v1.size() == v2.size());
    using namespace RI::Map_Operator;
    using namespace RI::Array_Operator;

    std::map<TkeyA, std::map<TkeyB, Tvalue>> dv;
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    while (it1 != v1.end() && it2 != v2.end()) {
        assert(it1->first == it2->first);
        const TkeyA& keyA = it1->first;
        const std::map<TkeyB, Tvalue>& map1 = it1->second;
        const std::map<TkeyB, Tvalue>& map2 = it2->second;
        dv[keyA] = map1 + map2;
        ++it1;
        ++it2;
    }
    return dv;
}

template <typename T, std::size_t N>
std::array<T, N> LRI_CV_Tools::negative(const std::array<T, N>& v_in) {
    std::array<T, N> v_out;
    for (std::size_t i = 0; i < N; ++i)
        v_out[i] = -v_in[i];
    return v_out;
}

template <typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::transpose12(const RI::Tensor<Tdata>& c_in) {
    RI::Tensor<Tdata> c_out({c_in.shape[0], c_in.shape[2], c_in.shape[1]});
    for (size_t i0 = 0; i0 < c_in.shape[0]; ++i0)
        for (size_t i1 = 0; i1 < c_in.shape[1]; ++i1)
            for (size_t i2 = 0; i2 < c_in.shape[2]; ++i2)
                c_out(i0, i2, i1) = c_in(i0, i1, i2);
    return c_out;
}

template <typename T, std::size_t N>
std::array<T, N> LRI_CV_Tools::transpose12(const std::array<T, N>& c_in) {
    std::array<T, N> c_out;
    for (size_t i = 0; i < N; ++i)
        c_out[i] = transpose12(c_in[i]);
    return c_out;
}

template <typename T, std::size_t N>
std::array<std::vector<T>, N>
    LRI_CV_Tools::change_order(std::vector<std::array<T, N>>&& ds_in) {
    std::array<std::vector<T>, N> ds;
    for (int ix = 0; ix < N; ++ix) {
        ds[ix].resize(ds_in.size());
        for (int iv = 0; iv < ds_in.size(); ++iv)
            ds[ix][iv] = std::move(ds_in[iv][ix]);
    }
    return ds;
}

template <typename T, std::size_t N>
std::vector<std::array<T, N>>
    LRI_CV_Tools::change_order(std::array<std::vector<T>, N>&& ds_in) {
    std::vector<std::array<T, N>> ds(ds_in[0].size());
    for (int ix = 0; ix < N; ++ix) {
        assert(ds.size() == ds_in[ix].size());
        for (int iv = 0; iv < ds.size(); ++iv)
            ds[iv][ix] = std::move(ds_in[ix][iv]);
    }
    return ds;
}

template <typename T, std::size_t N>
std::array<std::vector<std::vector<T>>, N> LRI_CV_Tools::change_order(
    std::vector<std::vector<std::array<T, N>>>&& ds_in) {
    std::array<std::vector<std::vector<T>>, N> ds;
    for (int ix = 0; ix < N; ++ix) {
        ds[ix].resize(ds_in.size());
        for (int i0 = 0; i0 < ds_in.size(); ++i0) {
            ds[ix][i0].resize(ds_in[i0].size());
            for (int i1 = 0; i1 < ds_in[i0].size(); ++i1)
                ds[ix][i0][i1] = std::move(ds_in[i0][i1][ix]);
        }
    }
    return ds;
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>
    LRI_CV_Tools::change_order(
        std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>&& ds_in) {
    std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N> ds;
    for (auto& ds_A: ds_in)
        for (auto& ds_B: ds_A.second)
            for (int ix = 0; ix < N; ++ix)
                ds[ix][ds_A.first][ds_B.first] = std::move(ds_B.second[ix]);
    return ds;
}

template <typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>>
    LRI_CV_Tools::change_order(
        std::array<std::map<TkeyA, std::map<TkeyB, Tvalue>>, N>&& ds_in) {
    std::map<TkeyA, std::map<TkeyB, std::array<Tvalue, N>>> ds;
    for (int ix = 0; ix < N; ++ix)
        for (auto& ds_A: ds_in[ix])
            for (auto& ds_B: ds_A.second)
                ds[ds_A.first][ds_B.first][ix] = std::move(ds_B.second);
    return ds;
}

template <typename Tcell>
std::array<Tcell, 3> LRI_CV_Tools::cal_latvec_range(const double& rcut_times,
							   const UnitCell &ucell,
							   const std::vector<double>& orb_cutoff) {
    double Rcut_max = 0;
    for(int T=0; T<ucell.ntype; ++T)
		Rcut_max = std::max(Rcut_max, orb_cutoff[T]);
	const ModuleBase::Vector3<double> proj = ModuleBase::Mathzone::latvec_projection(
		std::array<ModuleBase::Vector3<double>,3>{ucell.a1, ucell.a2, ucell.a3});
	const ModuleBase::Vector3<double> latvec_times = Rcut_max * rcut_times / (proj * ucell.lat0);
    const ModuleBase::Vector3<Tcell> latvec_times_ceil = {static_cast<Tcell>(std::ceil(latvec_times.x)),
                                                          static_cast<Tcell>(std::ceil(latvec_times.y)),
                                                          static_cast<Tcell>(std::ceil(latvec_times.z))};
    const ModuleBase::Vector3<Tcell> period = 2 * latvec_times_ceil + ModuleBase::Vector3<Tcell>{1,1,1};
	return std::array<Tcell,3>{period.x, period.y, period.z};
}

template<typename TA, typename Tcell, typename Tdata>
std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>>
LRI_CV_Tools::get_CVws(
	const UnitCell &ucell,
	const std::map<TA,std::map<std::pair<TA,std::array<Tcell,3>>,RI::Tensor<Tdata>>> &CVs)
{
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>> CVws;
	for(const auto &CVs_A : CVs)
	{
		const TA iat0 = CVs_A.first;
		const int it0 = ucell.iat2it[iat0];
		const int ia0 = ucell.iat2ia[iat0];
		const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];
		for(const auto &CVs_B : CVs_A.second)
		{
			const TA iat1 = CVs_B.first.first;
			const int it1 = ucell.iat2it[iat1];
			const int ia1 = ucell.iat2ia[iat1];
			const std::array<int,3> &cell1 = CVs_B.first.second;
			const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
			const Abfs::Vector3_Order<double> R_delta = -tau0+tau1+(RI_Util::array3_to_Vector3(cell1)*ucell.latvec);
			CVws[it0][it1][R_delta] = CVs_B.second;
		}
	}
	return CVws;
}

template <typename TA, typename Tcell, typename Tdata>
std::map<int, std::map<int, std::map<Abfs::Vector3_Order<double>, std::array<RI::Tensor<Tdata>, 3>>>> LRI_CV_Tools::
    get_dCVws(const UnitCell& ucell,
              const std::map<TA, std::map<std::pair<TA, std::array<Tcell, 3>>, std::array<RI::Tensor<Tdata>, 3>>>& dCVs)
{
    std::map<int, std::map<int, std::map<Abfs::Vector3_Order<double>, std::array<RI::Tensor<Tdata>, 3>>>> dCVws;
    for (const auto& dCVs_A: dCVs)
    {
        const TA iat0 = dCVs_A.first;
        const int it0 = ucell.iat2it[iat0];
        const int ia0 = ucell.iat2ia[iat0];
        const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];
        for (const auto& dCVs_B: dCVs_A.second)
        {
            const TA iat1 = dCVs_B.first.first;
            const int it1 = ucell.iat2it[iat1];
            const int ia1 = ucell.iat2ia[iat1];
            const std::array<int, 3>& cell1 = dCVs_B.first.second;
            const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
            const Abfs::Vector3_Order<double> R_delta
                = -tau0 + tau1 + (RI_Util::array3_to_Vector3(cell1) * ucell.latvec);
            dCVws[it0][it1][R_delta] = dCVs_B.second;
        }
    }
    return dCVws;
}

template <typename T, std::size_t N>
void LRI_CV_Tools::init_elem(std::array<RI::Tensor<T>, N>& data,
                             const size_t ndim0,
                             const size_t ndim1) {
    for (size_t i = 0; i < N; ++i) {
        data[i] = RI::Tensor<T>({ndim0, ndim1});
    }
}

template <typename T, std::size_t N>
void LRI_CV_Tools::add_elem(std::array<T, N>& data,
                            const T& val,
                            const T& frac) {
    for (size_t i = 0; i < N; ++i)
        data[i] += frac * val;
}

template <typename T, std::size_t N>
void LRI_CV_Tools::add_elem(std::array<RI::Tensor<T>, N>& data,
                            const int lmp,
                            const int lmq,
                            const std::array<T, N>& val,
                            const T& frac) {
    for (size_t i = 0; i < N; ++i) {
        data[i](lmp, lmq) += frac * val[i];
    }
}

template <typename T, std::size_t N>
void LRI_CV_Tools::add_elem(std::array<RI::Tensor<T>, N>& data,
                            const int lmp0,
                            const int lmq0,
                            const std::array<RI::Tensor<T>, N>& val,
                            const int lmp1,
                            const int lmq1,
                            const T& frac) {
    for (size_t i = 0; i < N; ++i) {
        data[i](lmp0, lmq0) += frac * val[i](lmp1, lmq1);
    }
}

template <typename Tout, typename Tin>
RI::Tensor<Tout> LRI_CV_Tools::convert(RI::Tensor<Tin>&& data) {
    return RI::Global_Func::convert<Tout>(data);
}

template <typename Tout, typename Tin, std::size_t N>
std::array<RI::Tensor<Tout>, N>
    LRI_CV_Tools::convert(std::array<RI::Tensor<Tin>, N>&& data) {
    std::array<RI::Tensor<Tout>, N> out;
    for (size_t i = 0; i != N; ++i)
        out[i] = RI::Global_Func::convert<Tout>(data[i]);
    return out;
}

// template <typename T>
// RI::Tensor<T> LRI_CV_Tools::check_zero(RI::Tensor<T>&& data) {
//     RI::Tensor<T> result(data.shape);

//     const std::size_t rows = data.shape[0];
//     const std::size_t cols = data.shape[1];

//     for (std::size_t i = 0; i < rows; ++i) {
//         for (std::size_t j = 0; j < cols; ++j) {
//             result(i, j) = LRI_CV_Tools::check_zero(data(i, j));
//         }
//     }

//     return result;
// }

// template <typename T, std::size_t N>
// std::array<RI::Tensor<T>, N>
//     LRI_CV_Tools::check_zero(std::array<RI::Tensor<T>, N>&& data) {
//     std::array<RI::Tensor<T>, N> result;

//     for (size_t i = 0; i != N; ++i)
//         result[i] = LRI_CV_Tools::check_zero(std::move(data[i]));

//     return result;
// }


// dMRs[ipos0][ipos1] = \nabla_{ipos0} M R_{ipos1}
template<typename TA, typename TC, typename Tdata>
std::array<std::array<std::map<TA,std::map<std::pair<TA,TC>,RI::Tensor<Tdata>>>,3>,3>
LRI_CV_Tools::cal_dMRs(
	const UnitCell &ucell,
	const std::array<std::map<TA,std::map<std::pair<TA,TC>,RI::Tensor<Tdata>>>,3> &dMs)
{
	auto get_R_delta = [&](const TA &iat0, const std::pair<TA,TC> &A1) -> std::array<Tdata,3>
	{
		const TA iat1 = A1.first;
		const TC &cell1 = A1.second;
		const int it0 = ucell.iat2it[iat0];
		const int ia0 = ucell.iat2ia[iat0];
		const int it1 = ucell.iat2it[iat1];
		const int ia1 = ucell.iat2ia[iat1];
		const ModuleBase::Vector3<double> tau0 = ucell.atoms[it0].tau[ia0];
		const ModuleBase::Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
		const Abfs::Vector3_Order<double> R_delta = -tau0+tau1+(RI_Util::array3_to_Vector3(cell1)*ucell.latvec);
		return std::array<Tdata,3>{R_delta.x, R_delta.y, R_delta.z};
	};
	constexpr int Npos = 3;
	std::array<std::array<std::map<TA,std::map<std::pair<TA,TC>,RI::Tensor<Tdata>>>,Npos>,Npos> dMRs;
	for(int ipos0=0; ipos0<Npos; ++ipos0)
	{
		for(int ipos1=0; ipos1<Npos; ++ipos1)
		{
			for(const auto &dMs_A : dMs[ipos0])
			{
				const TA iat0 = dMs_A.first;
				for(const auto &dMs_B : dMs_A.second)
				{
					const std::pair<TA,TC> A1 = dMs_B.first;
					const RI::Tensor<Tdata> &dM = dMs_B.second;
					const std::array<Tdata,3> R_delta = get_R_delta(iat0, A1);
					dMRs[ipos0][ipos1][iat0][A1] = dM * R_delta[ipos1];
				}
			}
		}
	}
	return dMRs;
}

#endif

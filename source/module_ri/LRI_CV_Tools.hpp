//=======================
// AUTHOR : Peize Lin
// DATE :   2022-10-24
//=======================

#ifndef LRI_CV_TOOLS_HPP
#define LRI_CV_TOOLS_HPP

#include "LRI_CV_Tools.h"
#include "Inverse_Matrix.h"
#include "../module_base/mathzone.h"

template<typename Tdata>
RI::Tensor<Tdata>
LRI_CV_Tools::cal_I( const RI::Tensor<Tdata> &m )
{
	Inverse_Matrix<Tdata> I;
	I.input(m);
	I.cal_inverse( Inverse_Matrix<Tdata>::Method::potrf );
	return I.output();
}

template<typename Tdata>
std::vector<std::vector<RI::Tensor<Tdata>>>
LRI_CV_Tools::cal_I( const std::vector<std::vector<RI::Tensor<Tdata>>> &ms )
{
	Inverse_Matrix<Tdata> I;
	I.input(ms);
	I.cal_inverse( Inverse_Matrix<Tdata>::Method::potrf );
	return I.output({ms[0][0].shape[0], ms[1][0].shape[0]}, {ms[0][0].shape[1], ms[0][1].shape[1]});
}



template<typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::transform_Rm(const RI::Tensor<Tdata> &V)
{
	return V.transpose();
}

template<typename Tdata>
std::array<RI::Tensor<Tdata>,3> LRI_CV_Tools::transform_Rm(const std::array<RI::Tensor<Tdata>,3> &dV)
{
	return std::array<RI::Tensor<Tdata>,3>{-dV[0].transpose(), -dV[1].transpose(), -dV[2].transpose()};
}

template<typename Tdata>
bool LRI_CV_Tools::exist(const RI::Tensor<Tdata> &V)
{
	return !V.empty();
}

template<typename T, std::size_t N>
bool LRI_CV_Tools::exist(const std::array<T,N> &dV)
{
	for(size_t i=0; i<3; ++i)
		if(!dV[i].empty())
			return true;
	return false;
}


template<typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::mul1(
	const RI::Tensor<Tdata> &t1,
	const RI::Tensor<Tdata> &t2)
{
	const size_t sa0=t1.shape[0], sa1=t2.shape[0], sl0=t2.shape[1], sl1=t2.shape[2];
	return (t1 * t2.reshape({sa1,sl0*sl1})).reshape({sa0,sl0,sl1});
}
template<typename T>
std::array<T,3> LRI_CV_Tools::mul1(
	const std::array<T,3> &t1,
	const T &t2)
{
	return std::array<T,3>{
		mul1(t1[0],t2), mul1(t1[1],t2), mul1(t1[2],t2) };
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

template<typename Tdata>
std::vector<RI::Tensor<Tdata>> LRI_CV_Tools::mul2(
	const std::vector<std::vector<RI::Tensor<Tdata>>> &mat,
	const std::vector<RI::Tensor<Tdata>> &vec)
{
	const size_t sa0=vec[0].shape[0], sa1=vec[1].shape[0], sl0=vec[0].shape[1], sl1=vec[0].shape[2];
	const RI::Tensor<Tdata> vec0=vec[0].reshape({sa0,sl0*sl1}), vec1=vec[1].reshape({sa1,sl0*sl1});
	return std::vector<RI::Tensor<Tdata>>
		{( mat[0][0]*vec0 + mat[0][1]*vec1 ).reshape({sa0,sl0,sl1}),
		 ( mat[1][0]*vec0 + mat[1][1]*vec1 ).reshape({sa1,sl0,sl1})};
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
template<typename T1, typename T2>
std::array<T2,3> LRI_CV_Tools::mul2(
	const T1 &t1,
	const std::array<T2,3> &t2)
{
	return std::array<T2,3>{
		mul2(t1,t2[0]), mul2(t1,t2[1]), mul2(t1,t2[2]) };
}

/*
template<typename T, std::size_t N>
std::array<T,N> LRI_CV_Tools::operator-(const std::array<T,N> &v1, const std::array<T,N> &v2)
{
	std::array<T,N> v;
	for(std::size_t i=0; i<N; ++i)
		v[i] = v1[i] - v2[i];
	return v;
}
template<typename T>
std::vector<T> LRI_CV_Tools::operator-(const std::vector<T> &v1, const std::vector<T> &v2)
{
	assert(v1.size()==v2.size());
	std::vector<T> v(v1.size());
	for(std::size_t i=0; i<v.size(); ++i)
		v[i] = v1[i] - v2[i];
	return v;
}
*/
template<typename T, std::size_t N>
std::vector<std::array<T,N>> LRI_CV_Tools::minus(
	const std::vector<std::array<T,N>> &v1,
	const std::vector<std::array<T,N>> &v2)
{
	assert(v1.size()==v2.size());
	std::vector<std::array<T,N>> v(v1.size());
	for(std::size_t i=0; i<v.size(); ++i)
		for(std::size_t j=0; j<N; ++j)
			v[i][j] = v1[i][j] - v2[i][j];
	return v;
}


template<typename T, std::size_t N>
std::array<T,N> LRI_CV_Tools::negative(const std::array<T,N> &v_in)
{
	std::array<T,N> v_out;
	for(std::size_t i=0; i<N; ++i)
		v_out[i] = -v_in[i];
	return v_out;
}


template<typename Tdata>
RI::Tensor<Tdata> LRI_CV_Tools::transpose12(const RI::Tensor<Tdata> &c_in)
{
	RI::Tensor<Tdata> c_out({c_in.shape[0], c_in.shape[2], c_in.shape[1]});
	for(size_t i0=0; i0<c_in.shape[0]; ++i0)
		for(size_t i1=0; i1<c_in.shape[1]; ++i1)
			for(size_t i2=0; i2<c_in.shape[2]; ++i2)
				c_out(i0,i2,i1) = c_in(i0,i1,i2);
	return c_out;
}

template<typename T, std::size_t N>
std::array<T,N> LRI_CV_Tools::transpose12(const std::array<T,N> &c_in)
{
	std::array<T,N> c_out;
	for(size_t i=0; i<N; ++i)
		c_out[i] = transpose12(c_in[i]);
	return c_out;
}


template<typename T, std::size_t N>
std::array<std::vector<T>,N>
LRI_CV_Tools::change_order(std::vector<std::array<T,N>> &&ds_in)
{
	std::array<std::vector<T>,N> ds;
	for(int ix=0; ix<N; ++ix)
	{
		ds[ix].resize(ds_in.size());
		for(int iv=0; iv<ds_in.size(); ++iv)
			ds[ix][iv] = std::move(ds_in[iv][ix]);
	}
	return ds;
}

template<typename T, std::size_t N>
std::vector<std::array<T,N>>
LRI_CV_Tools::change_order(std::array<std::vector<T>,N> &&ds_in)
{
	std::vector<std::array<T,N>> ds(ds_in[0].size());
	for(int ix=0; ix<N; ++ix)
	{
		assert(ds.size() == ds_in[ix].size());
		for(int iv=0; iv<ds.size(); ++iv)
			ds[iv][ix] = std::move(ds_in[ix][iv]);
	}
	return ds;
}

template<typename T, std::size_t N>
std::array<std::vector<std::vector<T>>,N>
LRI_CV_Tools::change_order(std::vector<std::vector<std::array<T,N>>> &&ds_in)
{
	std::array<std::vector<std::vector<T>>,N> ds;
	for(int ix=0; ix<N; ++ix)
	{
		ds[ix].resize(ds_in.size());
		for(int i0=0; i0<ds_in.size(); ++i0)
		{
			ds[ix][i0].resize(ds_in[i0].size());
			for(int i1=0; i1<ds_in[i0].size(); ++i1)
				ds[ix][i0][i1] = std::move(ds_in[i0][i1][ix]);
		}
	}
	return ds;
}

template<typename TkeyA, typename TkeyB, typename Tvalue, std::size_t N>
std::array<std::map<TkeyA,std::map<TkeyB,Tvalue>>,N>
LRI_CV_Tools::change_order(std::map<TkeyA,std::map<TkeyB,std::array<Tvalue,N>>> && ds_in)
{
	std::array<std::map<TkeyA,std::map<TkeyB,Tvalue>>,N> ds;
	for(auto &ds_A : ds_in)
		for(auto &ds_B : ds_A.second)
			for(int ix=0; ix<N; ++ix)
				ds[ix][ds_A.first][ds_B.first] = std::move(ds_B.second[ix]);
	return ds;
}


template<typename Tcell>
std::array<Tcell,3>
LRI_CV_Tools::cal_latvec_range(const double &rcut_times)
{
	double Rcut_max = 0;
	for(int T=0; T<GlobalC::ucell.ntype; ++T)
		Rcut_max = std::max(Rcut_max, GlobalC::ORB.Phi[T].getRcut());
	const ModuleBase::Vector3<double> proj = ModuleBase::Mathzone::latvec_projection(
		std::array<ModuleBase::Vector3<double>,3>{GlobalC::ucell.a1, GlobalC::ucell.a2, GlobalC::ucell.a3});
	const ModuleBase::Vector3<double> latvec_times = Rcut_max * rcut_times / (proj * GlobalC::ucell.lat0);
    const ModuleBase::Vector3<Tcell> latvec_times_ceil = {static_cast<Tcell>(std::ceil(latvec_times.x)),
                                                          static_cast<Tcell>(std::ceil(latvec_times.y)),
                                                          static_cast<Tcell>(std::ceil(latvec_times.z))};
    const ModuleBase::Vector3<Tcell> period = 2 * latvec_times_ceil + ModuleBase::Vector3<Tcell>{1,1,1};
	return std::array<Tcell,3>{period.x, period.y, period.z};
}

template<typename TA, typename Tcell, typename Tdata>
std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>>
LRI_CV_Tools::get_CVws(
	const std::map<TA,std::map<std::pair<TA,std::array<Tcell,3>>,RI::Tensor<Tdata>>> &CVs)
{
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,RI::Tensor<Tdata>>>> CVws;
	for(const auto &CVs_A : CVs)
	{
		const TA iat0 = CVs_A.first;
		const int it0 = GlobalC::ucell.iat2it[iat0];
		const int ia0 = GlobalC::ucell.iat2ia[iat0];
		const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it0].tau[ia0];
		for(const auto &CVs_B : CVs_A.second)
		{
			const TA iat1 = CVs_B.first.first;
			const int it1 = GlobalC::ucell.iat2it[iat1];
			const int ia1 = GlobalC::ucell.iat2ia[iat1];
			const std::array<int,3> &cell1 = CVs_B.first.second;
			const ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
			const Abfs::Vector3_Order<double> R_delta = -tau0+tau1+(RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec);
			CVws[it0][it1][R_delta] = CVs_B.second;
		}
	}
	return CVws;
}

template<typename TA, typename Tcell, typename Tdata>
std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,std::array<RI::Tensor<Tdata>,3>>>>
LRI_CV_Tools::get_dCVws(
	const std::array<std::map<TA,std::map<std::pair<TA,std::array<Tcell,3>>,RI::Tensor<Tdata>>>,3> &dCVs)
{
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,std::array<RI::Tensor<Tdata>,3>>>> dCVws;
	for(int ix=0; ix<3; ++ix)
	{
		for(const auto &dCVs_A : dCVs[ix])
		{
			const TA iat0 = dCVs_A.first;
			const int it0 = GlobalC::ucell.iat2it[iat0];
			const int ia0 = GlobalC::ucell.iat2ia[iat0];
			const ModuleBase::Vector3<double> tau0 = GlobalC::ucell.atoms[it0].tau[ia0];
			for(const auto &dCVs_B : dCVs_A.second)
			{
				const TA iat1 = dCVs_B.first.first;
				const int it1 = GlobalC::ucell.iat2it[iat1];
				const int ia1 = GlobalC::ucell.iat2ia[iat1];
				const std::array<int,3> &cell1 = dCVs_B.first.second;
				const ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[it1].tau[ia1];
				const Abfs::Vector3_Order<double> R_delta = -tau0+tau1+(RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec);
				dCVws[it0][it1][R_delta][ix] = dCVs_B.second;
			}
		}
	}
	return dCVws;
}

#endif
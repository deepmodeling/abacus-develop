//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef RI_UTIL_H
#define RI_UTIL_H

#include <RI/global/Array_Operator.h>
#include <RI/global/Global_Func-2.h>

#include <array>

namespace RI_Util
{
	inline extern std::array<int,3>
	get_Born_vonKarmen_period();

	template<typename Tcell, size_t Ndim>
	extern std::vector<std::array<Tcell,Ndim>>
	get_Born_von_Karmen_cells( const std::array<Tcell,Ndim> &Born_von_Karman_period );

	template<typename Tcell>
	inline std::array<Tcell,3>
	Vector3_to_array3(const ModuleBase::Vector3<Tcell> &v)
	{
		return std::array<Tcell,3> {v.x, v.y, v.z};
	}
	template<typename Tcell>
	inline ModuleBase::Vector3<Tcell>
	array3_to_Vector3(const std::array<Tcell,3> &v)
	{
		return ModuleBase::Vector3<Tcell> {v[0], v[1], v[2]};
	}

	template<typename Tdata, typename Tmatrix>
	RI::Tensor<Tdata>
	Matrix_to_Tensor(const Tmatrix &m_old)
	{
		RI::Tensor<Tdata> m_new({static_cast<size_t>(m_old.nr), static_cast<size_t>(m_old.nc)});
		for(int ir=0; ir<m_old.nr; ++ir)
			for(int ic=0; ic<m_old.nc; ++ic)
				m_new(ir,ic) = RI::Global_Func::convert<Tdata>(m_old(ir,ic));
		return m_new;
    }

    template<typename Tdata>
    RI::Tensor<Tdata>
        Vector_to_Tensor(const std::vector<Tdata>& m_old, const int nr, const int nc)
    {
        assert(nr * nc == m_old.size());
        RI::Tensor<Tdata> m_new({ static_cast<size_t>(nr), static_cast<size_t>(nc) });
        for (int ir = 0; ir < nr; ++ir)
            for (int ic = 0; ic < nc; ++ic)
                m_new(ir, ic) = RI::Global_Func::convert<Tdata>(m_old[ir * nc + ic]);
        return m_new;
    }
}

#include "RI_Util.hpp"

#endif
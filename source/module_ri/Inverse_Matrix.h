//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#pragma once

#include <RI/global/Tensor.h>
#include <vector>

template<typename Tdata>
class Inverse_Matrix
{
public:
	enum class Method{potrf};	//, syev};
	void cal_inverse(const Method &method);

	void input(const RI::Tensor<Tdata> &m);
	void input(const std::vector<std::vector<RI::Tensor<Tdata>>> &ms);
	RI::Tensor<Tdata> output() const;
	std::vector<std::vector<RI::Tensor<Tdata>>> output(const std::vector<size_t> &n0, const std::vector<size_t> &n1) const;

private:
	void using_potrf();
	void copy_down_triangle();
	RI::Tensor<Tdata> A;
};

#include "Inverse_Matrix.hpp"
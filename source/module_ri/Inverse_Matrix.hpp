//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef INVERSE_MATRIX_HPP
#define INVERSE_MATRIX_HPP

#include "Inverse_Matrix.h"
#include "module_base/lapack_connector.h"

#include <cassert>

template<typename Tdata>
void Inverse_Matrix<Tdata>::cal_inverse( const Method &method )
{
	switch(method)
	{
		case Method::potrf:		using_potrf();			break;
//		case Method::syev:		using_syev(1E-6);		break;
	}
}

template<typename Tdata>
void Inverse_Matrix<Tdata>::using_potrf()
{
	int info;
	LapackConnector::potrf('U', A.shape[0], A.ptr(), A.shape[0], info);
	if(info)
		throw std::range_error("info="+std::to_string(info)+"\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

	LapackConnector::potri('U', A.shape[0], A.ptr(), A.shape[0], info);
	if(info)
		throw std::range_error("info="+std::to_string(info)+"\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

	copy_down_triangle();
}

/*
void Inverse_Matrix::using_syev( const double &threshold_condition_number )
{
	std::vector<double> eigen_value(A.nr);
	LapackConnector::dsyev('V','U',A,eigen_value.data(),info);

	double eigen_value_max = 0;
	for( const double &ie : eigen_value )
		eigen_value_max = std::max( ie, eigen_value_max );
	const double threshold = eigen_value_max * threshold_condition_number;

	ModuleBase::matrix eA( A.nr, A.nc );
	int ie=0;
	for( int i=0; i!=A.nr; ++i )
		if( eigen_value[i] > threshold )
		{
			BlasConnector::axpy( A.nc, sqrt(1.0/eigen_value[i]), A.c+i*A.nc,1, eA.c+ie*eA.nc,1 );
			++ie;
		}
	BlasConnector::gemm( 'T','N', eA.nc,eA.nc,ie, 1, eA.c,eA.nc, eA.c,eA.nc, 0, A.c,A.nc );
}
*/

template<typename Tdata>
void Inverse_Matrix<Tdata>::input( const RI::Tensor<Tdata> &m )
{
	assert(m.shape.size()==2);
	assert(m.shape[0]==m.shape[1]);
	this->A = m.copy();
}


template<typename Tdata>
void Inverse_Matrix<Tdata>::input(const std::vector<std::vector<RI::Tensor<Tdata>>> &ms)
{
	const size_t N0 = ms.size();
	assert(N0>0);
	const size_t N1 = ms[0].size();
	assert(N1>0);
	for(size_t Im0=0; Im0<N0; ++Im0)
		assert(ms[Im0].size()==N1);

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
			assert(ms[Im0][Im1].shape.size()==2);

	std::vector<size_t> n0(N0);
	for(size_t Im0=0; Im0<N0; ++Im0)
		n0[Im0] = ms[Im0][0].shape[0];
	std::vector<size_t> n1(N1);
	for(size_t Im1=0; Im1<N1; ++Im1)
		n1[Im1] = ms[0][Im1].shape[1];

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
			assert((ms[Im0][Im1].shape[0]==n0[Im0]) && (ms[Im0][Im1].shape[1]==n1[Im1]));

	const size_t n_all = std::accumulate(n0.begin(), n0.end(), 0);
	assert(n_all == std::accumulate(n1.begin(), n1.end(), 0));
	this->A = RI::Tensor<Tdata>({n_all, n_all});

	std::vector<size_t> n0_partial(N0+1);
	std::partial_sum(n0.begin(), n0.end(), n0_partial.begin()+1);
	std::vector<size_t> n1_partial(N1+1);
	std::partial_sum(n1.begin(), n1.end(), n1_partial.begin()+1);

	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
		{
			const RI::Tensor<Tdata> &m_tmp = ms.at(Im0).at(Im1);
			for(size_t im0=0; im0<m_tmp.shape[0]; ++im0)
				for(size_t im1=0; im1<m_tmp.shape[1]; ++im1)
					this->A(im0+n0_partial[Im0], im1+n1_partial[Im1]) = m_tmp(im0,im1);
		}
}


template<typename Tdata>
RI::Tensor<Tdata> Inverse_Matrix<Tdata>::output() const
{
	return this->A.copy();
}


template<typename Tdata>
std::vector<std::vector<RI::Tensor<Tdata>>>
Inverse_Matrix<Tdata>::output(const std::vector<size_t> &n0, const std::vector<size_t> &n1) const
{
	assert( std::accumulate(n0.begin(), n0.end(), 0) == this->A.shape[0] );
	assert( std::accumulate(n1.begin(), n1.end(), 0) == this->A.shape[1] );

	const size_t N0 = n0.size();
	const size_t N1 = n1.size();

	std::vector<size_t> n0_partial(N0+1);
	std::partial_sum(n0.begin(), n0.end(), n0_partial.begin()+1);
	std::vector<size_t> n1_partial(N1+1);
	std::partial_sum(n1.begin(), n1.end(), n1_partial.begin()+1);

	std::vector<std::vector<RI::Tensor<Tdata>>> ms(N0, std::vector<RI::Tensor<Tdata>>(N1));
	for(size_t Im0=0; Im0<N0; ++Im0)
		for(size_t Im1=0; Im1<N1; ++Im1)
		{
			RI::Tensor<Tdata> &m_tmp = ms[Im0][Im1] = RI::Tensor<Tdata>({n0[Im0], n1[Im1]});
			for(size_t im0=0; im0<n0[Im0]; ++im0)
				for(size_t im1=0; im1<n1[Im1]; ++im1)
					m_tmp(im0,im1) = this->A(im0+n0_partial[Im0], im1+n1_partial[Im1]);
		}
	return ms;
}


template<typename Tdata>
void Inverse_Matrix<Tdata>::copy_down_triangle()
{
	for( size_t i0=0; i0<A.shape[0]; ++i0 )
		for( size_t i1=0; i1<i0; ++i1 )
			A(i0,i1) = A(i1,i0);
}

#endif
#include "exx_abfs-inverse_matrix_double.h"
#include "module_base/lapack_connector.h"
#include <cstring>

#include "src_external/src_test/src_global/matrix-test.h"		// Peize Lin test

void Exx_Abfs::Inverse_Matrix_Double::init(const int &dim_in)
{
	this->dim = dim_in;
	assert(dim>0);

	this->info = 0;
	this->A.create(dim, dim,false);
}

void Exx_Abfs::Inverse_Matrix_Double::cal_inverse( const Method &method )
{
	#if TEST_EXX_LCAO==1
		static int i=0;
		ofstream ofs("inverse_matrix_"+TO_STRING(i));
		ofs<<A<<endl;
		ofs.close();
		++i;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	switch(method)
	{
		case Method::dpotrf:	using_dpotrf();			break;
		case Method::dsyev:		using_dsyev(1E-6);		break;
	}
}

void Exx_Abfs::Inverse_Matrix_Double::using_dpotrf()
{
	LapackConnector::dpotrf('U',dim,A,dim,&info);

	if(info!=0)
	{
		cout << "\n info_dpotrf = " << info <<endl;
		ofs_warning<<A<<endl;
		QUIT();
	}

	LapackConnector::dpotri('U',dim,A,dim,&info);
	
	if(info!=0)
	{
		cout << "\n info_dpotri = " << info <<endl;
		ofs_warning<<A<<endl;
		QUIT();
	}
	
	copy_down_triangle();
}

void Exx_Abfs::Inverse_Matrix_Double::using_dsyev( const double &threshold_condition_number )
{
	vector<double> eigen_value(A.nr);
	LapackConnector::dsyev('V','U',A,VECTOR_TO_PTR(eigen_value),info);
	
	#if TEST_EXX_LCAO==1
		for( const double &ie : eigen_value )
			cout<<ie<<"\t";
		cout<<endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	double eigen_value_max = 0;
	for( const double &ie : eigen_value )
		eigen_value_max = max( ie, eigen_value_max );
	const double threshold = eigen_value_max * threshold_condition_number;

	#if TEST_EXX_LCAO==1
		cout<<eigen_value_max<<"\t"<<threshold_condition_number<<"\t"<<threshold<<endl;
	#elif TEST_EXX_LCAO==-1
		#error "TEST_EXX_LCAO"
	#endif
	
	/*
	matrix eigen_value_inverse(A.nr,A.nc);
	for( size_t i=0; i!=A.nr; ++i )
		eigen_value_inverse(i,i) = (eigen_value[i]>threshold) ? (1.0/eigen_value[i]) : 0;
	A = transpose(A) * eigen_value_inverse * A;
	*/
	
	matrix eA( A.nr, A.nc );
	int ie=0;
	for( int i=0; i!=A.nr; ++i )
		if( eigen_value[i] > threshold )
		{
			LapackConnector::axpy( A.nc, sqrt(1.0/eigen_value[i]), A.c+i*A.nc,1, eA.c+ie*eA.nc,1 );
			++ie;
		}
	LapackConnector::gemm( 'T','N', eA.nc,eA.nc,ie, 1, eA.c,eA.nc, eA.c,eA.nc, 0, A.c,A.nc );
}


void Exx_Abfs::Inverse_Matrix_Double::input( const matrix &m )
{
	for( size_t ir=0; ir!=m.nr; ++ir )
		for( size_t ic=0; ic!=m.nc; ++ic )
			A(ir,ic) = m(ir,ic);
}

void Exx_Abfs::Inverse_Matrix_Double::input( 
	const matrix &m_00,
	const matrix &m_01,
	const matrix &m_10,
	const matrix &m_11)
{
	const size_t delta_nr = m_00.nr, delta_nc = m_00.nc;
	
	for( size_t ir=0; ir!=m_00.nr; ++ir )
		for( size_t ic=0; ic!=m_00.nc; ++ic )
			A(ir,ic) = m_00(ir,ic);	
	for( size_t ir=0; ir!=m_01.nr; ++ir )
		for( size_t ic=0; ic!=m_01.nc; ++ic )
			A(ir,ic+delta_nc) = m_01(ir,ic);	
	for( size_t ir=0; ir!=m_10.nr; ++ir )
		for( size_t ic=0; ic!=m_10.nc; ++ic )
			A(ir+delta_nr,ic) = m_10(ir,ic);	
	for( size_t ir=0; ir!=m_11.nr; ++ir )
		for( size_t ic=0; ic!=m_11.nc; ++ic )
			A(ir+delta_nr,ic+delta_nc) = m_11(ir,ic);
}

void Exx_Abfs::Inverse_Matrix_Double::output( matrix &m ) const
{
	for( size_t ir=0; ir!=m.nr; ++ir )
		for( size_t ic=0; ic!=m.nc; ++ic )
			m(ir,ic) = A(ir,ic);
}

void Exx_Abfs::Inverse_Matrix_Double::output( 
	matrix &m_00,
	matrix &m_01,
	matrix &m_10,
	matrix &m_11) const
{
	const size_t delta_nr = m_00.nr, delta_nc = m_00.nc;
	
	for( size_t ir=0; ir!=m_00.nr; ++ir )
		for( size_t ic=0; ic!=m_00.nc; ++ic )
			m_00(ir,ic) = A(ir,ic);	
	for( size_t ir=0; ir!=m_01.nr; ++ir )
		for( size_t ic=0; ic!=m_01.nc; ++ic )
			m_01(ir,ic) = A(ir,ic+delta_nc);	
	for( size_t ir=0; ir!=m_10.nr; ++ir )
		for( size_t ic=0; ic!=m_10.nc; ++ic )
			m_10(ir,ic) = A(ir+delta_nr,ic);	
	for( size_t ir=0; ir!=m_11.nr; ++ir )
		for( size_t ic=0; ic!=m_11.nc; ++ic )
			m_11(ir,ic) = A(ir+delta_nr,ic+delta_nc);
}

void Exx_Abfs::Inverse_Matrix_Double::copy_down_triangle()
{
	for( size_t ir=0; ir!=A.nr; ++ir )
		for( size_t ic=0; ic!=ir; ++ic )
			A(ir,ic) = A(ic,ir);
}
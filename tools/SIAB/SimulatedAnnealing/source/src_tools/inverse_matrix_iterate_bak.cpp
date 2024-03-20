#include "inverse_matrix_iterate.h"

void Inverse_Matrix_Iterate::init(size_t dim_in)
{
	this->matrix_info.dim = dim_in;
	this->X_store_1.create(dim_in,dim_in);
	this->X_store_2.create(dim_in,dim_in);
	
	this->X0_tmp.create(dim_in,dim_in);
	this->identity.create(dim_in,dim_in);
	this->identity.set_as_identity_matrix();
	this->error_matrix.create(dim_in,dim_in);
}

Inverse_Matrix_Iterate &Inverse_Matrix_Iterate::set_X0(const ComplexMatrix &X0_in)
{
	if( &X0_in == &this->X_store_1 )		// without this will cause a bug when istep==1
	{
		this->X0_ptr = &this->X0_tmp;
		this->X0_tmp = X0_in;
	}
	else
	{
		this->X0_ptr = &X0_in;
	}
	return *this;
}

void Inverse_Matrix_Iterate::iterate(const ComplexMatrix &old_X, ComplexMatrix &new_X)
{
	// I - A * X_n
	if(this->matrix_info.symmetry) 
	{ this->error_matrix = this->identity - multiply_special(this->matrix_info.symmetry, 'L', this->matrix_info.uplo, *(this->A_ptr), old_X); }
	else 
	{ this->error_matrix = this->identity - *(this->A_ptr) * old_X; }
	
	const double error_norm(sqrt(abs2(this->error_matrix)));
	if ( error_norm < this->tolerate_error)
	{
		this->inverse_ptr = &old_X;
		throw Little_Error();
	}
	else if ( error_norm > this->last_error + this->fluctrate_error)
	{	throw Unstable();	}
	else { this->last_error = error_norm; }
	
	// X_{n+1} = X_n (I + (I-A*X_n) ) = X_n ( 2I - A*X_n )
	for( size_t i(0); i<this->matrix_info.dim; ++i) { this->error_matrix(i,i) += 1.0; }
	new_X = old_X * this->error_matrix;
//	new_X = old_X * (this->identity + this->error_matrix);

	this->inverse_ptr = &new_X;
}

const ComplexMatrix &Inverse_Matrix_Iterate::cal_inverse()
{
	try
	{
		this->inverse_ptr = this->X0_ptr;
		this->last_error = std::numeric_limits<double>::max();
		if(this->iterate_num >= 1)
		{
			iterate( *this->X0_ptr, this->X_store_1 );
		}
		for( size_t istep(2); istep<=this->iterate_num; ++istep)
		{
			if( istep&1 )		// istep odd, X_store_2 -> X_store_1
			{
				iterate( this->X_store_2, this->X_store_1 );
			}
			else				// istep even, X_store_1 -> X_store_2
			{
				iterate( this->X_store_1, this->X_store_2);
			}	
		}		
	}
	catch(const Little_Error &e){}

	return get_inverse();
}
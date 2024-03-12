#ifndef GRAM_SCHMIDT_ORTH_H
#define GRAM_SCHMIDT_ORTH_H

#include<cstddef>
#include<exception>
using std::exception;

class GS_Orth_Coordinate: public exception{};

class Gram_Schmidt_Orth
{
public:

	~Gram_Schmidt_Orth();
	double **orth();
	
private:
	int coordinate;			// 0:cartesian, 1:sphere
	size_t func_num;
	size_t variable_num;
	double **func{nullptr};
	double *r{nullptr};
	double *r_2{nullptr};	
	double *rab{nullptr};

private:
	bool finish_init{false};
	
public:

	size_t get_coordinate() const { return coordinate; }
	Gram_Schmidt_Orth &set_coordinate(size_t coordinate_in)
	{
		if( 0 != coordinate_in && 1 != coordinate_in) { throw GS_Orth_Coordinate();}
		coordinate = coordinate_in; 
		return *this; 
	}

	size_t get_func_num() const { return func_num; }
	Gram_Schmidt_Orth &set_func_num(size_t func_num_in){ func_num = func_num_in; return *this; }

	size_t get_variable_num() const { return variable_num; }
	Gram_Schmidt_Orth &set_variable_num(size_t variable_num_in){ variable_num = variable_num_in; return *this; }

	double **get_func() const { return func; }
	Gram_Schmidt_Orth &set_func(double **func_in){ func = func_in; return *this; }

	double *get_r() const { return r; }
	Gram_Schmidt_Orth &set_r(double *r_in){ r = r_in; return *this; }

	double *get_r_2() const { return r_2; }

	double *get_rab() const { return rab; }
	Gram_Schmidt_Orth &set_rab_element(const double &rab_element_in)
	{
		if(NULL==rab)
			rab = new double [variable_num];
		for( size_t i(0); i<variable_num; ++i) 
			rab[i] = rab_element_in;
		return *this; 
	}
		
private:
	void init();
	double cal_norm(double *func_in);

};

#endif
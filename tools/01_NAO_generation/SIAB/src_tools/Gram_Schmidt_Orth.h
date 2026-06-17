#ifndef GRAM_SCHMIDT_ORTH_H
#define GRAM_SCHMIDT_ORTH_H

#include<cstddef>
#include<exception>
using std::exception;

class GS_Orth_Coordinate: public exception{};

template<typename Func_Type, typename R_Type=double>
class Gram_Schmidt_Orth
{
public:

	~Gram_Schmidt_Orth();
	Func_Type **orth();
	
private:
	int coordinate;			// 0:cartesian, 1:sphere
	size_t func_num;
	size_t variable_num;
	Func_Type **func{nullptr};
	R_Type *r{nullptr};
	R_Type *r_2{nullptr};	
	R_Type *rab{nullptr};

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

	Func_Type **get_func() const { return func; }
	Gram_Schmidt_Orth &set_func(Func_Type **func_in){ func = func_in; return *this; }

	R_Type *get_r() const { return r; }
	Gram_Schmidt_Orth &set_r(R_Type *r_in){ r = r_in; return *this; }

	R_Type *get_r_2() const { return r_2; }

	R_Type *get_rab() const { return rab; }
	Gram_Schmidt_Orth &set_rab_element(const R_Type &rab_element_in)
	{
		if(NULL==rab)
			rab = new R_Type [variable_num];
		for( size_t i(0); i<variable_num; ++i) 
			rab[i] = rab_element_in;
		return *this; 
	}
		
private:
	void init();
	Func_Type cal_norm(Func_Type *func_in);

};

#include "mathzone.h"
#include "lapack_connector.h"

template<typename Func_Type, typename R_Type>
void Gram_Schmidt_Orth<Func_Type,R_Type>::init()
{
	if(this->finish_init) return;
	if(nullptr == this->r_2)
	{
		this->r_2 = new R_Type[variable_num];
		Mathzone::Pointwise_Product( this->variable_num, this->r, this->r, this->r_2);		
	}
	if(nullptr == rab)
	{
		this->rab = new R_Type[variable_num];
		this->rab[0] = 0.0;
		for( size_t ir(1); ir<variable_num; ++ir)
		{
			this->rab[ir] = this->r[ir] - this->r[ir-1];
		}
	}
	this->finish_init = true;
	return;
}

template<typename Func_Type, typename R_Type>
Gram_Schmidt_Orth<Func_Type,R_Type>::~Gram_Schmidt_Orth()
{
	delete[]rab;
	delete[]r_2;
}

template<typename Func_Type, typename R_Type>
Func_Type **Gram_Schmidt_Orth<Func_Type,R_Type>::orth()
{
	init();

	//	Schmidt: hn to en
	//		e1 = h1 / ||h1||
	//		gn = hn - \sum{i=1 to n-1}(hn,ei)ei
	//		en = gn / ||gn||

	for( size_t if1(0); if1<func_num; ++if1)
	{
		for( size_t if2(0); if2<if1; ++if2)
		{
			// (hn,ei)
			Func_Type *mul_func = new Func_Type[this->variable_num];
			Mathzone::Pointwise_Product( this->variable_num, this->func[if1], this->func[if2], mul_func);		
			const Func_Type in_product = cal_norm(mul_func);
			delete[]mul_func; mul_func=nullptr;

			// hn -(hn,ei)ei
			LapackConnector::daxpy( this->variable_num, -in_product, this->func[if2], 1, this->func[if1], 1);
		}
		
		// ||gn||
		Func_Type *func_2 = new Func_Type[this->variable_num];
		Mathzone::Pointwise_Product( this->variable_num, this->func[if1], this->func[if1], func_2);
		const Func_Type norm = sqrt(cal_norm(func_2));
		delete[]func_2; func_2=nullptr;
		
		// en = gn / ||gn||
		LapackConnector::dscal( this->variable_num, 1.0/norm, this->func[if1], 1);
	}
	return func;
}

template<typename Func_Type, typename R_Type>
Func_Type Gram_Schmidt_Orth<Func_Type,R_Type>::cal_norm(Func_Type *func_in)
{
	Func_Type norm(0.0);
	switch(coordinate)
	{
		case 0:
			Mathzone::Simpson_Integral( this->variable_num, func_in, this->rab, norm);		
			break;
		case 1:	
			Func_Type *tmp_func = new Func_Type[this->variable_num];
			Mathzone::Pointwise_Product( this->variable_num, func_in, this->r_2, tmp_func);
			Mathzone::Simpson_Integral( this->variable_num, tmp_func, this->rab, norm);	
			delete[]tmp_func; tmp_func=nullptr;
			break;
	}
	return norm;
}
#endif

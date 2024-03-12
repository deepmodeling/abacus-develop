#include "Gram_Schmidt_Orth.h"
#include "mathzone.h"
#include "lapack_connector.h"

void Gram_Schmidt_Orth::init()
{
	if(this->finish_init) return;
	if(nullptr == this->r_2)
	{
		this->r_2 = new double[variable_num];
		Mathzone::Pointwise_Product( this->variable_num, this->r, this->r, this->r_2);		
	}
	if(nullptr == rab)
	{
		this->rab = new double[variable_num];
		this->rab[0] = 0.0;
		for( size_t ir(1); ir<variable_num; ++ir)
		{
			this->rab[ir] = this->r[ir] - this->r[ir-1];
		}
	}
	this->finish_init = true;
	return;
}

Gram_Schmidt_Orth::~Gram_Schmidt_Orth()
{
	delete[]rab;
	delete[]r_2;
}

double **Gram_Schmidt_Orth::orth()
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
			double *mul_func = new double[this->variable_num];
			Mathzone::Pointwise_Product( this->variable_num, this->func[if1], this->func[if2], mul_func);		
			const double in_product = cal_norm(mul_func);
			delete[]mul_func; mul_func=nullptr;

			// hn -(hn,ei)ei
			daxpy( this->variable_num, -in_product, this->func[if2], 1, this->func[if1], 1);
		}
		
		// ||gn||
		double *func_2 = new double[this->variable_num];
		Mathzone::Pointwise_Product( this->variable_num, this->func[if1], this->func[if1], func_2);
		const double norm = sqrt(cal_norm(func_2));
		delete[]func_2; func_2=nullptr;
		
		// en = gn / ||gn||
		dscal( this->variable_num, 1.0/norm, this->func[if1], 1);
	}
	return func;
}

double Gram_Schmidt_Orth::cal_norm(double *func_in)
{
	double norm(0.0);
	switch(coordinate)
	{
		case 0:
			Mathzone::Simpson_Integral( this->variable_num, func_in, this->rab, norm);		
			break;
		case 1:	
			double *tmp_func = new double[this->variable_num];
			Mathzone::Pointwise_Product( this->variable_num, func_in, this->r_2, tmp_func);
			Mathzone::Simpson_Integral( this->variable_num, tmp_func, this->rab, norm);	
			delete[]tmp_func; tmp_func=nullptr;
			break;
	}
	return norm;
}
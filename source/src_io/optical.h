#ifndef OPTICAL_H
#define OPTICAL_H

//This seems to be calculating epsilon for optical excitations
//But is not currently in use

class Optical
{
	public:

	Optical();
	~Optical();
	
	static bool opt_epsilon2;
	static int  opt_nbands;
	
	void cal_epsilon2(const int &nbands);

	private:
	double element_cvk(const int &ik, const int &iv, const int &ic);


};

#endif

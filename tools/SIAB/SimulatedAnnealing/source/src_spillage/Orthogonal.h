#ifndef ORTHOGONAL_H
#define ORTHOGONAL_H

#include "common.h"
#include "SpillageStep.h"

class Orthogonal
{
	public:
	Orthogonal();
	~Orthogonal();

	void start(SpillageStep &step1, SpillageStep &step2);

	int test;
};

#endif


#ifndef ORTHOGONAL_H
#define ORTHOGONAL_H

#include "SpillageStep.h"
#include "common.h"

class Orthogonal
{
  public:
    Orthogonal();
    ~Orthogonal();

    void start(SpillageStep& step1, SpillageStep& step2);

    int test;
};

#endif

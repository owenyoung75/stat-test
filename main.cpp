#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

#include "basic_funcs.hpp"

typedef std::vector<double> Vector;

using namespace STAT_TEST;

int main()
{
    Vector example {6.0,10.0,7.0,9.0,8.0,12.0,0.4};

    double mean = Variance<double> (example);
    
    std::cout << mean << std::endl;
    
    return 0;
}

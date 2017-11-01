//
//  distribution.cpp
//  
//
//  Created by YangTong on 11/1/17.
//
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <functional>

#include <math.h>
#include <vector>
#include <numeric>
#include <random>

#include "distribution.hpp"

int main()
{
    float mean = 0.5;
    float var = 1.0/12.0;
    
//    auto f = std::bind(&RAND_D::Uniform_p, std::placeholders::_1, mean, var );

    
    RAND_D::d_func Uniform;
        Uniform.mean = mean;
        Uniform.var = var;
        Uniform.pdf = std::bind(&RAND_D::Uniform_p, std::placeholders::_1, mean, var );
        Uniform.cdf = std::bind(&RAND_D::Uniform_c, std::placeholders::_1, mean, var );

    
    
    printf("%f\n", Uniform.cdf(0.8) );
    
    return 0;
}




//
//  distribution.hpp
//
//
//  Created by YangTong on 10/31/17.
//
//
//  This file defines a general random distribution characterized by
//  * param:  mean
//  * param:  var
//  * param:  PDF
//  * param:  CDF


#ifndef _DISTR_FUNC_H_
#define _DISTR_FUNC_H_


#include <stdlib.h>
#include <functional>

#include <math.h>
#include <numeric>
#include <random>

typedef std::function<float (float)> float_pointer;



namespace RAND_D
{
    class d_func
	{
    public:
		float mean;
		float var;
        float_pointer pdf;
		float_pointer cdf;
	};
	


    
    
	// uniform distribution
	inline float Uniform_p(float x, float mean, float var)
	{		
		float up_bound = ( 2*mean + sqrt(var*12) )/2.0;
		float low_bound = ( 2*mean - sqrt(var*12) )/2.0;
		float range = up_bound - low_bound;
		if (x >= low_bound && x <= up_bound)
			return 1.0/range;
		else
			return 0.0;
	}
	inline float Uniform_c(float x, float mean, float var)
	{
		float up_bound = ( 2*mean + sqrt(var*12) )/2.0;
		float low_bound = ( 2*mean - sqrt(var*12) )/2.0;
		float range = up_bound - low_bound;
		if (x < low_bound)
			return 0.0;
		else if (x > up_bound)
			return 1.0;
		else
			return (x - low_bound)/range;			
	}
	

}

#endif

//
//  ED.hpp
//  
//
//  Created by YangTong on 10/31/17.
//
//
//  This file defines empirical distribution function (EDF) from a given training-set
//  and then calculate EDF values for a given testing-set


#ifndef ED_hpp
#define ED_hpp

#include <fstream>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <numeric>
#include <random>


typedef std::vector<float> Vector;


namespace ED
{
      
    template<class T1, class T2> inline float EDF(T1 y, T2 X)
    {
        y = float(y);
        float n = 0;
        int train_set_size = X.size();
        float *xd = X.data();
        
        while ( *xd++ <= y && n < train_set_size) n += 1.0;
        return n/train_set_size;
    }

    template<class vec> inline vec Z(vec Y)
    {
        Vector Z;
        int sz = Y.size();
        float *yd = Y.data();
        while (--sz >= 0) Z.push_back(EDF<vec, Vector>(*yd++, X));
        return Z;
    }

}


#endif /* ED_hpp */

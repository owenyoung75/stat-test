//
//  ED.cpp
//  
//
//  Created by YangTong on 10/31/17.
//
//
//  This file produce training data set for EDF calculation
//  training data, which is not needed anywhere else, can be read from a MC sampling, or a test sample

#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "ED.hpp"

typedef std::vector<float> Vector;

void obtain_EDF()
{
    // read some data here, or use the example below
    
    inline Vector data_example()
    {
        static const float arr[] = {6.0,1.0,77.0,29.0};
        Vector example (arr, arr + sizeof(arr) / sizeof(arr[0]) );
        sort(example.begin(), example.end());
    }
    
    Vector X = data_example();
    
}

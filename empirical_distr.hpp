//
//  empirical_distr.hpp
//  
//
//  Created by YangTong on 10/31/17.
//
//
//  This file defines empirical distribution function (EDF) from a given training-set
//  and then calculate EDF values for a given testing-set
#pragma once

#ifndef empirical_dsitr_hpp
#define empirical_dsitr_hpp

#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

#include "basic_funcs.hpp"

namespace STAT_TEST     
{
    // A general structure of an EMPIRICAL-DISTRIBUTION
    template<typename _num_type>
    struct Empirical_Distribution
    {
        // empirical-data
        std::vector<_num_type> Data_sample;
        
        // empirical-data features
        int Data_size;
        _num_type Min;
        _num_type Max;
        double Mean;
        double Variance;
        
        // empirical CDF for a single numer
        inline double cumulative_distr_func (_num_type _data_test);
        
        // empirical CDF for a vector
        inline std::vector<double> cumulative_distr_func (std::vector<_num_type> _data_test);
   
        
        
        // Constructor of an EMPIRICAL-DISTRIBUTION from TRAINING DATA
        Empirical_Distribution (std::vector<_num_type> _data_sample);
        
    };

    
}

#include "empirical_distr.ipp"

#endif /* empirical_distr_hpp */

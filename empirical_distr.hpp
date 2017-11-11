//
//  empirical_dsitr.hpp
//  
//
//  Created by YangTong on 10/31/17.
//
//
//  This file defines empirical distribution function (EDF) from a given training-set
//  and then calculate EDF values for a given testing-set


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
    // example data group for test
    std::vector<double> Data_Example_1 {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};

    
    
    // A general structure for an EMPIRICAL-DISTRIBUTION
    template<typename _num_type = double>
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
        inline double
        cumulative_distr_func (_num_type _data_test)
        {
            double cdf = 0.0;
            _num_type *xd = Data_sample.data();
            
            while ( *xd ++ <= _data_test  &&  cdf < Data_size) cdf += 1.0;
            return cdf/Data_size;
        }

        // empirical CDF for a vector
        inline std::vector<double>
        cumulative_distr_func (std::vector<_num_type> _data_test)
        {
            std::vector<double> cdf_vector;
            int _test_size = _data_test.size();
            double *yd = _data_test.data();
            
            while (--_test_size >= 0)
            {
                cdf_vector.push_back( cumulative_distr_func(*yd++) );
            }
            return cdf_vector;
        }
   
        
        
        // Constructor of an EMPIRICAL-DISTRIBUTION from TRAINING DATA
        Empirical_Distribution (std::vector<_num_type> _data_sample)
        {
            Data_sample = _data_sample;
            Data_size = Data_sample.size();
            sort(Data_sample.begin(), Data_sample.end());
            Min = Data_sample[0];
            Max = Data_sample[Data_size-1];
            Mean = Mean_value<_num_type>(Data_sample);
            Variance = Variance_value<_num_type>(Data_sample);
        }
        
        // Default constructor of EMPIRICAL-DISTRIBUTION for TEST
        Empirical_Distribution ()
        {
            Data_sample = Data_Example_1;
            Data_size = Data_sample.size();
            sort(Data_sample.begin(), Data_sample.end());
            Min = Data_sample[0];
            Max = Data_sample[Data_size-1];
            Mean = Mean_value<_num_type>(Data_sample);
            Variance = Variance_value<_num_type>(Data_sample);
        }
        
    };

    
}


#endif /* empirical_dsitr_hpp */

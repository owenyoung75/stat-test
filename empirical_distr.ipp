//
//  empirical_dsitr.ipp
//  
//
//  Created by YangTong on 1/3/18.
//
//
//  This file defines empirical distribution function (EDF) from a given training-set
//  and then calculate EDF values for a given testing-set
#pragma once

namespace STAT_TEST     
{
    // A general structure of an EMPIRICAL-DISTRIBUTION
    // empirical CDF for a single numer
    template<typename _num_type> inline double
    Empirical_Distribution<_num_type>::cumulative_distr_func (_num_type _data_test)
    {
        double cdf = 0.0;
        _num_type *xd = Data_sample.data();
        
        while ( *xd ++ <= _data_test  &&  cdf < Data_size) cdf += 1.0;
        return cdf/Data_size;
    }
    
    // empirical CDF for a vector
    template<typename _num_type> inline std::vector<double>
    Empirical_Distribution<_num_type>::cumulative_distr_func (std::vector<_num_type> _data_test)
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
    template<typename _num_type>
    Empirical_Distribution<_num_type>::Empirical_Distribution (std::vector<_num_type> _data_sample)
    {
        Data_sample = _data_sample;
        Data_size = Data_sample.size();
        sort(Data_sample.begin(), Data_sample.end());
        Min = Data_sample[0];
        Max = Data_sample[Data_size-1];
        Mean = Mean_value<_num_type>(Data_sample);
        Variance = Variance_value<_num_type>(Data_sample);
    }

    
}

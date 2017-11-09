//
//  empirical_dsitr.cpp
//  
//
//  Created by YangTong on 10/31/17.
//
//
//  This file produce training data set for EDF calculation
//  training data, which is not needed anywhere else, can be read from a MC sampling, or a test sample


#include "empirical_disitr.hpp"



namespace STAT_TEST {

    // a data-example for program testing
    inline std::vector<double> Data_Example_1()
    {
        static const double arr[] = {6.0,1.0,77.0,29.0};
        std::vector<double> example (arr, arr + sizeof(arr) / sizeof(arr[0]) );
        return example;
    }
    
    
    
    template<typename _num_type>
    Empirical_Distribution<_num_type>::Empirical_Distribution (std::vector<_num_type> _data_sample)
    {
        Data_sample = _data_sample;
        Data_size = Data_sample.size();
        sort(Data_sample.begin(), Data_sample.end());
        Min = Data_sample[0];
        Max = Data_sample[data_size-1];
        Mean = Mean_values<_num_type>(Data_sample);
        Variance = Variance<_num_type>(Data_sample);
    }
    
    
    template<typename _num_type>
    Empirical_Distribution<_num_type>::Empirical_Distribution ()
    {
        Data_sample = Data_Example_1();
        Data_size = Data_sample.size();
        sort(Data_sample.begin(), Data_sample.end());
        Min = Data_sample[0];
        Max = Data_sample[data_size-1];
        Mean = Mean_values<_num_type>(Data_sample);
        Variance = Variance<_num_type>(Data_sample);
    }

    
    
    
    template<typename _num_type>
    template<class _test_type> inline _test_type 
    Empirical_Distribution<_num_type>::cumulative_distr_func (_test_type _data_test)
    {
        if (typeid(_test_type) == typeid(_num_type))
        {
            double cdf = 0.0;
            _num_type *xd = Data_sample.data();
            
            while ( *xd ++ <= _data_test  &&  cdf < Data_size) cdf += 1.0;
            return cdf/Data_size;
        }
        else if (typeid(_test_type) == typeid(std::vector<_num_type>))
        {
            std::vector<double> cdf_vector;
            int _test_size = _data_test.size();
            _num_type *yd = _data_test.data();
            
            while (--_test_size >= 0)
            {
                double cdf = cumulative_distr_func<_num_type>(*yd++);
                cdf_vector.push_back(cdf);
            }
            return cdf_vector;
        }

    }
    


}

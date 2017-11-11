//
//  resampling.hpp
//  
//
//  Created by YangTong on 11/7/17.
//
//

#ifndef resampling_hpp
#define resampling_hpp

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include <stdio.h>



namespace STAT_TEST {
    // a data-example for program testing
    std::vector<double> Data_Example_0 {6.0,1.0,7.0,9.0,8.0,2.0,0.4,0.9};

 
    
    // FISHER-YATES-SHUFFLE process for the training data
    template<typename _num_type>
    inline std::vector<_num_type>
    Fisher_Yates_shuffle(std::vector<_num_type> _data_total, int _shuffle_times)
    {
        for (auto _times = 0; _times < _shuffle_times; _times++)
        {
            std::random_device rd;
            std::mt19937 mt(rd());
            int currentIndexCounter = _data_total.size();
            _num_type *iter = _data_total.data();
        
            while (--currentIndexCounter >= 0)
            {
                std::uniform_int_distribution<> dis(0, currentIndexCounter);
                const int randomIndex = dis(mt);
                if (*iter != _data_total.at(randomIndex))
                    std::swap(_data_total.at(randomIndex), *iter++);
            }
        }
    
        return _data_total;
    }
    
    
    
    
    // RESAMPLING process for a given group of data
    // --parameter-> _data_total        ::  training data
    // --parameter-> _subsample_size    ::  size of the pick-out sample
    // --parameter-> _resample_interval ::  number of shuffling for each resampling
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total,
               int _subsample_size,
               int _resample_interval
               )
    {
        if (_subsample_size > (_data_total.size())/2)
            _subsample_size = _data_total.size() - _subsample_size;
        
        // shuffle once
        _data_total = Fisher_Yates_shuffle(_data_total, _resample_interval);
        
        std::size_t const smaller_size = _subsample_size;

        std::vector<_num_type> smaller_part(_data_total.begin(), _data_total.begin() + smaller_size);
        std::vector<_num_type> larger_part(_data_total.begin() + smaller_size, _data_total.end());

        return std::make_tuple(smaller_part, larger_part);
    }
   
    // RESAMPLING process with
    // default resampling_interval = 1
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total,
               int _subsample_size
               )
    {
        int _resample_interval = 1;
        
        return resampling(_data_total, _subsample_size, _resample_interval);
    }
    
    // RESAMPLING process with
    // default resampling_interval = 1
    // default pick-out sample size = half of total size
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total)
    {
        int _subsample_size = _data_total.size()/2;
        return resampling(_data_total, _subsample_size);
    }
    
    // TEST RESAMPLING process with
    // testing data example
    // default resampling_interval = 1
    // default pick-out sample size = half of total size
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling()
    {
        std::vector<_num_type> _data_total = Data_Example_0;
        return resampling(_data_total);
    }
    
    
}



#endif /* resampling_hpp */

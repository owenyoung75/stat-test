//
//  resampling.cpp
//  
//
//  Created by YangTong on 11/7/17.
//
//

#include "resampling.hpp"

namespace STAT_TEST {
  
    // a data-example for program testing
    inline std::vector<double> Data_Example_0()
    {
        std::vector<double> example {6.0,10.0,7.0,9.0,8.0,12.0,0.4};
        return example;
    }
 
    
    
    
    template<typename _num_type>
    inline std::vector<std::vector<Num_type>>
    resampling(std::vector<_num_type> _data_total,
               int _single_sample_size)
    {
        _data_total = Fisher_Yates_shuffle(_data_total);

        if (_single_sample_size > (_data_total.size())/2)
            _single_sample_size = data_size - _single_sample_size;
        
        std::size_t const smaller_size = _single_sample_size;

        std::vector<Num_type> smaller_part(_data_total.begin(), _data_total.begin() + smaller_size);
        std::vector<Num_type> larger_part(_data_total.begin() + smaller_size, _data_total.end());

        return {smaller_part, larger_part};
    }
    
    
    template<typename _num_type>
    inline std::vector<std::vector<Num_type>>
    resampling(std::vector<_num_type> _data_total)
    {
        int smaller_size = _data_total.size()/2;
        return resampling(_data_total, smaller_size);
    }
    
    
    template<typename _num_type>
    inline std::vector<std::vector<Num_type>>
    resampling()
    {
        std::vector<_num_type> _data_total = Data_Example();
        return resampling(_data_total);
    }
    
    
    
    
    
    template<typename _num_type>
    inline std::vector<_num_type>
    Fisher_Yates_shuffle(std::vector<_num_type> _data_total)
    {
        std::random_device rd;
        std::mt19937 mt(rd());
        
//        std::cout << "Before: ";
//        std::copy(_data_total.cbegin(), _data_total.cend(),
//                  std::ostream_iterator<int>(std::cout, " "));
        
        auto currentIndexCounter = _data_total.size();
        
        for (auto iter = _data_total.rbegin(); iter != _data_total.rend();
             ++iter, --currentIndexCounter)
        {
            // get int distribution with new range
            std::uniform_int_distribution<> dis(0, currentIndexCounter);
            const int randomIndex = dis(mt);
            
            if (*iter != _data_total.at(randomIndex))
            {
                std::swap(_data_total.at(randomIndex), *iter);
            }
        }
        
//        std::cout << "\nAfter: ";
//        std::copy(_data_total.cbegin(), _data_total.cend(),
//                  std::ostream_iterator<int>(std::cout, " "));
        
        return _data_total;
    }
    
}

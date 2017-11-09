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

    template<typename Num_type = double>
    std::vector<std::vector<Num_type>> resampling();

    template<typename Num_type = double>
    std::vector<std::vector<Num_type>> resampling(std::vector<Num_type> Data_total);
    
    template<typename Num_type = double>
    std::vector<std::vector<Num_type>> resampling(std::vector<Num_type> Data_total,
                                                  int _Single_sample_size);

    
    
    template<typename Num_type = double>
    std::vector<Num_type> Fisher_Yates_shuffle(std::vector<Num_type> Data_total);
    
}



#endif /* resampling_hpp */

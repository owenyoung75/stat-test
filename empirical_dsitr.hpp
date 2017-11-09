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

#include "basic_functions.hpp"


namespace STAT_TEST     
{
    template<typename Num_type = double> struct Empirical_Distribution
    {
        Empirical_Distribution();
        Empirical_Distribution(std::vector<Num_type> _data_sample);
        
        std::vector<Num_type> Data_sample;
        
        int Data_size;
        Num_type Min;
        Num_type Max;
        double Mean;
        double Variance;
        
        template<class Test_type> Test_type cumulative_distr_func(Test_type data_test);
//        Num_type Sample_generator(int number_of_sampling);
    }

}


#endif /* empirical_dsitr_hpp */

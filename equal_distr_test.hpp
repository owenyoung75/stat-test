//
//  equal_distr_test.hpp
//  
//
//  Created by YangTong on 11/6/17.
//
//

#ifndef equal_distr_test_hpp
#define equal_distr_test_hpp


#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/legendre.hpp>

#include "empirical_distr.hpp"


namespace EQUAL_TEST {
    
    template<typename Num_type>
    std::vector<double> Test_vector_Psi(STAT_TEST::Empirical_Distribution<Num_type> _ED_X,
                                    std::vector<Num_type> Data_test,
                                    int Highest_order
                                    )
    
    double Psi_square(std::vector<double> Psi_vector)
    
    template<class T> inline T Normalized_Legendre_Poly(int n, T Z);
    
}


#endif /* equal_distr_test_hpp */

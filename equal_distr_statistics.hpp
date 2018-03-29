//
//  equal_distr_statistics.hpp
//  
//
//  Created by YangTong on 11/6/17.
//
//
#pragma once

#ifndef equal_test_statistics_hpp
#define equal_test_statistics_hpp


#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>

#include <boost/math/special_functions/legendre.hpp>

#include "empirical_distr.hpp"


namespace EQUAL_TEST {
    typedef std::vector<double> Vector;
    
    // NORMALIZED-LEGENDRE-POLYNOMIAL for a number
    inline double Normalized_Legendre_Poly (int _order, double _unit_variable);
    
    // NORMALIZED-LEGENDRE-POLYNOMIAL for a vector
    inline Vector Normalized_Legendre_Poly (int _order, Vector _unit_variable_vector);
    
    
    
    
    // creat the TESTED-VECTOR PI for a given empirical-distibution
    template<typename _num_type>
    inline Vector Test_vector_PI (STAT_TEST::Empirical_Distribution<_num_type> _ED_X,
                                  std::vector<_num_type> _data_test,
                                  int _highest_order
                                  );
    
    
    
    
    // MODULO-SQUARE of PSI boost::VECTOR
    inline double PSI_square (boost::numeric::ublas::vector<double> _pi_vector);
    
    // MODULO-SQUARE of PSI std::VECTOR
    inline double PSI_square (Vector _pi_vector);

}

#include "equal_distr_statistics.ipp"

#endif /* equal_distr_statistics_hpp */

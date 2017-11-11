//
//  equal_distr_test.h
//  
//
//  Created by YangTong on 11/1/17.
//
//
//  This is the file for data generation from a chosen distribution
#pragma once

#ifndef equal_distr_test_h
#define equal_distr_test_h

#include <fstream>
#include <stdlib.h>

#include <math.h>
#include <vector>
#include <numeric>
#include <random>

#include "equal_test_statistics.hpp"
#include "resampling.hpp"

typedef std::vector<double> Vector;

namespace STAT_TEST {
    namespace EQUAL_TEST{
        using namespace boost::numeric::ublas;

        // construct TESTING PI-VECTORS of a GIVEN SAMPLE
        template<typename _num_type>
        inline matrix<double>
        PI_package (std::vector<_num_type> _total_sample,
                                          int _size_of_subsample,
                                          int _num_of_resampling,
                                          int _highest_order_poly,
                                          int _resampling_interval,
                                          int _wamup_steps)
        {
            matrix<double> PI_package(_highest_order_poly, _num_of_resampling);
            std::tuple<Vector, Vector> resampled;
            Vector pi_vec;
            
            // warmup shuffling
            _total_sample = STAT_TEST::Fisher_Yates_shuffle(_total_sample, _wamup_steps);
            
            // get PI-vectors for distribution test
            for (int i = 0; i < _num_of_resampling; i++)
            {
                resampled = STAT_TEST::resampling<double>(_total_sample,
                                                          _size_of_subsample,
                                                          _resampling_interval);
                
                // Use large group to construct the empirical distribution
                STAT_TEST::Empirical_Distribution<_num_type> _emp_distr( std::get<1>(resampled) );
                
                // Use small group to construct the equal-distribution-test
                pi_vec = Test_vector_PI<_num_type>(_emp_distr,
                                                   std::get<0>(resampled),
                                                   _highest_order_poly);

                vector<double, Vector> pi(pi_vec);
                column(PI_package, i) = pi;
            }
            return PI_package;
        }
        
        // construct TESTING PI-VECTORS of a sample example
        template<typename _num_type = double>
        inline matrix<double> PI_package ()
        {
            Vector _total_sample {1.0, 3.0, 4.0, 4.3, 4.6, 4.9, 5.0, 5.1, 5.4, 5.7, 6.2, 7.0, 9.1};
            int _size_of_subsample = 4;
            int _num_of_resampling = 10;
            int _highest_order_poly = 3;
            int _resampling_interval = 1;
            int _wamup_steps = 10;
            
            return PI_package<double>(_total_sample,
                                      _size_of_subsample,
                                      _num_of_resampling,
                                      _highest_order_poly,
                                      _resampling_interval,
                                      _wamup_steps);
        }
        
    }
}



#endif /* equal_distr_test_h */

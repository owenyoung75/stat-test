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


#include "equal_distr_statistics.hpp"
#include "resampling.hpp"

typedef std::vector<double> Vector;

namespace STAT_TEST {
    namespace EQUAL_TEST{
        using namespace boost::numeric::ublas;

        
        // construct a PI-Psi-vector for a single testing sample
        template<typename _num_type> inline Vector Pi_and_Psi_tested(std::vector<_num_type> _training_sample,
                                                                     std::vector<_num_type> _testing_sample,
                                                                     int _highest_order_poly
                                                                     );

        // construct the PI-package for a TRAINING SAMPLE
        template<typename _num_type> inline matrix<double> PI_package (std::vector<_num_type> _total_sample,
                                                                       int _subsample_size,
                                                                       int _num_of_resampling,
                                                                       int _highest_order_poly,
                                                                       int _resampling_interval,
                                                                       int _wamup_steps
                                                                       );
        
        
        
        
        
        
        // Calculate critical VECTORs of each psi component and psi^2 for a training sample
        template<typename _num_type> inline std::tuple< matrix<double>, matrix<double> >
        Pi_and_Psi_ciritical_values(std::vector<_num_type> _training_sample,
                                    int    _subsample_size,
                                    int    _num_of_resampling,
                                    int    _highest_order_poly,
                                    int    _resampling_interval,
                                    int    _wamup_steps,
                                    Vector _critical_portions,
                                    bool   _plot_psi_square_distr,
                                    std::string _file_name
                                    );

        
        // Calculate critical VALUEs of each psi component and psi^2 for a training sample
        template<typename _num_type> inline std::tuple< Vector, Vector >
        Pi_and_Psi_ciritical_values(std::vector<_num_type> _training_sample,
                                    int    _subsample_size,
                                    int    _num_of_resampling,
                                    int    _highest_order_poly,
                                    int    _resampling_interval,
                                    int    _wamup_steps,
                                    double _critical_portion,
                                    bool   _plot_psi_square_distr,
                                    std::string _file_name
                                    );
        
        
        
        
        // For un-plotting training cases
        template<typename _num_type> inline std::tuple< Vector, Vector >
        Pi_and_Psi_ciritical_values(std::vector<_num_type> _training_sample,
                                    int    _subsample_size,
                                    int    _num_of_resampling,
                                    int    _highest_order_poly,
                                    int    _resampling_interval,
                                    int    _wamup_steps,
                                    double _critical_portion
                                    );

        
        
    }
}

#include "equal_distr_test.ipp"

#endif /* equal_distr_test_h */

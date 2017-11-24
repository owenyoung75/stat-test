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

        // construct TESTING PI-VECTORS for a GIVEN SAMPLE
        template<typename _num_type>
        inline matrix<double>
        PI_package (std::vector<_num_type> _total_sample,
                    int _subsample_size,
                    int _num_of_resampling,
                    int _highest_order_poly,
                    int _resampling_interval,
                    int _wamup_steps)
        {
            matrix<double> PI_package(_highest_order_poly + 1, _num_of_resampling);
            std::tuple<Vector, Vector> resampled;
            Vector pi_vec;
            
            // warmup shuffling
            // _total_sample = STAT_TEST::Fisher_Yates_shuffle(_total_sample, _wamup_steps);
            _total_sample = STAT_TEST::standard_shuffle(_total_sample, _wamup_steps);
            
            // get PI-vectors for distribution test
            for (int i = 0; i < _num_of_resampling; i++)
            {
                resampled = STAT_TEST::resampling<double>(_total_sample,
                                                          _subsample_size,
                                                          _resampling_interval);
                
                // Use large group to construct the empirical distribution
                STAT_TEST::Empirical_Distribution<_num_type> _emp_distr( std::get<1>(resampled) );
                
                // Use small group to construct the equal-distribution-test
                pi_vec = Test_vector_PI<_num_type>(_emp_distr,
                                                   std::get<0>(resampled),
                                                   _highest_order_poly);
                pi_vec.push_back(PSI_square(pi_vec));

                vector<double, Vector> pi(pi_vec);
                column(PI_package, i) = pi;
            }
            return PI_package;
        }
        
        
        

        template<typename _num_type>
        inline std::tuple< Vector, Vector >
        Pi_and_Psi_ciritical_values(std::vector<_num_type> _training_sample,
                                    int    _subsample_size,
                                    int    _num_of_resampling,
                                    int    _highest_order_poly,
                                    int    _resampling_interval,
                                    int    _wamup_steps,
                                    double _critical_portion,
                                    bool   _plot_psi_square_distr
                                    )
        {
            matrix<double> PI = PI_package<double> (_training_sample,
                                                    _subsample_size,
                                                    _num_of_resampling,
                                                    _highest_order_poly,
                                                    _resampling_interval,
                                                    _wamup_steps);
            
            int critical_index = int(_training_sample.size() * (1-_critical_portion) / 2 );
            
            // construct two-sided test boundary
            Vector lower_bound;
            Vector upper_bound;
            for (int i = 0; i <= _highest_order_poly; i++)
            {
                vector<double> pi_i = row(PI, i);
                boost::sort(pi_i);
                lower_bound.push_back(pi_i(critical_index));
                upper_bound.push_back(pi_i(_num_of_resampling - critical_index - 1));
            }
            
            
            // Plot the Psi^2 distribution shape
            if (_plot_psi_square_distr)
            {
                std::ofstream ofile("trained_normal.txt",std::ios::out);
                int rescale = int(_num_of_resampling/200);
                vector<double> Psi = row(PI, _highest_order_poly);
                std::map<int, int> hist;
                for(int n=0; n<_num_of_resampling; ++n)
                    ++hist[std::round(Psi(n))];
                for(auto p : hist)
                {
                    ofile << std::fixed << std::setprecision(1) << std::setw(2)
                    << p.first << ' ' << std::string(p.second/rescale , '*') << '\n';
                }
            }
            
            return std::make_tuple(lower_bound, upper_bound);

        }

        
        
        
        template<typename _num_type>
        inline Vector Pi_and_Psi_tested(std::vector<_num_type> _training_sample,
                                        std::vector<_num_type> _testing_sample,
                                        int _highest_order_poly
                                        )
        {
            STAT_TEST::Empirical_Distribution<_num_type> ED_training (_training_sample);
            Vector pi_test = Test_vector_PI<_num_type> (ED_training,
                                                        _testing_sample,
                                                        _highest_order_poly
                                                        );
            pi_test.push_back(PSI_square(pi_test));
            return pi_test;
        }
        

        
    }
}



#endif /* equal_distr_test_h */

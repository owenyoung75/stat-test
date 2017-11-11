//
//  equal_test_statistics.hpp
//  
//
//  Created by YangTong on 11/6/17.
//
//

#ifndef equal_test_statistics_hpp
#define equal_test_statistics_hpp


#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/legendre.hpp>

#include "empirical_distr.hpp"


namespace STAT_TEST {
    namespace EQUAL_TEST {
        
        // NORMALIZED-LEGENDRE-POLYNOMIAL for a number
        inline double Normalized_Legendre_Poly
        (int _order, double _cdf)
        {
            double Pi = boost::math::legendre_p<double>(_order, _cdf);
            Pi *= sqrt(2*_order + 1);
            return Pi;
        }
        
        // NORMALIZED-LEGENDRE-POLYNOMIAL for a vector
        inline std::vector<double> Normalized_Legendre_Poly
        (int _order, std::vector<double> _cdf_vector)
        {
            std::vector<double> Pi;
            int _cdf_size = _cdf_vector.size();
            double *cdf_p = _cdf_vector.data();
            
            while (--_cdf_size >= 0)
            {
                Pi.push_back( Normalized_Legendre_Poly(_order, *cdf_p++) );
            }
            
            return Pi;
        }
        
        
        
        
        // creat the TESTED-VECTOR PI for a given empirical-distibution
        template<typename _num_type>
        inline std::vector<double>
        Test_vector_PI(STAT_TEST::Empirical_Distribution<_num_type> _ED_X,
                       std::vector<_num_type> _data_test,
                       int _highest_order
                       )
        {
            std::vector<double> cdf_vector = _ED_X.cumulative_distr_func(_data_test);
            std::vector<double> PI_vector;
            
            for (int i=1; i <= _highest_order; i++)
            {
                std::vector<double> pi_i = Normalized_Legendre_Poly(i, cdf_vector);
                double Pi = std::accumulate(pi_i.begin(), pi_i.end(), 0.0) / sqrt(cdf_vector.size());
                PI_vector.push_back(Pi);
            }
            
            return PI_vector;
        };
        
        
        
        
        // MODULO-SQUARE of PSI boost::VECTOR
        inline double
        PSI_square (boost::numeric::ublas::vector<double> _pi_vector)
        {
            double psi_square = 0.0;
            for (int i = 0; i < _pi_vector.size(); i++)
                psi_square += _pi_vector(i) * _pi_vector(i);
            
            return psi_square;
        }
        
        // MODULO-SQUARE of PSI std::VECTOR
        inline double
        PSI_square (std::vector<double> _pi_vector)
        {
            double psi_square = 0.0;
            int _order = _pi_vector.size();
            double *pi = _pi_vector.data();
            
            while (--_order >= 0)
            {
                psi_square += (*pi) * (*pi);
                pi += 1;
            }
            return psi_square;
        }
        
    }
}



#endif /* equal_test_statistics_hpp */

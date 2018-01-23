//
//  equal_test_statistics.hpp
//  
//
//  Created by YangTong on 11/6/17.
//
//
#pragma once

typedef std::vector<double> Vector;


namespace STAT_TEST {
    namespace EQUAL_TEST {
        
        // NORMALIZED-LEGENDRE-POLYNOMIAL for a number
        inline double Normalized_Legendre_Poly (int _order, double _unit_variable)
        {
            double Pi = boost::math::legendre_p<double>(_order, (2*_unit_variable-1) );
            Pi *= sqrt(2*_order + 1);
            return Pi;
        }
        
        // NORMALIZED-LEGENDRE-POLYNOMIAL for a vector
        inline Vector Normalized_Legendre_Poly (int _order,
                                                Vector _unit_variable_vector)
        {
            Vector Pi;
            int _unit_variable_size = _unit_variable_vector.size();
            double *unit_variable_p = _unit_variable_vector.data();
            
            while (--_unit_variable_size >= 0)
                Pi.push_back( Normalized_Legendre_Poly(_order, *unit_variable_p++) );
            
            return Pi;
        }
        
        
        
        
        // creat the TESTED-VECTOR PI for a given empirical-distibution
        template<typename _num_type>
        inline Vector Test_vector_PI (STAT_TEST::Empirical_Distribution<_num_type> _ED_X,
                                      std::vector<_num_type> _data_test,
                                      int _highest_order
                                      )
        
        {
            Vector cdf_vector = _ED_X.cumulative_distr_func(_data_test);
            Vector PI_vector;
            
            for (int i=1; i <= _highest_order; i++)
            {
                Vector pi_i = Normalized_Legendre_Poly(i, cdf_vector);
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
        PSI_square (Vector _pi_vector)
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


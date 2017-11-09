//
//  equal_distr_test.cpp
//  
//
//  Created by YangTong on 11/6/17.
//
//

#include "equal_distr_test.hpp"


namespace EQUAL_TEST {
    
    // Obtain the to-be-tested Psi vector for general cases
    // _num_type usually can eigher be integer or double
    template<typename _num_type>
    inline std::vector<double>
    Test_vector_Psi(STAT_TEST::Empirical_Distribution<_num_type> _ED_X,
                    std::vector<_num_type> _data_test,
                    int _highest_order
                    )
    {
        std::vector<double> Z_cdf_vector = _ED_X.cumulative_distr_func<std::vector<_num_type>>(_data_test);
        
        std::vector<double> Psi_vector;
        
        for (int i=0; i<_highest_order; i++)
        {
            std::vector<double> pi_i = Normalized_Legendre_Poly(i, Z_cdf_vector);
            double Psi_i = std::accumulate(pi_i.begin(), pi_i.end(), 0.0);
            Psi_vector.push_back(Psi_i);
        }
        
        return Psi_vector;
    }
    
    // Obtain the to-be-tested Psi^2 after obtain Psi
    // Note this is not asymptotically chi^2 distributed for case in Theorem-3
    inline double Psi_square(std::vector<double> _psi_vector)
    {
        return std::accumulate(_psi_vector.begin(), _psi_vector.end(), 0.0);
    }
    
    
    
    
    // Obtain normalized Legendre polynomials values for a vector Z
    template<class _cdf_tpe>
    inline _cdf_type
    Normalized_Legendre_Poly(int _order, _cdf_type _cdf)
    {
        _cdf_type Pi = boost::math::legendre_p<_cdf_type>(_order, _cdf);
        
        if (typeid(_cdf_type) == typeid(double))
            Pi *= sqrt(2*_order + 1);
        else if (typeid(_cdf_type) == typeid(std::vector<double>))
            transform(Pi.begin(), Pi.end(), Pi.begin(), _1 * sqrt(2*_order + 1));
        
        return Pi;
    }
    
}

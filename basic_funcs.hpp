//
//  basic_funcs.hpp
//  
//
//  Created by YangTong on 11/5/17.
//
//
#pragma once

#include <vector>
#include <stdio.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#ifndef basic_functions_hpp
#define basic_functions_hpp


namespace STAT_TEST {
    using namespace boost::numeric::ublas;
    
    //////////////////////////////////////////////////////////////////////////
    ///////////////////   Statistical Constants setting  /////////////////////
    //////////////////////////////////////////////////////////////////////////
    // Mathematical constants
    constexpr double PI     = boost::math::constants::pi<double>();
    constexpr double E      = boost::math::constants::e<double>();
    // Critical values
    constexpr double Nv0    = 0.999;
    constexpr double Nv1    = 0.995;
    constexpr double Nv2    = 0.990;
    constexpr double Nv3    = 0.975;
    constexpr double Nv4    = 0.950;
    constexpr double Nv5    = 0.900;
    
    
    

    //////////////////////////////////////////////////////////////////////////
    //////////////////   Statistical Characteristics   ///////////////////////
    //////////////////////////////////////////////////////////////////////////
    // calculate MEAN of a std::vector
    template<typename _num_type>
    inline double Mean_value(std::vector<_num_type> _data)
    {
        double rlt = accumulate( _data.begin(), _data.end(), 0.0)/ _data.size();
        return rlt;
    }
    // calculate MEAN of a boost::vector
    template<typename _num_type>
    inline double Mean_value(vector<_num_type> _data)
    {
        double rlt = accumulate( _data.begin(), _data.end(), 0.0)/ _data.size();
        return rlt;
    }

    
    // calculate VARIANCE of a std::vector  --> sigma^2
    template<typename _num_type>
    inline double Variance_value(std::vector<_num_type> _data_group)
    {
        int size_group = _data_group.size();
        double mean = Mean_value<_num_type>(_data_group);
        
        double sum = 0.0;
        for (int i = 0; i<size_group; i++)
            sum += (_data_group[i] - mean)*(_data_group[i] - mean);
        
        return sum/size_group;
    }
    // calculate VARIANCE of a boost::vector --> sigma^2
    template<typename _num_type>
    inline double Variance_value(vector<_num_type> _data_group)
    {
        int size_group = _data_group.size();
        double mean = Mean_value<_num_type>(_data_group);
        
        double sum = 0.0;
        for (int i = 0; i<size_group; i++)
            sum += (_data_group(i) - mean)*(_data_group(i) - mean);
        
        return sum/size_group;
    }
    

    // calculate COVARIANCE of a rectagular boost::matrix
    // _data_groups = (row_idx, col_idx)
    // row_idx :: variable 1, 2, 3 ......
    // col_idx :: copy     1, 2, 3 ...... --> have closer address
    template<typename _num_type>
    inline matrix<double> Covariance_Matrix (matrix<_num_type> _data_groups)
    {
        int num_of_varbs = _data_groups.size1();
        vector<double> means(num_of_varbs);
        matrix<double> Matx(num_of_varbs, num_of_varbs);
        
        for (int i = 0; i<num_of_varbs; i++)
        {
            means(i) = Mean_value<_num_type>( row(_data_groups, i) );
            Matx(i,i) = Variance_value<_num_type>( row(_data_groups, i) );
        }

        for (int i = 0; i<(num_of_varbs-1); i++)
        {
            for (int j = i+1; j<num_of_varbs; j++)
            {
                double sum = 0.0;
                for (int cpy = 0; cpy<_data_groups.size2(); cpy++)
                    sum += ( _data_groups(i,cpy) - means(i) ) * ( _data_groups(j,cpy) - means(j) );
                
                Matx(i,j) = sum/(_data_groups.size2());
                Matx(j,i) = Matx(i,j);
            }
        }
        
        return Matx;
    }
   
    
    
    
    
    //////////////////////////////////////////////////////////////////////////
    ////////////////////////   Normal Distributions   ////////////////////////
    //////////////////////////////////////////////////////////////////////////
    // calculate normal-distr PDF value for given parameters:
    // _x_scalar:: scalar-type input variable
    // _mean    :: mean value of distr
    // _variance:: variance value of distr
    inline double pdf_normal_distr(double _x_scalar,
                                   double _mean,
                                   double _variance)
    {
        return exp( -(_x_scalar - _mean)*(_x_scalar - _mean)/(2.0 * _variance) )/ sqrt(2 * PI * _variance);
    }

    // _x_vector::  std::vector-type input variable
    inline std::vector<double> pdf_normal_distr(std::vector<double> _x_vector,
                                                double _mean,
                                                double _variance)
    {
        std::vector<double> _pdf_vector;
        int _x_size = _x_vector.size();
        double *xd = _x_vector.data();
        
        while (--_x_size >= 0)
        {
            _pdf_vector.push_back( pdf_normal_distr(*xd++, _mean, _variance) );
        }
        return _pdf_vector;
    }
    // _x_vector::  boost::vector-type input variable
    inline vector<double> pdf_normal_distr( vector<double> _x_vector,
                                            double _mean,
                                            double _variance)
    {
        int _x_size = _x_vector.size();
        vector<double> _pdf_vector(_x_size);
        for (int i = 0; i<_x_size; i++)
        {
            _pdf_vector(i) = pdf_normal_distr ( _x_vector(i), _mean, _variance );
        }

        return _pdf_vector;
    }
    
    //calculate normal-distr PDF value for given distr
    // _x_scalr::  scalar-type input variable
    inline double pdf_normal_distr(double _x_scalar,
                                   std::normal_distribution<> _normal_distr)
    {
        return pdf_normal_distr( _x_scalar,
                                 _normal_distr.mean(),
                                 _normal_distr.stddev() * _normal_distr.stddev()
                               );
    }
    // _x_vector::  std::vector-type input variable
    inline std::vector<double>  pdf_normal_distr(std::vector<double> _x_vector,
                                                 std::normal_distribution<> _normal_distr)
    {
        std::vector<double> _pdf_vector;
        int _x_size = _x_vector.size();
        double *xd = _x_vector.data();
        
        while (--_x_size >= 0)
        {
            _pdf_vector.push_back(pdf_normal_distr( *xd++,
                                                    _normal_distr.mean(),
                                                    _normal_distr.stddev() * _normal_distr.stddev()
                                                  )
                                 );
        }
        return _pdf_vector;
    }
    // _x_vector::  boost::vector-type input variable
    inline vector<double>  pdf_normal_distr(vector<double> _x_vector,
                                            std::normal_distribution<> _normal_distr)
    {
        int _x_size = _x_vector.size();
        vector<double> _pdf_vector(_x_size);
        for (int i = 0; i<_x_size; i++)
        {
            _pdf_vector(i) = pdf_normal_distr( _x_vector(i),
                                               _normal_distr.mean(),
                                               _normal_distr.stddev() * _normal_distr.stddev()
                                              );
        }
        
        return _pdf_vector;
    }

    
    // calculate normal-distr cDF value for given parameters:
    // _x_scalar:: scalar-type input variable
    // _mean    :: mean value of distr
    // _variance:: variance value of distr
    inline double cdf_normal_distr(double _x_scalar,
                                   double _mean,
                                   double _variance)
    {
        return (1 + std::erf( (_x_scalar - _mean)/sqrt(_variance * 2.0) ) ) / 2.0;
    }
    
    // _x_vector::  std::vector-type input variable
    inline std::vector<double> cdf_normal_distr(std::vector<double> _x_vector,
                                                double _mean,
                                                double _variance)
    {
        std::vector<double> _cdf_vector;
        int _x_size = _x_vector.size();
        double *xd = _x_vector.data();
        
        while (--_x_size >= 0)
        {
            _cdf_vector.push_back( cdf_normal_distr(*xd++, _mean, _variance) );
        }
        return _cdf_vector;
    }
    // _x_vector::  boost::vector-type input variable
    inline vector<double> cdf_normal_distr( vector<double> _x_vector,
                                           double _mean,
                                           double _variance)
    {
        int _x_size = _x_vector.size();
        vector<double> _cdf_vector(_x_size);
        for (int i = 0; i<_x_size; i++)
        {
            _cdf_vector(i) = cdf_normal_distr ( _x_vector(i), _mean, _variance );
        }
        
        return _cdf_vector;
    }
    
    //calculate normal-distr PDF value for given distr
    // _x_scalr::  scalar-type input variable
    inline double cdf_normal_distr(double _x_scalar,
                                   std::normal_distribution<> _normal_distr)
    {
        return cdf_normal_distr( _x_scalar,
                                 _normal_distr.mean(),
                                 _normal_distr.stddev() * _normal_distr.stddev()
                                );
    }
    // _x_vector::  std::vector-type input variable
    inline std::vector<double>  cdf_normal_distr(std::vector<double> _x_vector,
                                                 std::normal_distribution<> _normal_distr)
    {
        std::vector<double> _cdf_vector;
        int _x_size = _x_vector.size();
        double *xd = _x_vector.data();
        
        while (--_x_size >= 0)
        {
            _cdf_vector.push_back(cdf_normal_distr( *xd++,
                                                    _normal_distr.mean(),
                                                    _normal_distr.stddev() * _normal_distr.stddev()
                                                   )
                                  );
        }
        return _cdf_vector;
    }
    // _x_vector::  boost::vector-type input variable
    inline vector<double>  cdf_normal_distr(vector<double> _x_vector,
                                            std::normal_distribution<> _normal_distr)
    {
        int _x_size = _x_vector.size();
        vector<double> _cdf_vector(_x_size);
        for (int i = 0; i<_x_size; i++)
        {
            _cdf_vector(i) = cdf_normal_distr( _x_vector(i),
                                               _normal_distr.mean(),
                                               _normal_distr.stddev() * _normal_distr.stddev()
                                              );
        }
        
        return _cdf_vector;
    }
    
    
}


#endif /* basic_funcs_hpp */










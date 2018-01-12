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
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#ifndef basic_funcs_hpp
#define basic_funcs_hpp


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
    constexpr double Nv5    = 0.9;
    
    
    

    //////////////////////////////////////////////////////////////////////////
    //////////////////   Statistical Characteristics   ///////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    // calculate MEAN of a std::vector
    template<typename _num_type>
    inline double Mean_value(std::vector<_num_type> _data);
    
    // calculate MEAN of a boost::vector
    template<typename _num_type>
    inline double Mean_value(vector<_num_type> _data);

    
    // calculate VARIANCE of a std::vector  --> sigma^2
    template<typename _num_type>
    inline double Variance_value(std::vector<_num_type> _data_group);
    
    // calculate VARIANCE of a boost::vector --> sigma^2
    template<typename _num_type>
    inline double Variance_value(vector<_num_type> _data_group);
    

    // calculate COVARIANCE of a rectagular boost::matrix
    // _data_groups = (row_idx, col_idx)
    // row_idx :: variable 1, 2, 3 ......
    // col_idx :: copy     1, 2, 3 ...... --> have closer address
    template<typename _num_type>
    inline matrix<double> Covariance_Matrix (matrix<_num_type> _data_groups);
   
    
    
    // obtain the coarse-grained histgram
    template<typename _num_type>
    inline std::map<int, int> coars_grained_distr (std::vector<_num_type> _data);
    
    
    
    // Plot the distribution
    template<typename _num_type>
    inline void Plot_distribution (std::vector<_num_type> _data);
    
    
}

#include "basic_funcs.ipp"

#endif /* basic_funcs_hpp */










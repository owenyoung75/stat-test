//
//  basic_funcs.hpp
//  
//
//  Created by YangTong on 11/5/17.
//
//
#pragma once

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
#include <random>
#include <string>
#include <iomanip>
#include <map>
#include <cmath>
#include <chrono>
#include <thread>
#include <unistd.h>


#include <boost/array.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/math/constants/constants.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>


#ifndef basic_funcs_hpp
#define basic_funcs_hpp


namespace STAT_TEST {
    using namespace boost::numeric::ublas;
    typedef vector<double>  Vector;
    typedef matrix<double>  Matrix;
    
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
    // return:: matrix(row_idx, row_idx)
    // row_idx :: # of variables  1, 2, 3 ......
    // col_idx :: # of copies     1, 2, 3 ...... --> have closer address
    template<typename _num_type>
    inline Matrix Covariance_Matrix (matrix<_num_type> _data_groups);
   
    
    // Calculate inverse matrix by LU decomp
    template<class T> bool InvertMatrix (const matrix<T>& input, matrix<T>& inverse);
    
    
    
    // PDF function of Normal distribution
    inline double NormalPDF (double x, double mu, double var);
    inline Vector NormalPDF (Vector X, double mu, double var);
    inline double NormalPDF (double x);
    inline Vector NormalPDF (Vector X);
    
    
    
    // obtain the coarse-grained histgram
    template<typename _num_type>
    inline std::map<int, int> coarse_grained_distr (std::vector<_num_type> _data);
    
    
    
    // Plot the distribution
    template<typename _num_type>
    inline void Plot_distribution (std::vector<_num_type> _data);
    
    
}

#include "basic_funcs.ipp"

#endif /* basic_funcs_hpp */










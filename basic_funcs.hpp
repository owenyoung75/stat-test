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

#ifndef basic_functions_hpp
#define basic_functions_hpp


namespace STAT_TEST {
    using namespace boost::numeric::ublas;

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
 
    
    
    
    // calculate VARIANCE of a std::vector
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
    // calculate VARIANCE of a boost::vector
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
    // col_idx :: copy 1, 2, 3 ...... --> closer address
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
   
    
    
}


#endif /* basic_funcs_hpp */










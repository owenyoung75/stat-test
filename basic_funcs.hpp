//
//  basic_funcs.hpp
//  
//
//  Created by YangTong on 11/5/17.
//
//
#include <vector>
#include <stdio.h>

#ifndef basic_functions_hpp
#define basic_functions_hpp


namespace STAT_TEST {
    
    template<typename _num_type>
    inline double Mean_value(std::vector<_num_type> _data_group)
    {
        int size_group = _data_group.size();
        
        double sum = 0.0;
        for (int i = 0; i<size_group; i++)
            sum += _data_group[i];
        
        return sum/size_group;
    }
 
    
    
    template<typename _num_type>
    inline double Variance(std::vector<_num_type> _data_group)
    {
        int size_group = _data_group.size();
        double mean = Mean_value<_num_type>(_data_group);
        
        double sum = 0.0;
        for (int i = 0; i<size_group; i++)
            sum += (_data_group[i] - mean)*(_data_group[i] - mean);
        
        return sum/size_group;
    }
    
    
    
    template<typename _num_type>
    inline std::vector<std::vector<double>> Covariance_Matrix
    (std::vector<std::vector<_num_type>> _data_groups)
    {
        int num_of_varbs = _data_groups.size();
        std::vector<double> means;
        std::vector<std::vector<double>> Matx;
        
        for (int i = 0; i<num_of_varbs; i++)
        {
            means[i] = Mean_value<_num_type>(_data_groups[i]);
            Matx[i][i] = Variance<_num_type>(_data_groups[i]);
        }
        
        for (int i = 0; i<(num_of_varbs-1); i++)
        {
            for (int j = i+1; j<num_of_varbs; j++)
            {
                int size_group1 = _data_groups[i].size();
                int size_group2 = _data_groups[j].size();
                double sum = 0.0;
                
                for (int p = 0; p<size_group1; p++)
                {
                    for (int q = 0; q<size_group2; j++)
                        sum += (_data_groups[i] - means[i])*(_data_groups[j] - means[j]);
                }
                
                Matx[i][j] = sum/(size_group1*size_group2);
                Matx[j][i] = Matx[i][j];
            }
        }
        
        return Matx;
    }
   
    
    
}


#endif /* basic_funcs_hpp */










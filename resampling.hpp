//
//  resampling.hpp
//  
//
//  Created by YangTong on 11/7/17.
//
//
#pragma once

#ifndef resampling_hpp
#define resampling_hpp

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>
#include <stdio.h>



namespace STAT_TEST {
    
    
    // FISHER-YATES-SHUFFLE process for the training data
    template<typename _num_type>
    inline std::vector<_num_type>
    Fisher_Yates_shuffle(std::vector<_num_type> _data_total, int _shuffle_times);
    
    // alternatively, one could use basic shuffling process defined in std namespace
    template<typename _num_type>
    inline std::vector<_num_type>
    standard_shuffle(std::vector<_num_type> _data_total, int _shuffle_times);
    
    
    
    
    
    
    // RESAMPLING process for a given group of data
    // --return   -> {part1, part2}     ::  one shuffled data divided into two groups: usually smaller ones are samples
    // --parameter-> _data_total        ::  training data
    // --parameter-> _subsample_size    ::  size of the pick-out sample
    // --parameter-> _resample_interval ::  number of shuffling in each resampling before pick out samples
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total,
               int _subsample_size,
               int _resample_interval
               );
   
    // RESAMPLING process with
    // default resampling_interval = 1
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total,
               int _subsample_size
               );
    
    
    // RESAMPLING process with
    // default resampling_interval = 1
    // default pick-out sample size = half of total size
    template<typename _num_type>
    inline std::tuple< std::vector<_num_type>, std::vector<_num_type> >
    resampling(std::vector<_num_type> _data_total);
    
    
}

#include "resampling.ipp"

#endif /* resampling_hpp */

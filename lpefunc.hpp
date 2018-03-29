//
//  lpefunc.hpp
//  
//
//  Created by YangTong on 2/18/18.
//
//
#pragma once

#ifndef lpefunc_hpp
#define lpefunc_hpp


#include "basic_funcs.hpp"

namespace LOCAL_EST{
    using namespace boost::numeric::ublas;
    using namespace STAT_TEST;
    
    inline Vector lpefunc (Vector   _vecX,
                           Vector   _vecY,
                           Vector   _weight,
                           double   _center,
                           double   _bandwidth,
                           int      _order
                           );
    
    inline Vector lpefunc (std::vector<double>   _vecX,
                           std::vector<double>   _vecY,
                           std::vector<double>   _weight,
                           double   _center,
                           double   _bandwidth,
                           int      _order
                           );
    
}


#include "lpefunc.ipp"
#endif /* lpefunc_hpp */

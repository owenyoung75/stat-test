//
//  TS_models.hpp
//  
//
//  Created by YangTong on 2/19/18.
//
//

#ifndef TS_models_hpp
#define TS_models_hpp


#include "basic_funcs.hpp"

namespace STAT_TEST
{
    using namespace boost::numeric::ublas;
    
    
    struct ARCH
    {
        double  Omega;
        Vector  Alpha;
        double  Var;
        int     Lag;
        int     Time;
        
        Vector  X_t;
        Vector  M_t;
        Vector  Sigma_t;
        Vector  U_t;
        Vector  Series;
        
        ARCH(double _omega, Vector _alpha, int _time);
        
        void generator(double (*Func)(double));
        void generator();
  
    };
    
    
    
    
    struct GARCH
    {
        double  Omega;
        Vector  Alpha;
        Vector  Beta;
        double  Var;
        int     Lag1;
        int     Lag2;
        int     Time;
        
        Vector  X_t;
        Vector  M_t;
        Vector  Sigma_t;
        Vector  U_t;
        Vector  Series;
        
        GARCH(double _omega, Vector _alpha, Vector _beta, int _time);
        
        void generator(double (*Func)(double));
        void generator();
    };
    
    
    
}




#include "TS_models.ipp"
#endif /* TS_models_hpp */

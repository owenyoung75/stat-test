//
//  estimator_comparison.cpp
//  
//
//  Created by YangTong on 2/21/18.
//
//

#include <stdio.h>


#include "basic_funcs.hpp"
#include "TS_models.hpp"

#include "local_estimator.hpp"

#define UTC (-5)

using namespace boost::numeric::ublas;
using namespace LOCAL_EST;


double Func_tobe_Tested(double x0){
    return x0 * x0;
}

Vector Func_tobe_Tested(Vector x0){
    for (int i = 0; i < x0.size(); i++)     x0(i)  =  Func_tobe_Tested(x0(i));
    return x0;
}




int main()
{
    // GARCH parameters
    int     Time    =   1000;
    double  Omega   =   1.0;
    Vector  Alpha(1);       Alpha(0) = 0.4;
    Vector  Beta(1);        Beta(0)  = 0.3;
    
    GARCH    Garch(Omega, Alpha, Beta, Time);
    
    // Estimation parameters
    int     truncation  =   int(3 * std::log(Time));
    int     order       =   3;
    double  winsor      =   0.01;
    double  d0          =   5.5;    double  h0;
    double  d1          =   3.5;    double  h1;
    
    
    // Testing parameters
    int     repetition  =   100;
    double  x0          =   0.0;
    double  m0          =   Func_tobe_Tested(x0);
    
    
    
    // Estimation comparison
    Matrix  m_value(4, repetition);
    Matrix  m_bias(4, repetition);
    Vector  est_rlt;
    
    for (int i = 0; i < repetition; i++)
    {
        time_t rawtime0;
        time ( &rawtime0 );
        
        Garch.generator(Func_tobe_Tested);
        double  Var_Xt = Variance_value(Garch.X_t);
        h0 = d0 * Var_Xt * std::pow(Time*1.0, -1.0/6.0);
        h1 = d1 * Var_Xt * std::pow(Time*1.0, -1.0/9.0);
        
        est_rlt = InfExact_LOWESS(Garch, x0, m0, h0, order);
        m_value(0, i) = est_rlt(0);    m_bias(0, i) = est_rlt(1);
        
        est_rlt = LOESS(Garch, x0, m0, h0, order);
        m_value(1, i) = est_rlt(0);    m_bias(1, i) = est_rlt(1);
        
        est_rlt = LOWESS_G11(Garch, x0, m0, h0, h1, winsor, order, truncation);
        m_value(2, i) = est_rlt(0);    m_bias(2, i) = est_rlt(1);
        
        est_rlt = InfApprx_LOWESS_G11(Garch, x0, m0, h0, winsor, order, truncation);
        m_value(3, i) = est_rlt(0);    m_bias(3, i) = est_rlt(1);
        
        
        time_t rawtime1;
        time( &rawtime1 );
        double t = difftime(rawtime1, rawtime0);
        std::cout<< "    This estimation took " << t << " seconds." << '\n' << '\n';
    }
    
    Vector empBias(4);
    Vector empVar(4);
    Vector empMES(4);

    for (int i  = 0; i < 4; i++)
    {
        Vector bias = row(m_bias, i);
        empBias(i) = Mean_value(bias);
        empVar(i)  = Variance_value(bias);
        empMES(i)  = empBias(i) * empBias(i) + empVar(i);
    }
    
    
    std::cout <<  empBias << std::endl;
    std::cout <<  empMES << std::endl;
    
    
    
 
    return 0;
}

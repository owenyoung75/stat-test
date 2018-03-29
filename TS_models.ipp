//
//  TS_models.ipp
//  
//
//  Created by YangTong on 2/19/18.
//
//

#pragma once


namespace STAT_TEST
{
    ARCH::ARCH(double _omega,
               Vector _alpha,
               int _time
               ):
    Omega(_omega),
    Alpha(_alpha),
    Time(_time),
    Lag(Alpha.size()),
    X_t(Time, 0.0),
    M_t(Time, 0.0),
    Sigma_t(Time, 0.0),
    U_t(Time, 0.0),
    Series(Time, 0.0)

    {
        std::random_device rd{};
        std::mt19937_64 gen{rd()};
        std::normal_distribution<>          normal{0.0,1.0};
        std::uniform_real_distribution<>    uniform(0.0,1.0);
        
        Var     =   Omega/(1.0 - accumulate(Alpha.begin(), Alpha.end(), 0.0));
        if (Time <= Lag) Time = Lag * 2;
        
        
        Vector u_lag (Lag,  Var);
        
        for (int i = 0; i < Time*2; i++)
        {
            Vector  ulag_sq =   element_prod(u_lag, u_lag);
            double  sigma   =   sqrt( Omega + inner_prod(Alpha, ulag_sq) );
            
            if (Lag > 1)
                subrange(u_lag, 1, u_lag.size())  =   subrange(u_lag, 0, u_lag.size()-1);
            u_lag(0) = sigma * normal(gen);
        }
        
        for (int i = 0; i < Time; i++)
        {
            Vector  ulag_sq = element_prod(u_lag, u_lag);
            double sigma  =   sqrt( Omega + inner_prod(Alpha, ulag_sq) );
            
            if (Lag > 1)
                subrange(u_lag, 1, u_lag.size())  =   subrange(u_lag, 0, u_lag.size()-1);
            u_lag(0) = sigma * normal(gen);
            
            Sigma_t(i) = sigma;
            U_t(i) = u_lag(0);
            X_t(i) = uniform(gen);
        }
    }

    
    void ARCH::generator(double (*Func)(double))
    {
        Vector  _m_t(Time);
        for (int i = 0; i < Time; i++) _m_t(i) = (*Func)(X_t(i));
        M_t = _m_t;
        Series = M_t + U_t;
    }

    
    void ARCH::generator()
    {
        M_t = X_t;
        Series = M_t + U_t;
    }
    
    
    
    
    GARCH::GARCH(double _omega,
                 Vector _alpha,
                 Vector _beta,
                 int _time
                 ):
    Omega(_omega),
    Alpha(_alpha),
    Beta(_beta),
    Time(_time),
    Lag1(Alpha.size()),
    Lag2(Beta.size()),
    X_t(Time, 0.0),
    M_t(Time, 0.0),
    Sigma_t(Time, 0.0),
    U_t(Time, 0.0),
    Series(Time, 0.0)
    {
        std::random_device rd{};
        std::mt19937_64 gen{rd()};
        std::normal_distribution<>          normal{0.0,1.0};
        std::uniform_real_distribution<>    uniform(0.0,1.0);
        
        Var =   Omega/(1.0 - accumulate(Alpha.begin(), Alpha.end(), 0.0)- accumulate(Beta.begin(), Beta.end(), 0.0));
        if (Time < Lag1 || Time < Lag2) Time = Lag1 + Lag2;

        
      
        Vector u_lag (Lag1, Var);
        Vector s_lag (Lag2, Var);
        
        for (int i = 0; i < Time*2; i++)
        {
            Vector  ulag_sq = element_prod(u_lag, u_lag);
            Vector  slag_sq = element_prod(s_lag, s_lag);
            
            double sigma  =   sqrt( Omega + inner_prod(Alpha, ulag_sq) + inner_prod(Beta, slag_sq) );
            
            if (Lag1 > 1)
                subrange(u_lag, 1, Lag1)  =   subrange(u_lag, 0, Lag1-1);
            if (Lag1 > 2)
                subrange(s_lag, 1, Lag2)  =   subrange(s_lag, 0, Lag2-1);
            
            s_lag(0) = sigma;
            u_lag(0) = sigma * normal(gen);
        }
        
        for (int i = 0; i < Time; i++)
        {
            Vector  ulag_sq = element_prod(u_lag, u_lag);
            Vector  slag_sq = element_prod(s_lag, s_lag);
            
            double sigma  =   sqrt( Omega + inner_prod(Alpha, ulag_sq) + inner_prod(Beta, slag_sq) );
            
            if (Lag1 > 1)
                subrange(u_lag, 1, Lag1)  =   subrange(u_lag, 0, Lag1-1);
            if (Lag1 > 2)
                subrange(s_lag, 1, Lag2)  =   subrange(s_lag, 0, Lag2-1);
           
            s_lag(0) = sigma;
            u_lag(0) = sigma * normal(gen);
            
            Sigma_t(i) = sigma;
            U_t(i) = u_lag(0);
            X_t(i) = uniform(gen);
        }
    }
    
    
    void GARCH::generator(double (*Func)(double))
    {
        Vector  _m_t(Time);
        for (int i = 0; i < Time; i++) _m_t(i) = (*Func)(X_t(i));
        M_t = _m_t;
        Series = M_t + U_t;
    }
    
    
    void GARCH::generator()
    {
        M_t = X_t;
        Series = M_t + U_t;
    }
    
    
    
}

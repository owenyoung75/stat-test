//
//  local_estimator.ipp
//  
//
//  Created by YangTong on 2/20/18.
//
//

namespace LOCAL_EST
{
    void WinsorG11(double *omega, double *alpha, double *beta, double _winsor)
    {
        if (*omega < 0) *omega = 0;
        
        if (*alpha < 0) *alpha = _winsor;
        else if (*alpha > 1)    *alpha = 1 - _winsor;
        
        if (*beta < 0) *beta = _winsor;
        else if (*beta > 1)    *beta = 1 - _winsor;
    }
    
    template<class T>
    inline Vector GH11_Two_Step_Estimator(T      _ts_model,
                                          Vector _Usq_t,
                                          double _x0,
                                          double _m0,
                                          double _winsor,
                                          double _bandwidth,
                                          int    _order,
                                          int    _trunc
                                          )
    {
        Vector rlt(2);
        
        int    time = _ts_model.Time;
        double var  = _ts_model.Var;
        Vector vecX = _ts_model.X_t;
        Vector vecY = _ts_model.Series;
        
        Vector ones(time, 1.0);
        Vector SigmahatSq_t(time);
        double omegahat;
        double alphahat;
        double betahat;
        
        // First step
        Matrix U_lag(_trunc+1, time - _trunc);
        for (int i=0; i<(_trunc+1); i++)
            row(U_lag, i) = subrange(_Usq_t, i, i + time - _trunc );
        Matrix CorrMat = Covariance_Matrix<double>(U_lag);
        Vector rho = abs(subrange(row(CorrMat, 0), 1, _trunc+1));
        Vector phis = element_div(subrange(rho, 1, _trunc), subrange(rho, 0, _trunc-1));
        double phi = Mean_value( phis );
        if (phi<0)
            phi = 0;
        else if(phi > 1)
            phi = 1.0;
        
        double b = (phi * phi + 1 - 2 * rho(0) * phi)/(phi - rho(0));   if (b < (2.0+_winsor))  b = 2.0 + _winsor;
        double theta = -std::abs((-b+sqrt(b*b-4.0))/2.0);
        alphahat = theta + phi;
        betahat = -theta;
        omegahat = Mean_value(_Usq_t)*(1-phi);
        WinsorG11(&omegahat, &alphahat, &betahat, _winsor);
        
        SigmahatSq_t(0) = omegahat + betahat*var + alphahat*var;
        for (int i = 1; i < time; i++)
            SigmahatSq_t(i) = omegahat + alphahat*_Usq_t(i-1) + betahat*SigmahatSq_t(i-1);

        // Second step
        Vector UDep_t = subrange(element_div(_Usq_t, SigmahatSq_t) , 1, time);
        Matrix UIndep(3, time -1);
        row(UIndep, 0) = element_div( subrange(ones,         1, time),   subrange(SigmahatSq_t, 1, time) );
        row(UIndep, 1) = element_div( subrange(_Usq_t,       0, time-1), subrange(SigmahatSq_t, 1, time) );
        row(UIndep, 2) = element_div( subrange(SigmahatSq_t, 0, time-1), subrange(SigmahatSq_t, 1, time) );
        
        UDep_t = prod(UIndep, UDep_t);
        UIndep = prod(UIndep, trans(UIndep));
        Matrix inverse(3, 3);   InvertMatrix<double>(UIndep, inverse);
        UDep_t = prod( inverse, UDep_t);
        alphahat = UDep_t(1);
        betahat = UDep_t(2);
        omegahat = UDep_t(0);
        WinsorG11(&omegahat, &alphahat, &betahat, _winsor);
        
        SigmahatSq_t(0) = omegahat + betahat*var + alphahat*var;
        for (int i = 1; i < time; i++)
            SigmahatSq_t(i) = omegahat + alphahat*_Usq_t(i-1) + betahat*SigmahatSq_t(i-1);
        
        // Have weights now, start lpefunc
        rlt(0) = lpefunc( vecX, vecY, SigmahatSq_t, _x0, _bandwidth, _order )(0);
        rlt(1) = rlt(0) - _m0;
        return rlt;
    }
    
    
    
    
    
    template<class T>
    inline Vector LOESS(T      _ts_model,
                        double _x0,
                        double _m0,
                        double _bandwidth,
                        int    _order)
    {
        Vector rlt(2);
        
        Vector vecX = _ts_model.X_t;
        Vector vecY = _ts_model.Series;
        Vector weight(vecX.size(), 1.0);
        
        rlt(0) = lpefunc( vecX, vecY, weight, _x0, _bandwidth, _order )(0);
        rlt(1) = rlt(0) - _m0;
        return rlt;
    }
    
    
    template<class T>
    inline Vector InfExact_LOWESS(T      _ts_model,
                                  double _x0,
                                  double _m0,
                                  double _bandwidth,
                                  int    _order)
    {
        Vector rlt(2);
        
        Vector vecX = _ts_model.X_t;
        Vector vecY = _ts_model.Series;
        Vector weight = element_prod(_ts_model.Sigma_t, _ts_model.Sigma_t);
        
        rlt(0) = lpefunc( vecX, vecY, weight, _x0, _bandwidth, _order )(0);
        rlt(1) = rlt(0) - _m0;
        return rlt;
    }
    
    
    
    template<class T>
    inline Vector LOWESS_G11(T      _ts_model,
                             double _x0,
                             double _m0,
                             double _bandwidth0,
                             double _bandwidth1,
                             double _winsor,
                             int    _order,
                             int    _trunc)
    {
        int    time = _ts_model.Time;
        Vector vecX = _ts_model.X_t;
        Vector vecY = _ts_model.Series;
        Vector weight(time, 1.0);
        Vector Uhat_t(time, 0.0);
        
        for (int i=0; i<time; i++){
            Uhat_t(i) = lpefunc( vecX, vecY, weight, vecX(i), _bandwidth1, _order )(0);
        }
        
        Uhat_t = element_prod(Uhat_t, Uhat_t);
        
        return GH11_Two_Step_Estimator(_ts_model, Uhat_t, _x0, _m0, _winsor, _bandwidth0, _order, _trunc);
    }
    
    
    template<class T>
    inline Vector InfApprx_LOWESS_G11(T      _ts_model,
                                      double _x0,
                                      double _m0,
                                      double _bandwidth,
                                      double _winsor,
                                      int    _order,
                                      int    _trunc
                                      )
    {
        Vector Usq_t(_ts_model.Time);
        
        Usq_t = element_prod(_ts_model.U_t, _ts_model.U_t);
        
        return GH11_Two_Step_Estimator(_ts_model, Usq_t, _x0, _m0, _winsor, _bandwidth, _order, _trunc);
    }
    
}

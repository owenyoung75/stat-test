//
//  local_estimator.hpp
//  
//
//  Created by YangTong on 2/20/18.
//
//

#ifndef local_estimator_hpp
#define local_estimator_hpp

#include <stdio.h>

#include "lpefunc.hpp"
#include "basic_funcs.hpp"

namespace LOCAL_EST {
    
    void WinsorG11(double *omega, double *alpha, double *beta, double _winsor);
 
    template<class T>
    inline Vector GH11_Two_Step_Estimator(T      _ts_model,
                                          Vector _Usq_t,
                                          double _x0,
                                          double _m0,
                                          double _winsor,
                                          double _bandwidth,
                                          int    _order,
                                          int    _trunc
                                          );

    template<class T>
    inline Vector LOESS(T      _ts_model,
                        double _x0,
                        double _m0,
                        double _bandwidth,
                        int    _order
                        );
    
    
    template<class T>
    inline Vector InfExact_LOWESS(T      _ts_model,
                                  double _x0,
                                  double _m0,
                                  double _bandwidth,
                                  int    _order
                                  );
    
    template<class T>
    inline Vector LOWESS_G11(T      _ts_model,
                             double _x0,
                             double _m0,
                             double _bandwidth,
                             double _winsor,
                             int    _order,
                             int    _trunc
                             );

    
    
    template<class T>
    inline Vector InfApprx_LOWESS_G11(T      _ts_model,
                                      double _x0,
                                      double _m0,
                                      double _bandwidth,
                                      double _winsor,
                                      int    _order,
                                      int    _trunc
                                      );

}


#include "local_estimator.ipp"
#endif /* loca_estimator_hpp */

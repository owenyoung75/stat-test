//
//  lpefunc.ipp
//  
//
//  Created by YangTong on 2/18/18.
//
//

#pragma once


namespace LOCAL_EST
{
    inline Vector lpefunc (Vector   _vecX,
                           Vector   _vecY,
                           Vector   _weight,
                           double   _center,
                           double   _bandwidth,
                           int      _order
                           )
    {
        int data_size = _vecX.size();
        
        Matrix  dX(data_size, _order+1, 1.0);
        unbounded_array<double>  diag (data_size);
        Vector  BetaVec (_order+1);
        
        for (int i = 0; i < data_size; i++){
            _vecX[i] = (_vecX[i] - _center)/_bandwidth;
            diag[i] = NormalPDF(_vecX[i])/_weight[i];
        }
        diagonal_matrix<double>  kernel_weighted(data_size, diag);
        
        for (int i = 1; i <= _order; i++){
            column(dX, i) = element_prod(column(dX, i-1), _vecX);
        }
        
        Matrix dXW = prod(trans(dX), kernel_weighted);
        Matrix inverse(_order+1, _order+1);
        BetaVec = prod(dXW, _vecY);

        bool invs = InvertMatrix<double>(prod(dXW, dX), inverse);
        BetaVec = prod(inverse, BetaVec);
        
        return BetaVec;
    }
    
    
    inline Vector lpefunc (std::vector<double>   _vecX,
                           std::vector<double>   _vecY,
                           std::vector<double>   _weight,
                           double   _center,
                           double   _bandwidth,
                           int      _order
                           )
    {
        vector<double, std::vector<double>> vecX(_vecX);
        vector<double, std::vector<double>> vecY(_vecX);
        vector<double, std::vector<double>> weight(_vecX);
        
        return  lpefunc (vecX, vecY, weight, _center, _bandwidth, _order);
    }
    
    
    
}

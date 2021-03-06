//
//  basic_funcs.tpp
//
//
//  Created by YangTong on 1/15/18.
//
//
#pragma once

namespace STAT_TEST
{
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
    
    
    // calculate VARIANCE of a std::vector  --> sigma^2
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
    // calculate VARIANCE of a boost::vector --> sigma^2
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
    // return:: matrix(row_idx, row_idx)
    // row_idx :: # of variables  1, 2, 3 ......
    // col_idx :: # of copies     1, 2, 3 ...... --> have closer address
    template<typename _num_type>
    inline Matrix Covariance_Matrix (matrix<_num_type> _data_groups)
    {
        int num_of_varbs = _data_groups.size1();
        Vector means(num_of_varbs);
        Matrix Matx(num_of_varbs, num_of_varbs);
        
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
    
    
    
    // Calculate inverse matrix by LU decomp
    template<class T> bool InvertMatrix (const matrix<T>& input,
                                         matrix<T>& inverse
                                         )
    {
        typedef permutation_matrix<std::size_t> pmatrix;
        matrix<T> A(input);
        pmatrix pm(A.size1());
        
        int res = lu_factorize(A,pm);
        if( res != 0 ) return false;
        
        inverse.assign(identity_matrix<T>(A.size1()));
        lu_substitute(A, pm, inverse);
        return true;
    }
    
    
    // PDF function of Normal distribution
    inline double NormalPDF (double x, double mu, double var){
        return exp(- (x-mu) * (x-mu) /2.0)/( sqrt(2 * STAT_TEST::PI) * var );
    }
    inline double NormalPDF (double x){
        return exp(- x * x /2.0)/sqrt(2 * STAT_TEST::PI);
    }
    inline Vector NormalPDF (Vector X, double mu, double var){
        Vector Y(X.size());
        for (int i = 0; i < X.size(); i++)
            Y(i) = exp(- (X(i)-mu) * (X(i)-mu) /2.0)/( sqrt(2 * STAT_TEST::PI) * var );
        return Y;
    }
    inline Vector NormalPDF (Vector X){
        Vector Y(X.size());
        for (int i = 0; i < X.size(); i++)
            Y(i) = exp(- X(i)* X(i) /2.0)/sqrt(2 * STAT_TEST::PI);
        return Y;
    }
    
    
    // obtain the coarse-grained histgram
    template<typename _num_type>
    inline std::map<int, int> coarse_grained_distr (std::vector<_num_type> _data)
    {
        std::map<int, int> hist{};
        for(int n=0; n<_data.size(); ++n)
            ++hist[std::round(_data[n])];
        for(auto p : hist)
            hist[p.first] = p.second/int(_data.size()/200);
        return hist;
    }
    
    
    
    // Plot the distribution
    template<typename _num_type>
    inline void Plot_distribution (std::vector<_num_type> _data)
    {
        std::map<int, int> hist = coarse_grained_distr<_num_type>(_data);
        
        std::cout << '\n' << '\n' << '\n' << '\n';
        for(auto p : hist)
            std::cout << std::setw(2) << p.first << ' ' << std::string(p.second, '*') << '\n';
        std::cout << '\n' << '\n' << '\n' << '\n';
    }
    
    
}











#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <string>
#include <iomanip>
#include <map>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>


#include "basic_funcs.hpp"
#include "empirical_distr.hpp"
#include "resampling.hpp"
#include "equal_test_statistics.hpp"
#include "equal_distr_test.hpp"


typedef std::vector<double> Vector;

using namespace STAT_TEST::EQUAL_TEST;
int main()
{
    using namespace boost::numeric::ublas;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> Normal_D(0,0.3);
    Vector _total_sample;
    for (int n = 0; n < 1000; n++)
        _total_sample.push_back(Normal_D(gen));
    
    int _size_of_subsample = 50;
    int _num_of_resampling = 1000;
    int _highest_order_poly = 3;
    int _resampling_interval = 2;
    int _wamup_steps = 20;
    

    for (auto i = 0; i < 50; i++)
    {
        
    
    matrix<double> PI = PI_package<double> (_total_sample,
                                            _size_of_subsample,
                                            _num_of_resampling,
                                            _highest_order_poly,
                                            _resampling_interval,
                                            _wamup_steps);

    matrix<double> Cov = STAT_TEST::Covariance_Matrix(PI);
    
    std::cout << Cov << std::endl;
//    for (int r = 0; r<=Cov.size1(); r++)
//        std::cout << row(Cov, r) << std::endl;
    
    }
    return 0;
}

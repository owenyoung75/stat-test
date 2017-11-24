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
#include <boost/range/algorithm.hpp>

#include "basic_funcs.hpp"
#include "equal_distr_test.hpp"

using namespace STAT_TEST::EQUAL_TEST;
typedef std::vector<double> Vector;

int main()
{
    using namespace boost::numeric::ublas;
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        Construct random machines       /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::normal_distribution<> training_D(0.0,1.0);
//    std::normal_distribution<>        testing_D(0.0,1.0);
//    std::student_t_distribution<>     testing_D(2);
//    std::chi_squared_distribution<>   testing_D(2);
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           Training parameters           ////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    int     training_sample_size = 1000;
    int     subsample_size = 400;
    int     num_of_resampling = 10000;
    int     highest_order_poly = 4;
    int     resampling_interval = 1;
    int     wamup_steps = 10;
    double  critical_portion = STAT_TEST::Nv1;
    bool    plot_or_not = true;
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////           Testing parameters           /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    int     testing_sample_size = int(training_sample_size * subsample_size / (training_sample_size - subsample_size));
    int     testing_times = 1000;
    

    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////// Construct empirical PSI for given sample ////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    
    Vector training_sample;
    for (int n = 0; n < training_sample_size; n++)
        training_sample.push_back(training_D(gen));


    std::tuple< Vector, Vector >  pair = Pi_and_Psi_ciritical_values<double> (training_sample,
                                                                              subsample_size,
                                                                              num_of_resampling,
                                                                              highest_order_poly,
                                                                              resampling_interval,
                                                                              wamup_steps,
                                                                              critical_portion,
                                                                              plot_or_not
                                                                              );

    Vector lower_bound = std::get<0>(pair);
    Vector upper_bound = std::get<1>(pair);
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////// Obtain testing PSI for testing data groups ///////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    
    STAT_TEST::Empirical_Distribution<double> ED_training(training_sample);

    std::vector<int> DOF {1, 3, 5, 10, 15, 20, 30, 50};
    for (int d= 0; d < DOF.size(); d++)
    {
    std::ofstream ofile("tested_student_t_" + std::to_string(DOF[d]) + ".txt",std::ios::out);
    std::student_t_distribution<> testing_D(DOF[d]);
    
    std::map<int, int> hist;
    int violations = 0;
    for (int k = 0; k < testing_times; k++)
    {
        Vector testing_sample = {};
        for (int n = 0; n < testing_sample_size; n++)
            testing_sample.push_back(testing_D(gen));
        
        Vector pi_test = Test_vector_PI<double> (ED_training,
                                                 testing_sample,
                                                 highest_order_poly);
        pi_test.push_back(PSI_square(pi_test));
        ++hist[std::round(pi_test[highest_order_poly])];

        if (pi_test[highest_order_poly] > upper_bound[highest_order_poly]
            || pi_test[highest_order_poly] < lower_bound[highest_order_poly])
            violations += 1;
    }
    
    std::cout<< DOF[d] << int(testing_times * (1-critical_portion) )<< "    " << violations << std::endl;
    ofile<<'\n'<<'\n'<<'\n'<< testing_times<< "    " << int(testing_times * (1-critical_portion) )<< "    " <<violations;
    for(auto p : hist)
    {
        ofile << std::fixed << std::setprecision(1) << std::setw(2)
        << p.first << ' ' << std::string(p.second/int(testing_times/200) , '*') << '\n';
    }

    }
    
    return 0;
}






//    // Testing under N>>M condition
//    // which converges to a chi_square_distribution
//    std::map<int, int> hist;
//    for(int n=0; n<10000; ++n)
//    {
//        Vector cdf = {};
//        for (int n = 0; n< 1000; n++)
//            cdf.push_back(Uniform(gen));
//
//        double Pi = 0;
//        for (int i = 1; i<= 6; i++)
//        {
//            Vector pi_i = Normalized_Legendre_Poly (i, cdf);
//
//            double pi2 = std::accumulate(pi_i.begin(), pi_i.end(), 0.0) / sqrt(cdf.size());
//            Pi += pi2 * pi2;
//        }
//
//        ++hist[std::round(Pi)];
//    }
//    for(auto p : hist) {
//        std::cout << std::fixed << std::setprecision(1) << std::setw(2)
//        << p.first << ' ' << std::string(p.second/200, '*') << '\n';
//    }






//    // test mutiple times and count H0-violating-events
//    int violations = 0;
//    std::map<int, int> hist;
//    for (int k = 0; k<testing_times; k++)
//    {
//        Vector testing_sample = {};
//        for (int n = 0; n < testing_sample_size; n++)
//            testing_sample.push_back(testing_D(gen));
//
//        Vector pi_test = Pi_and_Psi_tested(training_sample,
//                                           testing_sample,
//                                           testing_times
//                                           );
//        ++hist[std::round(pi_test[highest_order_poly])];
//
//        vector<bool> test_rlt(highest_order_poly + 1);
//        for (int i = 0; i < highest_order_poly; i++)
//        {
//            if (pi_test[i] < upper_bound[i] && pi_test[i] > lower_bound[i])
//                test_rlt(i) = true;
//            else
//                test_rlt(i) = false;
//        }
//
//        if (pi_test[highest_order_poly] < upper_bound[highest_order_poly]
//            && pi_test[highest_order_poly] > lower_bound[highest_order_poly])
//        {test_rlt(highest_order_poly) = true;}
//        else
//        {test_rlt(highest_order_poly) = false;  violations += 1;}
//
//        std::cout<< test_rlt << std::endl;
//    }
//    std::cout<< testing_times<< "    " << violations << std::endl;

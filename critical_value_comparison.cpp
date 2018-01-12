/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////                                                          //////////////////
/////////////////    Main program 1 for Equal_distribution_test problem    //////////////////
/////////////////                                                          //////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////      critical_value_comparison.cpp      ////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////    Created by YangTong on 11/10/17.     ////////////////////////////
////////////////////////                                         ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
#include <random>
#include <string>
#include <iomanip>
#include <map>
#include <cmath>
#include <chrono>
#include <thread>
#include <unistd.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/algorithm.hpp>

#include "basic_funcs.hpp"
#include "equal_distr_test.hpp"

using namespace STAT_TEST::EQUAL_TEST;
typedef std::vector<double> Vector;
#define UTC (+19)

int main()
{
    std::cout << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "    Computation started ......" << '\n'<< '\n'<< '\n'<< '\n';

    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        Construct random machines       /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::random_device rd{};
    std::mt19937_64 gen0{rd()};
    
    std::normal_distribution<> training_D{0.0,1.0};
//    std::normal_distribution<> testing_D {0.0,5.0};
    std::student_t_distribution<>     testing_D{15};
//    std::chi_squared_distribution<>   testing_D{5};
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    //////////////////           Training & Testing parameters           ///////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Parameters Setting:" << '\n';
    
    int     training_sample_size = 5000;
    int     testing_sample_size  = 0.1 * training_sample_size;
    int     subsample_size = int(testing_sample_size * training_sample_size / (testing_sample_size + training_sample_size));
    int     testing_times = 5000;
    
    int     num_of_resampling = 10000;
    int     highest_order_poly = 4;
    int     resampling_interval = 1;
    int     wamup_steps = 10;
    double  critical_portion = STAT_TEST::Nv5;
    bool    plot_or_not = true;

    
    std::cout << "    training sample size     = " << training_sample_size << '\n';
    std::cout << "    testing  sample size     = " << testing_sample_size  << '\n';
    std::cout << "    subsampling size         = " << subsample_size << '\n';
    std::cout << "    resampling times         = " << num_of_resampling << '\n';
    std::cout << "    testing times            = " << testing_times << '\n';
    std::cout << "    critical value           = " << (1-critical_portion)*100 << "%" << '\n';
    std::cout << "    highest polynomial-order = " << highest_order_poly << '\n';
    std::cout << "    ........" << '\n';
    std::cout << "    Preference Setting:" << '\n';
    std::cout << "    shuffle " << wamup_steps << " times for warm-up before resampling." << '\n';
    std::cout << "    shuffle " << resampling_interval << " times for each resampling." << '\n';
    std::cout << "    training sample forms a normal distribution as N(0,5)." << '\n';
    if (plot_or_not)
        std::cout << "    empirical distribution of the training sample would be plotted." << '\n';
    else
        std::cout << "    empirical distribution of the training sample would NOT be plotted." << '\n';
    std::cout << '\n' << '\n' << '\n';
    
    
    
    

    for (int it = 1; it <= 5; it++)
    {
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////    Record the current time as file names    //////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    time_t rawtime0;
    struct tm * ptm;
    time ( &rawtime0 );
    ptm = gmtime ( &rawtime0 );
    std::string c_time = std::to_string((ptm->tm_hour)%24+UTC) + ":" + std::to_string(ptm->tm_min)+ ":" + std::to_string(ptm->tm_sec);
    std::cout << "    Current time: " << c_time << '\n';
    
    std::string training_file = c_time + "_trained_" + std::to_string(int(testing_sample_size*10/training_sample_size)) + ":10" + ".txt";
    std::string tesing_file   = c_time + "_student15_"  + std::to_string(int(testing_sample_size*10/training_sample_size)) + ":10" + ".txt";
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////      Construct empirical PSI for given sample      ///////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    Vector training_sample;
        
    for (int n = 0; n < training_sample_size; n++)
        training_sample.push_back(training_D(gen0));
    std::map<int, int> hist1 = STAT_TEST::coars_grained_distr<double>(training_sample);
    

    std::tuple< Vector, Vector >  pair = Pi_and_Psi_ciritical_values<double> (training_sample,
                                                                              subsample_size,
                                                                              num_of_resampling,
                                                                              highest_order_poly,
                                                                              resampling_interval,
                                                                              wamup_steps,
                                                                              critical_portion,
                                                                              plot_or_not,
                                                                              training_file
                                                                              );

    Vector lower_bound = std::get<0>(pair);
    Vector upper_bound = std::get<1>(pair);
    
    
    time_t rawtime1;
    time( &rawtime1 );
    double t = difftime(rawtime1, rawtime0);
    std::cout<< "    Training part took " << t << " seconds." << '\n' << '\n';
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////      Obtain testing PSI for testing data groups      //////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Testing the chosen distribution for "<< testing_times << " times ......" << std::endl;
    
    std::map<int, int> hist{};
    int violation_l = 0;
    int violation_r = 0;

    int sum = 0;
    for (int k = 0; k < testing_times; k++)
    {
        // obtain a testing sample
        Vector testing_sample = {};
        for (int n = 0; n < testing_sample_size; n++)
            testing_sample.push_back(testing_D(gen0));
        
        // roughly estimate the distribution difference
        std::map<int, int> hist2 = STAT_TEST::coars_grained_distr<double>(testing_sample);
        for (auto p : hist1)
        {
            if (hist2.find(p.first) != hist2.end())
                sum += abs(p.second - hist2[p.first]);
            else
                sum += p.second;
        }
        for (auto p : hist2)
            if (hist1.find(p.first) == hist1.end())    sum += p.second;
    
        // coarse-grained histogram of testing samples
        Vector pi_test = Pi_and_Psi_tested<double>(training_sample, testing_sample, highest_order_poly);
        ++hist[std::round(pi_test[highest_order_poly])];

        if (pi_test[highest_order_poly] > upper_bound[highest_order_poly]) violation_r += 1;
        if (pi_test[highest_order_poly] < lower_bound[highest_order_poly]) violation_l += 1;
    }
    
    std::cout << "    Testing finished." << '\n' << "    Critical violations recorded as below:" << std::endl;
    std::cout << "    " << testing_times << "    " << int(testing_times * (1-critical_portion) + 0.0001)
              << "    " << violation_l + violation_r << "    " << sum/testing_times << std::endl;
   
        
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////            Plot Psi^2 distributions          /////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Plotting the Psi-square distribution of testing samples ......" << std::endl;
    
    std::ofstream ofile(tesing_file, std::ios::out);
    for(auto p : hist)
        ofile << std::fixed << std::setprecision(1) << std::setw(2)
              << p.first << ' ' << std::string(p.second/int(testing_times/200) , '*') << '\n';
    
    ofile << '\n' << '\n' << '\n';
    ofile << testing_times<< "    "
          << int(testing_times * (1-critical_portion) + 0.0001) << "    "
          << violation_l << "   "<< violation_r << "   "<< violation_l + violation_r;

    
    time_t rawtime2;
    time( &rawtime2 );
    t = difftime(rawtime2, rawtime1);
    std::cout<< "    Testing part took " << t << " seconds." << '\n' << '\n' << '\n';
    }
    
    
    
    std::cout << '\n'<< "    Computation completed." << '\n';
    std::cout << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////                        ////////////////////////////////////
/////////////////////////////////   End of main program  ////////////////////////////////////
/////////////////////////////////                        ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////











/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////   Appendix part for supplementary  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////
////
////
////
////////    // Testing under N>>M condition
////////    // which converges to a chi_square_distribution
////////    std::map<int, int> hist;
////////    for(int n=0; n<10000; ++n)
////////    {
////////        Vector cdf = {};
////////        for (int n = 0; n< 1000; n++)
////////            cdf.push_back(Uniform(gen));
////////
////////        double Pi = 0;
////////        for (int i = 1; i<= 6; i++)
////////        {
////////            Vector pi_i = Normalized_Legendre_Poly (i, cdf);
////////
////////            double pi2 = std::accumulate(pi_i.begin(), pi_i.end(), 0.0) / sqrt(cdf.size());
////////            Pi += pi2 * pi2;
////////        }
////////
////////        ++hist[std::round(Pi)];
////////    }
////////    for(auto p : hist) {
////////        std::cout << std::fixed << std::setprecision(1) << std::setw(2)
////////        << p.first << ' ' << std::string(p.second/200, '*') << '\n';
////////    }
////
////
////
////
////
////
////////    // test mutiple times and count H0-violating-events
////////    int violations = 0;
////////    std::map<int, int> hist;
////////    for (int k = 0; k<testing_times; k++)
////////    {
////////        Vector testing_sample = {};
////////        for (int n = 0; n < testing_sample_size; n++)
////////            testing_sample.push_back(testing_D(gen));
////////
////////        Vector pi_test = Pi_and_Psi_tested(training_sample,
////////                                           testing_sample,
////////                                           testing_times
////////                                           );
////////        ++hist[std::round(pi_test[highest_order_poly])];
////////
////////        vector<bool> test_rlt(highest_order_poly + 1);
////////        for (int i = 0; i < highest_order_poly; i++)
////////        {
////////            if (pi_test[i] < upper_bound[i] && pi_test[i] > lower_bound[i])
////////                test_rlt(i) = true;
////////            else
////////                test_rlt(i) = false;
////////        }
////////
////////        if (pi_test[highest_order_poly] < upper_bound[highest_order_poly]
////////            && pi_test[highest_order_poly] > lower_bound[highest_order_poly])
////////        {test_rlt(highest_order_poly) = true;}
////////        else
////////        {test_rlt(highest_order_poly) = false;  violations += 1;}
////////
////////        std::cout<< test_rlt << std::endl;
////////    }
////////    std::cout<< testing_times<< "    " << violations << std::endl;

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////                                                          //////////////////
/////////////////    Main program 2 for Equal_distribution_test problem    //////////////////
/////////////////                                                          //////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////      critical_vector_comparison.cpp     ////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////                                         ////////////////////////////
////////////////////////    Created by YangTong on 01/10/18.     ////////////////////////////
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
using namespace boost::numeric::ublas;
typedef std::vector<double> Vector;


int main()
{
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "    Computation started ......" << '\n'<< '\n'<< '\n'<< '\n';

    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        Construct random machines       /////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::random_device seed{};
    std::mt19937_64 gen0{seed()};
    

    std::normal_distribution<> training_D{0.0,1.0};
    std::normal_distribution<> testing_D {0.0,1.0};
//    std::student_t_distribution<>     testing_D{5};
//    std::chi_squared_distribution<>   testing_D{5};
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    //////////////////           Training & Testing parameters           ///////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Parameters Setting:" << '\n';
    
    int     training_sample_size = 5000;
    int     testing_sample_size  = 0.1 * training_sample_size;
    int     testing_times = 5000;
    
    int     num_of_resampling = 10000;
    int     highest_order_poly = 4;
    int     resampling_interval = 1;
    int     wamup_steps = 10;
    Vector  critical_portions = {STAT_TEST::Nv2, STAT_TEST::Nv4, STAT_TEST::Nv5};
    bool    plot_or_not = true;

    
    
    std::cout << "    training sample size     = " << training_sample_size << '\n';
    std::cout << "    testing  sample size     = " << testing_sample_size  << '\n';
    std::cout << "    subsampling size         = " << subsample_size << '\n';
    std::cout << "    resampling times         = " << num_of_resampling << '\n';
    std::cout << "    testing times            = " << testing_times << '\n';
    std::cout << "    critical values          = ";
    for (int i = 0; i< critical_portions.size(); i++)
        std::cout << (1-critical_portions[i])*100 << "%" << "   ";
    std::cout << '\n';
    std::cout << "    highest polynomial-order = " << highest_order_poly << '\n';
    std::cout << "    ........" << '\n';
    std::cout << "    Preference Setting:" << '\n';
    std::cout << "    shuffle " << wamup_steps << " times for warm-up before resampling." << '\n';
    std::cout << "    shuffle " << resampling_interval << " times for each resampling." << '\n';
    std::cout << "    training sample forms a normal distribution as N(0,1)." << '\n';
    if (plot_or_not)
        std::cout << "    empirical distribution of the training sample would be plotted." << '\n';
    else
        std::cout << "    empirical distribution of the training sample would NOT be plotted." << '\n';
    std::cout << '\n' << '\n';
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////    Record the current time as file names    //////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    time_t rawtime0;
    struct tm * ptm;
    time ( &rawtime0 );
    ptm = localtime ( &rawtime0 );
    std::string c_time = std::to_string((ptm->tm_hour)%24) + ":" + std::to_string(ptm->tm_min)+ ":" + std::to_string(ptm->tm_sec);
    std::cout << "    Current time: " << c_time << '\n';
    
    std::string training_file = c_time + "_trained_" + std::to_string(int(testing_sample_size*10/training_sample_size)) + ":10" + ".txt";
    std::string tesing_file   = c_time + "_student5_"  + std::to_string(int(testing_sample_size*10/training_sample_size)) + ":10" + ".txt";
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////      Construct empirical PSI for given sample      ///////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    Vector training_sample;
    for (int n = 0; n < training_sample_size; n++)
        training_sample.push_back(training_D(gen0));

    int subsample_size = int(testing_sample_size * training_sample_size / (testing_sample_size + training_sample_size));
    
    std::tuple< matrix<double>, matrix<double> >  pair = Pi_and_Psi_ciritical_vectors<double> (training_sample,
                                                                                              subsample_size,
                                                                                              num_of_resampling,
                                                                                              highest_order_poly,
                                                                                              resampling_interval,
                                                                                              wamup_steps,
                                                                                              critical_portions,
                                                                                              plot_or_not,
                                                                                              training_file
                                                                                              );
                                                                              
    
    matrix<double> lower_bounds = std::get<0>(pair);
    matrix<double> upper_bounds = std::get<1>(pair);
    
    
    time_t rawtime1;
    time( &rawtime1 );
    double t = difftime(rawtime1, rawtime0);
    std::cout<< "    Training part took " << t << " seconds." << '\n' << '\n';
    
    
    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////      Obtain testing PSI for testing data groups      //////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Testing the chosen distribution for "<< testing_times << " times ......" << std::endl;
    
    std::map<int, int> hist{};
    std::map<double, int> violation_lower{};
    std::map<double, int> violation_upper{};
    for (int i = 0; i < critical_portions.size(); i++)
    {
        violation_lower[critical_portions[i]] = 0;
        violation_upper[critical_portions[i]] = 0;
    }
    
    for (int k = 0; k < testing_times; k++)
    {
        // obtain a testing sample
        Vector testing_sample = {};
        for (int n = 0; n < testing_sample_size; n++)
            testing_sample.push_back( testing_D(gen0) );
        
        // coarse-grained histogram of testing samples
        Vector pi_test = Pi_and_Psi_tested<double>(training_sample, testing_sample, highest_order_poly);
        ++hist[std::round(pi_test[highest_order_poly])];
        
        for (int i = 0; i < critical_portions.size(); i++)
        {
            if (pi_test[highest_order_poly] > upper_bounds(i, highest_order_poly)) violation_upper[critical_portions[i]] += 1;
            if (pi_test[highest_order_poly] < lower_bounds(i, highest_order_poly)) violation_lower[critical_portions[i]] += 1;
        }
    }
    
    std::cout << "    Testing finished." << '\n' << "    Critical violations recorded as below:" << std::endl;
    for (int i = 0; i < critical_portions.size(); i++)
    {
        std::cout << "    " << "    " << (1-critical_portions[i])*100 << "%";
        std::cout << "    " << testing_times << "    " << int(testing_times * (1-critical_portions[i]) + 0.0001);
        std::cout << "    " << violation_upper[critical_portions[i]] + violation_lower[critical_portions[i]] << std::endl;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////            Plot Psi^2 distributions          /////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << "    Plotting the Psi-square distribution of testing samples ......" << std::endl;
    
    std::ofstream ofile(tesing_file, std::ios::out);
    for(auto p : hist)
        ofile << std::fixed << std::setprecision(1) << std::setw(2)
              << p.first << ' ' << std::string(p.second/int(testing_times/200) , '*') << '\n';
    
    ofile << '\n' << '\n' << '\n';
    for (int i = 0; i < critical_portions.size(); i++)
    {
        ofile << "    " << (1-critical_portions[i])*100 << "%";
        ofile << "    " << testing_times << "    " << int(testing_times * (1-critical_portions[i]) + 0.0001);
        ofile << "    " << violation_upper[critical_portions[i]] + violation_lower[critical_portions[i]] << std::endl;
    }
    ofile.close();
    
    
    time_t rawtime2;
    time( &rawtime2 );
    t = difftime(rawtime2, rawtime1);
    std::cout<< "    Testing part took " << t << " seconds." << '\n' << '\n' << '\n';

    
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    std::cout << '\n'<< "    Computation completed." << '\n';
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    return 0;
}



/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////                        ////////////////////////////////////
/////////////////////////////////   End of main program  ////////////////////////////////////
/////////////////////////////////                        ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


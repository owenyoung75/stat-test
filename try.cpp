//
//  try.cpp
//  
//
//  Created by YangTong on 2/18/18.
//
//
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

#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/assign/list_inserter.hpp>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>



using namespace boost::numeric::ublas;


#include "local_estimator.hpp"
#include "TS_models.hpp"

void chgInt(int *p){
    if (*p > 1)
    *p = 1;
}

int chg(int p){
    return p+1;
}


int myFunction(int (*otherFunction)(int ), int ipt) {
    return (*otherFunction)(ipt);
}

int f(int a)
{
    a += 1;
    return a;
}


int main()
{
    std::vector<double> v = {0.2};
    vector<double, std::vector<double> > v2(std::vector<double> {0.2});

    
    STAT_TEST::ARCH    model(1.0, v2, 10);
    
    
    std::cout << model.X_t << std::endl;
    
    
    return 0;
}

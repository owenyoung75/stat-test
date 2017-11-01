#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>

#include "num_ex.hpp"

typedef std::vector<float> Vector;

using namespace std;

int main()
{

//    static const float arr[] = {16.0, 3.0, 100.0, 120.0, 45.0};
//    Vector Y (arr, arr + sizeof(arr) / sizeof(arr[0]) );
//
////    Vector Z0 = Z<Vector>(Y);
//    printf( "%f\n" , edf_estimator::EDF<float, Vector>(45.0, edf_estimator::X));
    
    Num_ob n;
    cout << n.getNum() << endl;

    return 0;
}

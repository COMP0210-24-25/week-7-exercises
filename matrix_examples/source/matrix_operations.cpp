#include "matrix_maths.h"
#include <random>
#include <vector>
#include <iostream>
#include <chrono>

using std::vector;
typedef std::chrono::high_resolution_clock myClock;

double t_seconds(std::chrono::time_point<myClock> t1, std::chrono::time_point<myClock> t2)
{
    return (t2 - t1).count() / 1e9;
}

int main()
{
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> dist(0,1);

    const int N = 1000;

    // Some options for matrix representations!
    // double *X = new double[N*N];
    // vector<double> X(N*N)
    // double **X = new double *X[N];
    // vector<vector<double>> X(N);   

    return 0;
}
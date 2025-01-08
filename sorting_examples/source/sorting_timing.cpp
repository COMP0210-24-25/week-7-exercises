#include <chrono>
#include "sorting.h"
#include <iostream>
#include <random>

typedef std::chrono::high_resolution_clock myClock;

double t_seconds(std::chrono::time_point<myClock> t1, std::chrono::time_point<myClock> t2)
{
    return (t2 - t1).count() / 1e9;
}

int main()
{

    return 0;
}
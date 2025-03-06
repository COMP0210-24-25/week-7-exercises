#include "Summing.h"
#include <iostream>
#include <vector> 
#include <chrono>
#include <numeric>

using std::vector;

typedef std::chrono::high_resolution_clock myClock;

double t_seconds(std::chrono::time_point<myClock> t1, std::chrono::time_point<myClock> t2)
{
    return (t2 - t1).count() / 1e9;
}

int main()
{
    long N = 160'000'000l;
    float val = 1e-6;
    float *xs = new (std::align_val_t(32)) float[N];
    vector<float> v(N);
    
    // initialise data
    for(size_t i = 0; i < N; i++)
    {
        v[i] = val;
        xs[i] = val;
    }

    std::cout << "Addresses = " << v.data() << " " << xs << std::endl;

    float triv = trivial_sum(v);

    auto tk1 = myClock::now();
    float k = Kahan_sum(v);
    auto tk2 = myClock::now();

    auto ts1 = myClock::now();
    float S = SIMD_sum(xs, N);
    auto ts2 = myClock::now();

    auto tl1 = myClock::now();
    float STD = std::accumulate(v.begin(), v.end(), 0.0f);
    auto tl2 = myClock::now();

    printf("trivial sum = %f, Kahan sum = %f, SIMD sum = %f; expected 160.0\n", triv, k, S);

    printf("Timing for regular Kahan = %f, and SIMD Kahan = %f, and ratio = %f\n", t_seconds(tk1, tk2), t_seconds(ts1, ts2), t_seconds(tk1, tk2)/t_seconds(ts1, ts2));
    printf("std::accumulate = %f, Timing for std = %f\n", STD, t_seconds(tl1, tl2));
    return 0;
}
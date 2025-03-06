#ifndef SUMMING_H
#define SUMMING_H

// x86 without AVX
//    #include <xmmintrin.h>

// x86 with AVX
    #include <immintrin.h>

// ARM
//    #include <arm_neon.h>


#include <vector>
#include <iostream>
using std::vector;

template<typename T>
T trivial_sum(const vector<T> &vals)
{
    T sum = 0;
    for(auto &x : vals)
    {
        sum = sum + x;
    }
    return sum;
}

template<typename FPtype>
FPtype Kahan_sum(const vector<FPtype> &vals)
{
    FPtype sum = 0;
    FPtype c = 0;
    FPtype y, t;
    for (auto &x : vals)
    {
        y = x - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

//x86 with AVX
    const size_t W = 8;
    using VEC = __m256;
    inline VEC SET_ZERO() { return _mm256_setzero_ps(); }
    inline VEC ADD(VEC a, VEC b) { return _mm256_add_ps(a, b); }
    inline VEC SUB(VEC a, VEC b) { return _mm256_sub_ps(a, b); }
    inline VEC LOAD(const float *p) { return _mm256_load_ps(p); }
    inline void STORE(float* p, VEC a) { return _mm256_store_ps(p, a); }

//x86 without AVX
//    const size_t W = 4;
//    using VEC = __m128;
//    inline VEC SET_ZERO() { return _mm_setzero_ps(); }
//    inline VEC ADD(VEC a, VEC b) { return _mm_add_ps(a, b); }
//    inline VEC SUB(VEC a, VEC b) { return _mm_sub_ps(a, b); }
//    inline VEC LOAD(const float *p) { return _mm_load_ps(p); }
//    inline void STORE(float* p, VEC a) { return _mm_store_ps(p, a); }

//ARM
//    const size_t W = 4;
//    using VEC = float32x4_t;
//    inline VEC SET_ZERO() { return vdupq_n_f32(0); }
//    inline VEC ADD(VEC a, VEC b) { return vaddq_f32(a, b); }
//    inline VEC SUB(VEC a, VEC b) { return vsubq_f32(a, b); }
//    inline VEC LOAD(const float *p) { return vld1q_f32(p); }
//    inline void STORE(float* p, VEC a) { return vst1q_f32(p, a); }


// SIMD Kahan
// size of vals should be multiple of W; can be padded with 0s at the end
// better if vals is appropriately aligned
float SIMD_sum(const float *vals, const size_t N)
{
    size_t N_vec = N / W;
    
    float *sums = new (std::align_val_t(32)) float[W];

    VEC c_vals = SET_ZERO();
    VEC t_vals = SET_ZERO();
    for (size_t i = 0; i < N_vec; i++)
    {
        VEC xs = LOAD(vals + (i*W));  // pointer arithmetic to get the data address
        VEC y_vals = SUB(xs, c_vals);  // y = x - c. Note we don't need to store y_vals in RAM.
        VEC s_vals = t_vals;
        t_vals = ADD(s_vals, y_vals);  // t = sum + y. 
        c_vals = SUB(SUB(t_vals, s_vals), y_vals);  // c = (t-sum) - y
    }
    STORE(sums, t_vals);
    float tot_sum = 0;
    for(int i = 0; i < W; i++)
    {
        tot_sum += sums[i];  // probably don't need compensated mini-sum here because partial sums should be similar size
    }

    return tot_sum;
}

#endif
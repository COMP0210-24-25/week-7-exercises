#ifndef MATRIX_MATHS_H
#define MATRIX_MATHS_H

template<typename T>
void rowMajor_by_rowMajor(const T &A, const T &B, T &C, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                C[i][j] = A[i][k] * B[k][j];
            }
        }
    }
}

template<typename T>
void rowMajor_by_colMajor(const T &A, const T &B, T &C, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                C[i][j] = A[i][k] * B[j][k];
            }
        }
    }
}

template<typename T>
void trivial_transpose(const T &A, T &B, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            B[i][j] = A[j][i];
        }
    }
}

template<typename T>
void trivial_transpose_flat(T &A, T &B, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            B[i*N + j] = A[j*N + i];
        }
    }
}

template<typename T>
void rowMajor_by_rowMajor_flat(const T &A, const T &B, T &C, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                C[i*N + j] = A[i*N + k] * B[k*N + j];
            }
        }
    }
}

template<typename T>
void rowMajor_by_colMajor_flat(const T &A, const T &B, T &C, int N)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                C[i*N + j] = A[i*N + k] * B[j*N + k];
            }
        }
    }
}

#endif
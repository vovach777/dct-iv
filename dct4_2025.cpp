#include <iostream>
#include <fmt/core.h>
#include <cmath>
#include <algorithm>
#include <vector>

const double PI = 3.14159265358979323846;

// Naive DCT4
void dct4(const float * __restrict__   input, float * __restrict__  output, int N) {
    int M = N;  // Выходной размер будет в два раза меньше

    // Преобразование
    for (int k = 0; k < M; ++k) {
        float sum = 0.0f;

        // Косинусное преобразование
        for (int n = 0; n < N; ++n) {
            sum += input[n] * cos(PI / N * (n+0.5) * (k + 0.5) );
        }

        output[k] = sum * sqrt(2./N);
    }
}

// DCT type II, unscaled. Algorithm by Byeong Gi Lee, 1984.
static void forwardTransform(float vec[], float temp[], size_t len) {
    if (len == 1)
        return;
    size_t halfLen = len / 2;
    for (size_t i = 0; i < halfLen; i++) {
        float x = vec[i];
        float y = vec[len - 1 - i];
        temp[i] = x + y;
        temp[i + halfLen] = (x - y) / (std::cos((i + 0.5) * M_PI / len) * 2);
    }
    forwardTransform(temp, vec, halfLen);
    forwardTransform(&temp[halfLen], vec, halfLen);
    for (size_t i = 0; i < halfLen - 1; i++) {
        vec[i * 2 + 0] = temp[i];
        vec[i * 2 + 1] = temp[i + halfLen] + temp[i + halfLen + 1];
    }
    vec[len - 2] = temp[halfLen - 1];
    vec[len - 1] = temp[len - 1];
}

// Пример реализации DCT-IV с использованием DCT-II
void dct4_by_dct2(float* inout, int N) {
    std::vector<float> r(N);
    const float sqrt2_N = std::sqrt(2.0f / N);

    // (a) Pre-processing: r(n) = 2 * u(n) * cos(pi * (2 * n + 1) / (4 * N))
    for (int n = 0; n < N; ++n) {
        inout[n] = 2.0f * inout[n] * std::cos(M_PI * (2 * n + 1) / (4 * N));
    }

    // (b) Вычисление DCT-II (на месте)
        
    forwardTransform(inout, r.data(), N);

    inout[0] /= 2;
    
    // (c) Post-processing: Y4(k) = Y2(k) - Y4(k-1), Y4(-1) = Y4(0)
    r[0] = inout[0];
    for (int k = 1; k < N; ++k) {
        r[k] = inout[k] - r[k - 1];
    }

    // Копируем результат обратно
    for (int k = 0; k < N; ++k) {
        inout[k] = r[k] * sqrt2_N;
    }    
}



template <typename T>
void print_mat(T&& t) {
    for (int y = 0; y < 8; ++y, fmt::print("\n"))
    for (int x = 0; x < 8; ++x)
    {
        //fmt::print(" {:4}",int( t[y*8+x]  + 0.5f - (t[y*8+x] < 0.0f) ));
        fmt::print(" {:4}", std::lround( t[y*8+x]  ));
    }
    fmt::print("\n");
}

template <typename T>
void print_mat(T&& t, int SZ) {
    for (int y = 0; y < SZ; ++y, fmt::print("\n"))
    for (int x = 0; x < SZ; ++x)
    {
        //fmt::print(" {:4}",int( t[y*8+x]  + 0.5f - (t[y*8+x] < 0.0f) ));
        fmt::print(" {:4}", std::lround( t[y*SZ+x]  ));
    }
    fmt::print("\n");
}

void fill_mat(float* S, int SZ){
    for (int y=0; y < SZ; ++y)
    for (int x=0; x < SZ; ++x)
       S[y*SZ+x] = y*10+x;
}

int main() {

    float s[8 * 8];
    float S[8 * 8];
    fill_mat(S,8);
    print_mat(S,8);

    dct4(S,s,64);
    print_mat(s,8);

    fill_mat(s,8);
    dct4_by_dct2(s,64);
    print_mat(s,8);

    dct4_by_dct2(s,64);
    print_mat(s,8);

    return 0;
}

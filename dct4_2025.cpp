#include <iostream>
#include <fmt/core.h>
#include <cmath>
#include <algorithm>
#include <vector>

const double PI = 3.14159265358979323846;

const float G[64] = {
    -415.38, -30.19, -61.20, 27.24, 56.12, -20.10, -2.39, 0.46,
    4.47, -21.86, -60.76, 10.25, 13.15, -7.09, -8.54, 4.88,
    -46.83, 7.37, 77.13, -24.56, -28.91, 9.93, 5.42, -5.65,
    -48.53, 12.07, 34.10, -14.76, -10.24, 6.30, 1.83, 1.95,
    12.12, -6.55, -13.20, -3.95, -1.87, 1.75, -2.79, 3.14,
    -7.73, 2.91, 2.38, -5.94, -2.38, 0.94, 4.30, 1.85,
    -1.03, 0.18, 0.42, -2.42, -0.88, -3.02, 4.12, -0.66,
    -0.17, 0.14, -1.07, -4.19, -1.17, -0.10, 0.50, 1.68
};



// Функция для расчета коэффициента C
inline float C(int k) {
    return (k == 0) ? 1.0f / std::sqrt(2.0f) : 1.0f;
}

// Обратное DCT для матрицы 8x8, представленной одномерным массивом
static void inverseDCT(float* S) {
    float s[64];
    for (int y = 0; y < 8; ++y) { // Строки
        for (int x = 0; x < 8; ++x) { // Колонки
            float sum = 0.0f;
            for (int u = 0; u < 8; ++u) { // Координата частоты по строкам
                for (int v = 0; v < 8; ++v) { // Координата частоты по колонкам
                    float Cu = C(u), Cv = C(v);
                    float Suv = S[u * 8 + v]; // Доступ к элементу входной матрицы
                    sum += Cu * Cv * Suv *
                           std::cos((2 * y + 1) * u * PI / 16) * // Угол для строк
                           std::cos((2 * x + 1) * v * PI / 16);  // Угол для колонок
                }
            }
            s[y * 8 + x] = 0.25f * sum; // Сохраняем результат
        }
    }
    std::copy(s,s+64,S);    
}

#include <cmath>
#include <algorithm>

// // Обратное DCT для матрицы NxN, представленной одномерным массивом
// template<int N>
// static void inverseDCT(float* S) {
//     float s[N * N]; // Временный буфер для хранения результата
//     constexpr float scale = 0.25f; // Общий масштаб, фиксированный для всех N

//     for (int y = 0; y < N; ++y) { // Строки
//         for (int x = 0; x < N; ++x) { // Колонки
//             float sum = 0.0f;
//             for (int u = 0; u < N; ++u) { // Координата частоты по строкам
//                 for (int v = 0; v < N; ++v) { // Координата частоты по колонкам
//                     float Cu = C(u), Cv = C(v);
//                     float Suv = S[u * N + v]; // Доступ к элементу входной матрицы
//                     sum += Cu * Cv * Suv *
//                            std::cos((2 * y + 1) * u * PI / (2 * N)) * // Угол для строк
//                            std::cos((2 * x + 1) * v * PI / (2 * N));  // Угол для колонок
//                 }
//             }
//             s[y * N + x] = scale * sum; // Сохраняем результат
//         }
//     }
//     std::copy(s, s + N * N, S); // Копируем результат обратно в исходный массив
// }


static void transpose(float* o, int N=8) {
    for (int y = 0; y < N; ++y) {
        for (int x = y + 1; x < N; ++x) {
            float temp = o[y * N + x];
            o[y * N + x] = o[x * N + y];
            o[x * N + y] = temp;
        }
    }    
}


static void inverseTransform(float *vec, float *temp, size_t len) {
	if (len == 1)
		return;
	size_t halfLen = len / 2;
	temp[0] = vec[0];
	temp[halfLen] = vec[1];
	for (size_t i = 1; i < halfLen; i++) {
		temp[i] = vec[i * 2];
		temp[i + halfLen] = vec[i * 2 - 1] + vec[i * 2 + 1];
	}
	inverseTransform(temp, vec, halfLen);
	inverseTransform(&temp[halfLen], vec, halfLen);
	for (size_t i = 0; i < halfLen; i++) {
		float x = temp[i];
		float y = temp[i + halfLen] / float( (std::cos((i + 0.5f) * M_PI / len) * 2) );        
		vec[i] = x + y;
		vec[len - 1 - i] = x - y;
	}
}

// Функция для MDCT
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


static void normalize(float * vec, int stride=8)
{
    for (int j = 0; j < 8; ++j)
    for (int i = 0; i < 8; ++i)
    {        
        // auto jc = C() (j == 0) ? (1 / (std::sqrt(2.0f)*2)) : (1/2.0f);
        // auto ic =  (i == 0) ? (1 / (std::sqrt(2.0f)*2)) : (1/2.0f);
        // vec[j*stride+i] *= jc * ic;
        vec[j*stride+i] *= C(i) * C(j) * (1.0f / 4);
    }
}

// static void inverseTransform(float* vec, int N )
// {
//     alignas(32) float temp[N];
//     normalize(vec,N);
    
//     for (int dim = 0; dim < 2; ++dim)
//     {   
//         for (int i = 0; i < N; ++i)
//         {
//             inverseTransform(vec+i*N, temp,N);
//         }
//         transpose(vec,N);

//     }
// }


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


// // Функция для MDCT
// void mdct(const float * __restrict__   input, float * __restrict__  output, int N) {
//     int M = N / 2;  // Выходной размер будет в два раза меньше

//     // Преобразование
//     for (int k = 0; k < M; ++k) {
//         float sum = 0.0f;

//         // Косинусное преобразование
//         for (int n = 0; n < N; ++n) {
//             sum += input[n] * cos(PI / N * (n + 0.5f) * (k + 0.5f));
//         }

//         output[k] = sum;
//     }
// }

// // Функция для IMDCT
// void imdct(const float * __restrict__   input, float * __restrict__  output, int N) {
//     int M = N / 2;  // Входной размер будет в два раза меньше

//     // Обратное преобразование
//     for (int n = 0; n < N; ++n) {
//         float sum = 0.0f;

//         // Косинусное преобразование
//         for (int k = 0; k < M; ++k) {
//             sum += input[k] * cos(PI / N * (n + 0.5f) * (k + 0.5f));
//         }

//         output[n] = sum;
//     }
// }

void fill_mat(float* S, int SZ){
    for (int y=0; y < SZ; ++y)
    for (int x=0; x < SZ; ++x)
       S[y*SZ+x] = y*10+x;
}

int main() {

    float s[8 * 8]; // Результирующий массив
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

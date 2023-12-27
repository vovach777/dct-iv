#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

namespace dct_iv_implementation {
template <typename complex_t>
void fft0(int n, int s, bool eo, complex_t* x, complex_t* y)
// n  : sequence length
// s  : stride
// eo : x is output if eo == 0, y is output if eo == 1
// x  : input sequence(or output sequence if eo == 0)
// y  : work area(or output sequence if eo == 1)
{
    const int m = n / 2;
    const decltype(x->real()) theta0 = 2 * M_PI / n;

    if (n == 1) {
        if (eo)
            for (int q = 0; q < s; q++) y[q] = x[q];
    } else {
        for (int p = 0; p < m; p++) {
            const complex_t wp = complex_t(cos(p * theta0), -sin(p * theta0));
            for (int q = 0; q < s; q++) {
                const complex_t a = x[q + s * (p + 0)];
                const complex_t b = x[q + s * (p + m)];
                y[q + s * (2 * p + 0)] = a + b;
                y[q + s * (2 * p + 1)] = (a - b) * wp;
            }
        }
        fft0(n / 2, 2 * s, !eo, y, x);
    }
}

template <typename T>
void fft(std::vector<std::complex<T>> x)  // Fourier transform
// n : sequence length
// x : input/output sequence
{
    auto n = std::distance(begin, end);
    std::vector<std::complex<T>> y(n);
    fft0(n, 1, 0, x.data(), y.data());
    // for (auto & v : x)
    //    v /= n;
}

template <typename T>
void dct4(std::vector<T>& audio_signal) {
    int window_length = audio_signal.size();

    std::vector<std::complex<T>> audio_dct(8 * window_length,
                                           std::complex<T>(0));

    for (int i = 1, j = 0; i < 2 * window_length && j < audio_signal.size();
         i += 2, j += 1) {
        audio_dct[i] = audio_signal[j];
    }

    for (int i = 2 * window_length + 1, j = audio_signal.size() - 1;
         i < 4 * window_length && j >= 0; i += 2, j -= 1) {
        audio_dct[i] = -audio_signal[j];
    }

    for (int i = 4 * window_length + 1, j = 0;
         i < 6 * window_length && j < audio_signal.size(); i += 2, j += 1) {
        audio_dct[i] = -audio_signal[j];
    }

    for (int i = 6 * window_length + 1, j = audio_signal.size() - 1;
         i < 8 * window_length && j >= 0; i += 2, j -= 1) {
        audio_dct[i] = audio_signal[j];
    }

    fft(audio_dct);

    auto scale = T(std::sqrt(2.0 / window_length) / 4.0);
    for (int i = 1, j = 0; i < 2 * window_length && j < window_length; i += 2, j += 1) {
        audio_signal[j] = std::real(audio_dct[i]) * scale;
    }

}
}

using dct_iv_implementation::dct4;

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

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

template <typename Iter>
auto fft(Iter begin, Iter end)  // Fourier transform
// n : sequence length
// x : input/output sequence
{
    auto n = std::distance(begin, end);
    using T = typename std::remove_reference<decltype(*begin)>::type;
    std::vector<T> x(begin, end);
    std::vector<T> y(n);
    fft0(n, 1, 0, x.data(), y.data());
    // for (auto & v : x)
    //    v /= n;
    return x;
}

template <typename Container>
auto fft(Container&& container)  // Fourier transform
{
    return fft(container.begin(), container.end());
}

template <typename container>
void print_line(container&& C) {
    for (auto c : C) {
        std::cout << std::fixed << std::setprecision(4) << std::setw(8) << double(c);
    }
    std::cout << std::endl;
}

template <typename T>
auto dctIV(const std::vector<T>& audio_signal) {
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

    audio_dct = fft(audio_dct);

    std::vector<T> result;
    result.reserve(window_length);
    auto scale = T(std::sqrt(2.0 / window_length) / 4.0);
    for (int i = 1; i < 2 * window_length; i += 2) {
        result.push_back(std::real(audio_dct[i]) * scale);
    }

    return result;
}

int main() {
    std::vector<float> audio_signal{1, 2, 3, 4, 5, 6, 7, 8};

    auto tmp = audio_signal;
    print_line(tmp);
    tmp = dctIV(tmp);
    print_line(tmp);
    tmp = dctIV(tmp);
    print_line(tmp);

    return 0;
}

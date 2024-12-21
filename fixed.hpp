#pragma once

#include <bits/stdc++.h>

namespace types {

template <std::size_t N, std::size_t K, bool Fast = false>
struct Fixed {
    using type = std::conditional_t<
        (N == 8 && Fast == false), int8_t,
        std::conditional_t<
            (N == 16 && Fast == false), int16_t,
            std::conditional_t<
                (N == 32 && Fast == false), int32_t,
                std::conditional_t<
                    (N == 64 && Fast == false), int64_t,
                        std::conditional_t<
                            (N <= 8 && Fast == true), int_fast8_t,
                            std::conditional_t<
                                (N <= 16 && Fast == true), int_fast16_t,
                                std::conditional_t<
                                    (N <= 32 && Fast == true), int_fast32_t,
                                    std::conditional_t<
                                        (N <= 64 && Fast == true), int_fast64_t,
                                        void
                                >
                            >
                        >
                    >
                >
            >
        >
    >;
    static constexpr std::size_t kNValue = N;
    static constexpr std::size_t kKValue = K;
    static constexpr bool kFast          = Fast;

    constexpr Fixed(int v) : v(static_cast<int64_t>(v) << K) {}
    constexpr Fixed(float f) : v(static_cast<int64_t>(f * (type(1) << K))) {}
    constexpr Fixed(double f) : v(static_cast<int64_t>(f * (type(1) << K))) {}
    constexpr Fixed() : v(0) {}
    template <int _n, int _m>
    constexpr Fixed(Fixed<_n, _m> a)
    {
        v = a.v;
        if (K > _m)
        {
            v >> (K - _m);
        }
        else if (K < _m)
        {
            v << (_m - K);
        }
    }

    static constexpr Fixed from_raw(type x){
        Fixed ret;
        ret.v = x;
        return ret;
    }
    type v;
    constexpr Fixed operator+(const Fixed& other) const {
        return Fixed::from_raw(v + other.v);
    }

    constexpr Fixed operator-(const Fixed& other) const {
        return Fixed::from_raw(v - other.v);
    }

    constexpr Fixed operator*(const Fixed& other) const {
        return Fixed::from_raw((v * other.v) >> K);
    }

    constexpr Fixed operator/(const Fixed& other) const {
        return Fixed::from_raw((v << K) / other.v);
    }

    constexpr Fixed& operator+=(const Fixed& other) {
        v += other.v;
        return *this;
    }

    constexpr const Fixed& operator-=(const Fixed& other) {
        v -= other.v;
        return *this;
    }

    constexpr Fixed& operator*=(const Fixed& other) {
        v = (v * other.v) >> K; // Умножение с учетом масштаба
        return *this;
    }

    constexpr Fixed& operator/=(const Fixed& other) {
        v = (v << K) / other.v; // Деление с учетом масштаба
        return *this;
    }

    // Перегрузка операторов для работы с целыми числами
    constexpr Fixed operator+(int const &other) const {
        return *this + Fixed(other);
    }

    constexpr Fixed operator-(int const &other) const {
        return *this - Fixed(other);
    }

    constexpr Fixed operator*(int const &other) const {
        return Fixed::from_raw(v * other);
    }

    constexpr Fixed operator/(int const &other) const {
        return Fixed::from_raw(v / other);
    }

    // Перегрузка операторов для работы с плавающими числами
    constexpr Fixed operator+(float const &other) const {
        return *this + Fixed(other);
    }

    constexpr Fixed operator-(float const &other) const {
        return *this - Fixed(other);
    }

    constexpr Fixed operator*(float const &other) const {
        return Fixed::from_raw(v * Fixed(other).v);
    }

    constexpr Fixed operator/(float const &other) const {
        return Fixed::from_raw(v / Fixed(other).v);
    }


    constexpr Fixed operator+(double const &other) const {
        return *this + Fixed(other);
    }

    constexpr Fixed operator-(double const &other) const {
        return *this - Fixed(other);
    }

    constexpr Fixed operator*(double const &other) const {
        return Fixed::from_raw(v * Fixed(other).v);
    }

    constexpr Fixed operator/(double const &other) const {
        return Fixed::from_raw(v / Fixed(other).v);
    }
    constexpr Fixed operator+=(double const &other){
        return (*this) = *this + Fixed(other);
    }

    constexpr Fixed operator-=(double const &other){
        return (*this) = *this - Fixed(other);
    }

    constexpr Fixed operator*=(double const &other){
        return (*this) = Fixed::from_raw(v * Fixed(other).v);
    }

    constexpr Fixed operator/=(double const &other){
        return (*this) = Fixed::from_raw(v / Fixed(other).v);
    }



    template <int _n, int _m>
    auto operator<=>(const Fixed<_n, _m>& a) const
    {
        return (*this).v <=> equate(*this, a);
    }
    template <int _n, int _m>
    bool operator==(const Fixed<_n, _m>& a) const
    {
        return (*this).v == equate(*this, a);
    }
    auto operator<=>(const Fixed& a) const
    {
        return (*this).v <=> a.v;
    }
    bool operator==(const Fixed& a) const
    {
        return (*this).v == a.v;
    }
    bool operator==(const double& f) const
    {
        return (*this).v == f * (type(1) << K);
    }
    auto operator<=>(const double& f) const
    {
        return (*this).v <=> f * (type(1) << K);
    }
    bool operator==(const int& f) const
    {
        return (*this).v == f * (type(1) << K);
    }
    auto operator<=>(const int& f) const
    {
        return (*this).v <=> f * (type(1) << K);
    }
    bool operator==(const float& f) const
    {
        return (*this).v == f * (type(1) << K);
    }
    auto operator<=>(const float& f) const
    {
        return (*this).v <=> f * (type(1) << K);
    }
    explicit operator double() const {return v / (double) ((int64_t)(1) << K);}
    explicit operator float() const {return v / (float) ((int64_t)(1) << K);}
    explicit operator int() const {return (int64_t)v >> K;}
};
template<size_t N, size_t M, size_t K, size_t D>
int64_t equate(Fixed<N, M> const &a, Fixed<K, D> const &b)
{
    int64_t _v = b.v;
    if (M > D)
    {
        _v >> (M - D);
    }
    else if (M < D)
    {
        _v << (D - M);
    }
    return _v;
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> operator+(Fixed<N, M> const &a, Fixed<K, D> const &b) {
    return a+Fixed<N, M>(b);
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> operator-(Fixed<N, M> const &a, Fixed<K, D> const &b) {
    return a-Fixed<N, M>(b);
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> operator*(Fixed<N, M> const &a, Fixed<K, D> const &b) {
    return a*Fixed<N, M>(b);
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> operator/(Fixed<N, M> const &a, Fixed<K, D> const &b) {
    return a/Fixed<N, M>(b);
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> &operator+=(Fixed<N, M> &a, Fixed<K, D> b) {
    return a = a + b;
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> &operator-=(Fixed<N, M> &a, Fixed<K, D> b) {
    return a = a - b;
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> &operator*=(Fixed<N, M> &a, Fixed<K, D> b) {
    return a = a * b;
}
template<size_t N, size_t M, size_t K, size_t D>
Fixed<N, M> &operator/=(Fixed<N, M> &a, Fixed<K, D> b) {
    return a = a / b;
}
template<size_t N, size_t M>
Fixed<N, M> abs(Fixed<N, M> x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}
template<size_t N, size_t M>
ostream &operator<<(ostream &out, Fixed<N, M> x) {
    return out << x.v / (double) ((int64_t)1 << M);
}
template<size_t N, size_t M>
double operator+(double b, Fixed<N, M> a) {
    return ((double)a)+b;
}
template<size_t N, size_t M>
double operator-(double b, Fixed<N, M> a) {
    return b-((double)a);
}
template<size_t N, size_t M>
double operator*(double b, Fixed<N, M> a) {
    return ((double)a)*b;
}
template<size_t N, size_t M>
double operator/(double b, Fixed<N, M> a) {
    return b/((double)a);
}
template<size_t N, size_t M>
double operator+=(double &b, Fixed<N, M> a) {
    return b=b+((double)a);
}
template<size_t N, size_t M>
double operator-=(double &b, Fixed<N, M> a) {
    return b=b-((double)a);
}
template<size_t N, size_t M>
double operator*=(double &b, Fixed<N, M> a) {
    return b=b*((double)a);
}
template<size_t N, size_t M>
double operator/=(double &b, Fixed<N, M> a) {
    return b=b/((double)a);
}
template<size_t N, size_t M>
float operator+(float b, Fixed<N, M> a) {
    return ((float)a)+b;
}
template<size_t N, size_t M>
float operator-(float b, Fixed<N, M> a) {
    return b-((float)a);
}
template<size_t N, size_t M>
float operator*(float b, Fixed<N, M> a) {
    return ((float)a)*b;
}
template<size_t N, size_t M>
float operator/(float b, Fixed<N, M> a) {
    return b/((float)a);
}
template<size_t N, size_t M>
float operator+=(float &b, Fixed<N, M> a) {
    return b=b+((float)a);
}
template<size_t N, size_t M>
float operator-=(float &b, Fixed<N, M> a) {
    return b=b-((float)a);
}
template<size_t N, size_t M>
float operator*=(float &b, Fixed<N, M> a) {
    return b=b*((float)a);
}
template<size_t N, size_t M>
float operator/=(float &b, Fixed<N, M> a) {
    return b=b/((float)a);
}



template <std::size_t N, std::size_t K>
using FastFixed = Fixed<N, K, /* Fast */ true>;

}  // namespace types

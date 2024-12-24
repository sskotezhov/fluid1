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
    Fixed operator+(Fixed const &b) {
        return Fixed::from_raw(v+b.v);
    }
    Fixed operator-(Fixed const &b) {
        return Fixed::from_raw(v-b.v);
    }
    Fixed operator*(Fixed const &b) {
        return Fixed::from_raw(((int64_t)v*b.v) >> K);
    }
    Fixed operator/(Fixed const &b) {
        return Fixed::from_raw(((int64_t)v << K)/b.v);
    }
    Fixed& operator+=(Fixed const &b) {
        return (*this) = (*this) + b;
    }
    Fixed& operator-=(Fixed const &b) {
        return (*this) = (*this) - b;
    }
    Fixed& operator*=(Fixed const &b) {
        return (*this) = (*this) * b;
    }
    Fixed& operator/=(Fixed const &b) {
        return (*this) = (*this) / b;
    }
    constexpr Fixed(int v) : v(static_cast<int64_t>(v) << K) {}
    constexpr Fixed(float f) : v(static_cast<int64_t>(f * (type(1) << K))) {}
    constexpr Fixed(double f) : v(static_cast<int64_t>(f * (type(1) << K))) {}
    constexpr Fixed() : v(0) {}
    template <size_t _n, size_t _m, bool FAST>
    explicit constexpr Fixed(Fixed<_n, _m, FAST> a)
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
    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
    explicit operator double() const {return v / (double) ((int64_t)(1) << K);}
    explicit operator float() const {return v / (float) ((int64_t)(1) << K);}
    explicit operator int() const {return (int64_t)v >> K;}
};
template<size_t N, size_t M, size_t _N, size_t _M, bool FAST>
auto operator<=>(const Fixed<N, M, FAST>& a, const Fixed<_N,_M, FAST>& b)
{
    return a<=>Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M, bool FAST>
bool operator==(const Fixed<N, M, FAST>& a, const Fixed<_N,_M, FAST>& b)
{
    return a==Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M>
auto operator<=>(const Fixed<N, M>& a, const Fixed<N,M,true>& b)
{
    return a<=>Fixed<N,M>(b);
}
template<size_t N, size_t M>
bool operator==(const Fixed<N, M>& a, const Fixed<N,M,true>& b)
{
    return a==Fixed<N,M>(b);
}
template<size_t N, size_t M>
auto operator<=>(const Fixed<N, M, true>& a, const Fixed<N,M>& b)
{
    return a<=>Fixed<N,M,true>(b);
}
template<size_t N, size_t M>
bool operator==(const Fixed<N, M, true>& a, const Fixed<N,M>& b)
{
    return a==Fixed<N,M,true>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M>
auto operator<=>(const Fixed<N, M>& a, const Fixed<_N,_M,true>& b)
{
    return a<=>Fixed<N,M>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M>
bool operator==(const Fixed<N, M>& a, const Fixed<_N,_M,true>& b)
{
    return a==Fixed<N,M>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M>
auto operator<=>(const Fixed<N, M, true>& a, const Fixed<_N,_M>& b)
{
    return a<=>Fixed<N,M,true>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M>
bool operator==(const Fixed<N, M, true>& a, const Fixed<_N,_M>& b)
{
    return a==Fixed<N,M,true>(b);
}
template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator+(Fixed<N, M, FAST1> a, Fixed<_N, _M, FAST2> b) {
    return a+Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator+(Fixed<N, M, FAST1> a, Fixed<N, M, FAST2> b) {
    return a+Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> operator+(Fixed<N, M, FAST> a, type_b b) {
    return a+Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b operator+(type_b b, Fixed<N, M, FAST> a) {
    return type_b(Fixed<N,M, FAST>(b)+a);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator-(Fixed<N, M, FAST1> a, Fixed<_N, _M, FAST2> b) {
    return a-Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator-(Fixed<N, M, FAST1> a, Fixed<N, M, FAST2> b) {
    return a-Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> operator-(Fixed<N, M, FAST> a, type_b b) {
    return a-Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b operator-(type_b b, Fixed<N, M, FAST> a) {
    return type_b(Fixed<N,M,FAST>(b)-a);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator*(Fixed<N, M, FAST1> a, Fixed<_N, _M, FAST2> b) {
    return a*Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator*(Fixed<N, M, FAST1> a, Fixed<N, M, FAST2> b) {
    return a*Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> operator*(Fixed<N, M, FAST> a, type_b b) {
    return a*Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b operator*(type_b b, Fixed<N, M, FAST> a) {
    return type_b(a*Fixed<N,M, FAST>(b));
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator/(Fixed<N, M, FAST1> a, Fixed<_N, _M, FAST2> b) {
    return a/Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> operator/(Fixed<N, M, FAST1> a, Fixed<N, M, FAST2> b) {
    return a/Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> operator/(Fixed<N, M, FAST> a, type_b b) {
    return a/Fixed<N,M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b operator/(type_b b, Fixed<N, M, FAST> a) {
    return type_b(Fixed<N,M, FAST>(b)/a);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator+=(Fixed<N, M, FAST1> &a, Fixed<_N, _M, FAST2> b) {
    return a = a + Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator+=(Fixed<N, M, FAST1> &a, Fixed<N, M, FAST2> b) {
    return a = a + Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b &operator+=(type_b &b, Fixed<N, M, FAST> a) {
    return b = b + type_b(a);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> &operator+=(Fixed<N, M, FAST> &a, type_b b) {
    return a = a + Fixed<N, M, FAST>(b);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator-=(Fixed<N, M, FAST1> &a, Fixed<_N, _M, FAST2> b) {
    return a = a - Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator-=(Fixed<N, M, FAST1> &a, Fixed<N, M, FAST2> b) {
    return a = a - Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> &operator-=(Fixed<N, M, FAST> &a, type_b b) {
    return a = a - Fixed<N, M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b &operator-=(type_b &b, Fixed<N, M, FAST> a) {
    return b = b - type_b(a);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator*=(Fixed<N, M, FAST1> &a, Fixed<_N, _M, FAST2> b) {
    return a = a * Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator*=(Fixed<N, M, FAST1> &a, Fixed<N, M, FAST2> b) {
    return a = a * Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST>& operator*=(Fixed<N, M, FAST> &a, type_b b) {
    return a = a * Fixed<N, M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b &operator*=(type_b &b, Fixed<N, M, FAST> a) {
    return b = b * type_b(a);
}

template<size_t N, size_t M, size_t _N, size_t _M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator/=(Fixed<N, M, FAST1> &a, Fixed<_N, _M, FAST2> b) {
    return a = a / Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST1, bool FAST2>
Fixed<N, M, FAST1> &operator/=(Fixed<N, M, FAST1> &a, Fixed<N, M, FAST2> b) {
    return a = a / Fixed<N,M, FAST1>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
Fixed<N, M, FAST> &operator/=(Fixed<N, M, FAST> &a, type_b b) {
    return a = a / Fixed<N, M, FAST>(b);
}
template<size_t N, size_t M, bool FAST, typename type_b>
type_b &operator/=(type_b &b, Fixed<N, M, FAST> a) {
    return b = b / type_b(a);
}


template<size_t N, size_t M, bool FAST>
std::ostream &operator<<(std::ostream &out, Fixed<N, M, FAST> x) {
    return out << x.v / (double) (1 << M);
}
template<size_t N, size_t M, bool FAST>
Fixed<N, M> abs(Fixed<N, M, FAST> x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}
template <std::size_t N, std::size_t K>
using FastFixed = Fixed<N, K, /* Fast */ true>;

}  // namespace types
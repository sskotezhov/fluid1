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
    auto operator<=>(const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
    explicit operator double() const {return v / (double) ((int64_t)(1) << K);}
    explicit operator float() const {return v / (float) ((int64_t)(1) << K);}
    explicit operator int() const {return (int64_t)v >> K;}
};
template<size_t N, size_t M>
Fixed<N, M> operator-(Fixed<N, M> const &a, Fixed<N, M> const &b) {
    return Fixed<N, M>::from_raw(a.v-b.v);
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> operator-(Fixed<N, M> a, type_b b) {
    return a-Fixed<N,M>(b);
}
template<size_t N, size_t M, typename type_b>
type_b operator-(type_b b, Fixed<N, M> a) {
    return type_b(Fixed<N,M>(b)-a);
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> operator/(Fixed<N, M> a, type_b b) {
    return a/Fixed<N,M>(b);
}
template<size_t N, size_t M>
Fixed<N, M> operator+(Fixed<N, M> const &a, Fixed<N, M> const &b) {
    return Fixed<N, M>::from_raw(a.v+b.v);
}
template<size_t N, size_t M>
Fixed<N, M> operator*(Fixed<N, M> const &a, Fixed<N, M> const &b) {
    return Fixed<N, M>::from_raw(((int64_t)a.v*b.v) >> M);
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> operator*(Fixed<N, M> a, type_b b) {
    return a*Fixed<N,M>(b);
}
template<size_t N, size_t M, typename type_b>
type_b operator*(type_b b, Fixed<N, M> a) {
    return type_b(a*Fixed<N,M>(b));
}
template<size_t N, size_t M>
Fixed<N, M> operator/(Fixed<N, M> a, Fixed<N, M> b) {
    return Fixed<N, M>::from_raw(((int64_t)a.v<<M)/b.v);
}
template<size_t N, size_t M>
Fixed<N, M> &operator+=(Fixed<N, M> &a, Fixed<N, M> b) {
    return a = a + b;
}
template<size_t N, size_t M, typename type_b>
type_b &operator+=(type_b &b, Fixed<N, M> a) {
    return b = b + type_b(a);
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> &operator+=(Fixed<N, M> &a, type_b b) {
    return a = a + Fixed<N, M>(b);
}
template<size_t N, size_t M>
Fixed<N, M> &operator-=(Fixed<N, M> &a, Fixed<N, M> b) {
    return a = a - b;
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> &operator-=(Fixed<N, M> &a, type_b b) {
    return a = a - Fixed<N, M>(b);
}
template<size_t N, size_t M, typename type_b>
type_b &operator-=(type_b &b, Fixed<N, M> a) {
    return b = b - type_b(a);
}
template<size_t N, size_t M, typename type_b>
Fixed<N, M> operator*=(Fixed<N, M> &a, type_b b) {
    return a = a * Fixed<N, M>(b);
}
template<size_t N, size_t M>
Fixed<N, M> &operator*=(Fixed<N, M> &a, Fixed<N, M> b) {
    return a = a * b;
}
template<size_t N, size_t M>
Fixed<N, M> &operator/=(Fixed<N, M> &a, Fixed<N, M> b) {
    return a = a / b;
}
template<size_t N, size_t M>
ostream &operator<<(ostream &out, Fixed<N, M> x) {
    return out << x.v / (double) (1 << 16);
}
template<size_t N, size_t M>
Fixed<N, M> abs(Fixed<N, M> x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}
template <std::size_t N, std::size_t K>
using FastFixed = Fixed<N, K, /* Fast */ true>;

}  // namespace types

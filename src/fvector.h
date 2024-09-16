#pragma once
#include <cmath>
#include <iostream>
#include <vector>

// Fixed size vector.  Size of vector must be known at compile time.
template <int N>
class FVector {
public:
    FVector() {}
    FVector(std::initializer_list<double> x);
    double  operator()(int i) const { return _x[i]; }
    double& operator()(int i) { return _x[i]; }
    double* begin() { return _x; }
    const double* begin() const { return _x; }
    double* end() { return _x+N; }
    const double* end() const { return _x+N; }
    int size() const { return N; }

    // Adds scaled-vector (v*s) to this vector.
    void scaled_add(const FVector &v, double s=1.0) {
        for (int i=0; i<N; ++i) _x[i] += v(i)*s;
    }

private:
    double _x[N];
};

template <int N>
FVector<N>::FVector(std::initializer_list<double> x) {
  std::copy(x.begin(), x.end(), _x);
}

//! Adds to fixed size vectors together.
template<int N>
FVector<N> operator+(FVector<N> x, const FVector<N> &y) {
    x += y;
    return x;
}

//! Adds a double value into vector.
template<int N>
FVector<N>& operator+=(FVector<N> &x, double y) {
    x += y;
    return x;
}        

//! In place addition operator.
template<int N>
FVector<N>& operator+=(FVector<N> &x, const FVector<N> &y) {
    for (int i=0; i<N; ++i) x(i) += y(i);
    return x;
}

template<int N>
FVector<N> operator-(FVector<N> x, const FVector<N> &y) {
    x -= y;
    return x;
}

template<int N>
FVector<N>& operator-=(FVector<N> &x, const FVector<N> &y) {
    for (int i=0; i<N; ++i) x(i) -= y(i);
    return x;
}

template<int N>
FVector<N> operator-=(FVector<N> &x, double y) {
    for (int i=0; i<N; ++i) x(i) -= y;
    return x;
}

template<int N>
FVector<N> operator-(FVector<N> x, double y) {
    x -= y;
    return x;
}

// Scalar-vector product. 
template<int I>
FVector<I> operator*(const double a, FVector<I> x) {
    for (auto &v: x) v *= a;
    return x;
}

template<int I>
FVector<I> operator*(const FVector<I> &x, const double a) {
    return a*x;
}

template<int I>
FVector<I> operator/(FVector<I> x, double a) {
    x *= (1.0/a);
    return x;
}

template<int I>
FVector<I>& operator*=(FVector<I> &x, double a) {
    for (int i=0; i<I; ++i) x(i) *= a;
    return x;
}

template <int N>
double dot(const FVector<N> &x, const FVector<N> &y) {
    double r = 0.0;
    for (int i=0; i<N; ++i) r += x(i)*y(i);
    return r;
}

// Returns a statically sized vector of evenly spaced points (like MATLAB).
template <int N>
FVector<N> linspace(double start, double stop) {
    double dt = (stop-start)/double(N-1);
    FVector<N> x;
    for (int i=0; i<N; ++i) {
        x(i) = start + double(i)*dt;
    }
    x(N-1) = stop;
    return x;
}

template <int I>
double sum(const FVector<I> &x) {
    double s = 0.0;
    for (auto _x: x) s += _x;
    return s;
}

template <int I>
double norm(const FVector<I> &x) {
    double s = 0.0;
    for (auto _x: x) s += _x*_x;
    return sqrt(s);
}

// Returns the sum of an std::vector of doubles.
inline double sum(const std::vector<double> &x) {
    double s = 0.0;
    for (auto _x: x) s += _x;
    return s;
}

// Returns sign function of a vector.
template <int N>
FVector<N> sign(FVector<N> v) {
    for (auto &x: v) {
        x = x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
    }
    return v;
}

// Returns heaviside function of a vector.
template <int I>
inline FVector<I> step(const FVector<I> &v) {
    FVector<4> sv;
    for (int i=0; i<I; i++) {
        if (v(i)>0)      sv(i) = 1.0;
        else             sv(i) = 0.0;
    }
    return sv;
}

inline double dot(const std::vector<double> &x, const std::vector<double> &y) {
    double r = 0;
    for (int i=0; i<x.size(); ++i) r += x[i]*y[i];
    return r;
}

//! Computes a cross product for R3 vectors.
inline FVector<3> cross(const FVector<3> &x, const FVector<3> &y) {
    return FVector<3> {x(1)*y(2)-x(2)*y(1),
                       x(2)*y(0)-x(0)*y(2),
                       x(0)*y(1)-x(1)*y(0)};
}

//! Returns a vector full of zeros.
template <int N>
FVector<N> zeros() {
    FVector<N> v;
    for (auto &x: v) x = 0.0;
    return v;
}

//! Each element of vector generated from a uniform random distribution [0,1].
template <int N>
FVector<N> rand() {
    FVector<N> v;
    for (auto &x: v) x = double(rand()) / double(RAND_MAX);
    return v;
}

//! Generic template for a N-D unit vector with random orientation.
template<int N> FVector<N> rand_unit();
//! Generates a 2D unit vector with random orientation.
template <> inline FVector<2> rand_unit() {
    FVector<2> v;
    const auto pi = acos(-1.0);
    auto theta = 2.0*pi*double(rand())/double(RAND_MAX);
    return {cos(theta), sin(theta)};
}
//! Generates a 3D unit vector with random orientation.
template <> inline FVector<3> rand_unit() {
    FVector<3> v;
    const auto pi = acos(-1.0);
    auto theta = 2.0*pi*double(rand())/double(RAND_MAX);
    auto phi   = 2.0*pi*double(rand())/double(RAND_MAX);
    return {cos(theta)*cos(phi), sin(theta)*cos(phi), sin(phi)};
}

// Output vector to a stream (like cout).
template<int N>
std::ostream& operator<<(std::ostream &o, const FVector<N> &x) {
    if (N == 0) return o << "()";
    o << "(" << x(0);
    for (int i=1; i<N; ++i) {
        o << ", " << x(i);
    }
    return o << ")";
}

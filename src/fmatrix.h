#pragma once
#include "fvector.h"
#include <algorithm>

// Fixed size matrix class: rows and cols must be known at compile time.
template<int I, int J>
class FMatrix {
public:
    double operator()(int i, int j) const { return _x[i+j*I]; } 
    double& operator()(int i, int j) { return _x[i+j*I]; } 

    double* begin() { return _x; }
    const double* begin() const { return _x; }
    double* end() { return _x+I*J; }
    const double* end() const { return _x+I*J; }

    constexpr int rows() const { return I; }
    constexpr int cols() const { return J; }

private:
    double _x[I*J];
};

//! Performs matrix-vector product.
template<int I, int J>
FVector<I> operator*(const FMatrix<I,J> &A, const FVector<J> &x) {
    FVector<I> y;
    for (int i=0; i<I; ++i) {
        y(i) = A(i,0) * x(0);
        for (int j=1; j<J; ++j) {
            y(i) += A(i,j)*x(j);
        }
    }
    return y;
}

//! Performs vector-matrix product.
template<int I, int J>
FVector<I> operator*(const FVector<J> &x, const FMatrix<J,I> &A) {
    FVector<I> y;
    for (int i=0; i<I; ++i) {
        y(i) = A(0,i) * x(0);
        for (int j=1; j<J; ++j) {
            y(i) += A(j,i)*x(j);
        }
    }
    return y;
}

template<int I, int J>
FMatrix<I,J> operator*(double x, const FMatrix<I,J> &A) {
    FMatrix<I,J> B;
    for (int i=0; i<I; i++) {
        for (int j=0; j<J; j++) {
            B(i,j) = x * A(i,j);
        }
    }
    return B;
}

template<int I, int J>
FMatrix<I,J> operator*(const FMatrix<I,J> &A, double x) {
    FMatrix<I,J> B;
    for (int i=0; i<I; i++) {
        for (int j=0; j<J; j++) {
            B(i,j) = x * A(i,j);
        }
    }
    return B;
}

template<int I, int J>
double scalar_product(const FMatrix<I,J> &A, const FMatrix<I,J> &B) {
    double prod = 0.0;
    for (int i=0; i<I; ++i) {
        for (int j=0; j<J; ++j) {
            prod += A(i,j)*B(i,j);
        }
    }
    return prod;
}

template<int I, int J>
FMatrix<I,J>& operator+=(FMatrix<I,J> &x, const FMatrix<I,J> &y) {
    for (int i=0; i<I; i++) {
        for (int j=0; j<J; j++) {
            x(i,j) += y(i,j);
        }
    }
    return x; 
}

template<int I, int J>
FMatrix<I,J>& operator-=(FMatrix<I,J> &x, const FMatrix<I,J> &y) {
    for (int i=0; i<I; i++) {
        for (int j=0; j<J; j++) {
            x(i,j) -= y(i,j);
        }
    }
    return x; 
}
template<int I, int J>
FMatrix<I,J> operator+(const FMatrix<I,J> &x, const FMatrix<I,J> &y) {
    FMatrix<I,J> z(x);
    z += y;
    return z; 
}
template<int I, int J>
FMatrix<I,J> operator-(const FMatrix<I,J> &x, const FMatrix<I,J> &y) {
    FMatrix<I,J> z(x);
    z -= y;
    return z; 
}

// Returns vectors outer/tensor product.
template <int M, int N>
FMatrix<M,N> outer(const FVector<M> &x, const FVector<N> &y) {
    FMatrix<M, N> z;
    for (int i=0; i<M; ++i) {
        for (int j=0; j<N; ++j) {
            z(i,j) = x(i)*y(j);
        }
    }
    return z; 
}
 
// Matrix-matrix product.
template <int I, int J, int K>
FMatrix<I,J> operator*(const FMatrix<I,K> &A, const FMatrix<K,J> &B) {
    FMatrix<I,J> C;
    for (int i=0; i<I; ++i) {
        for (int j=0; j<J; ++j) {
            C(i,j) = A(i,0)*B(0,j);
            for (int k=1; k<K; ++k) {
                C(i,j) += A(i,k)*B(k,j);
            }
        }
    }
    return C;
}

template <int N>
struct LUDecomposition {
    FMatrix<N,N> LU;
    int pivot[N];
};

// Performs an LU decomposition of the Matrix A.
template<int N> LUDecomposition<N> lu(FMatrix<N,N> A) {
    LUDecomposition<N> lud;
    int *pivot = lud.pivot;
    for (int i=0; i<N; ++i) pivot[i] = i;
    for (int j=0; j<N; ++j) {
        for (int i=0; i<j+1; ++i) 
            for (int k=0; k<i; ++k) 
                A(i,j) -= A(i,k)*A(k,j);
        int n = j;
        double col_max = fabs(A(j,j));
        for (int i=j+1; i<N; ++i) {
            for (int k=0; k<j; ++k) {
                A(i,j) -= A(i,k)*A(k,j);
            }
            if (fabs(A(i,j) > col_max)) {
                col_max = fabs(A(i,j));
                n = i;
            }
        }
        if (n != j) {
            for (int i=0; i<N; ++i) std::swap(A(n,i), A(j,i));
            std::swap(pivot[j], pivot[n]);
        }
        if (A(j,j) == 0.0) {
            std::cerr << "FMatrix is singular.\n";
        }
        for (int i=j+1; i<N; ++i) A(i,j) /= A(j,j);
    }
    lud.LU = std::move(A);
    return lud;
}

template<int N>
FVector<N> back_solve(const LUDecomposition<N> &lud, const FVector<N> &b) {
    FVector<N> x;
    for (int i=0; i<N; ++i) {
        x(i) = b(lud.pivot[i]);
        for (int j=0; j<i; ++j)   x(i) -= lud.LU(i,j)*x(j);
    }
    for (int i=N-1; i>=0; --i) {
        for (int j=i+1; j<N; ++j) x(i) -= lud.LU(i,j)*x(j);
        x(i) /= lud.LU(i,i);
    }
    return x;
}

inline double det(const FMatrix<2,2> &A) {
    return A(0,0)*A(1,1)-A(0,1)*A(1,0);
}
inline double det(const FMatrix<3,3> &A) {
    return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1)) 
           + A(0,1)*(A(1,2)*A(2,0)-A(1,0)*A(2,2))
           + A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
}
inline FMatrix<2,2> inv(const FMatrix<2,2> &A) {
    auto Z = 1.0/det(A);
    FMatrix<2,2> Ai;
    Ai(0,0) = Z*A(1,1);  Ai(0,1) =-Z*A(0,1);
    Ai(1,0) =-Z*A(1,0);  Ai(1,1) = Z*A(0,0);
    return Ai;
}
inline FMatrix<3,3> inv(const FMatrix<3,3> &A) {
    FMatrix<3,3> Ai;
    double invDetA = 1.0 / det(A);
    Ai(0,0) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))*invDetA;
    Ai(0,1) = (A(0,2)*A(2,1)-A(0,1)*A(2,2))*invDetA;
    Ai(0,2) = (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invDetA;
    Ai(1,0) = (A(1,2)*A(2,0)-A(1,0)*A(2,2))*invDetA;
    Ai(1,1) = (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invDetA;
    Ai(1,2) = (A(0,2)*A(1,0)-A(0,0)*A(1,2))*invDetA;
    Ai(2,0) = (A(1,0)*A(2,1)-A(1,1)*A(2,0))*invDetA;
    Ai(2,1) = (A(0,1)*A(2,0)-A(0,0)*A(2,1))*invDetA;
    Ai(2,2) = (A(0,0)*A(1,1)-A(0,1)*A(1,0))*invDetA;
    return Ai;
}

template<int N>
inline FMatrix<N,N> inv(const FMatrix<N,N> &M) {
    constexpr int N1 = N/2;
    constexpr int N2 = N-N1;

    // Using blockwise inversion, [A,B;C,D]-1.
    FMatrix<N1,N2> A,B,C,D;
    for (int i=0; i<N1; ++i) {
        for (int j=0; j<N2; ++j) {
            A(i,j) = M(i,j);
            B(i,j) = M(i,j+N2);
            C(i,j) = M(i+N1,j);
            D(i,j) = M(i+N1,j+N2);
        }
    }

    auto T = D-C*inv(A)*B;
    auto A_prime = inv(A) + inv(A)*B*inv(T)*C*inv(A);
    auto B_prime = -1.0*inv(A)*B*inv(T);
    auto C_prime = -1.0*inv(T)*C*inv(A);
    auto D_prime = inv(T);

    FMatrix<N,N> M_inv;
    for (int i=0; i<N1; ++i) {
        for (int j=0; j<N2; ++j) {
            M_inv(i, j)     = A_prime(i,j);
            M_inv(i, j+N2)   = B_prime(i,j);
            M_inv(i+N1, j)   = C_prime(i,j);
            M_inv(i+N1, j+N2) = D_prime(i,j); 
        }
    }
     
    return M_inv;
}
template<int N>
inline FMatrix<N,N> inv_lud(const FMatrix<N,N> &A) {
    auto lud = lu(A);
    // Invert each col of M.
    FMatrix<N,N> Ai;
    for (int j=0; j<N; ++j) {
        auto y = zeros<N>();
        y(j) = 1.0;
        y = back_solve(lud, y);
      for (int i=0; i<N; ++i) Ai(i,j) = y(i);
    }
    return Ai;
}
template <int I, int J>
inline FMatrix<J,I> transpose(const FMatrix<I,J> &A) {
    FMatrix<J,I> At;
    for (int i=0; i<I; ++i) {
        for (int j=0; j<J; ++j) {
            At(j,i) = A(i,j);
        }
    }
    return At;
}

template <int M, int N>
FMatrix<M,N> zeros() {
    FMatrix<M,N> A;
    for (auto &a: A) a = 0.0;
    return A;
}

template<int N>
FMatrix<N,N> eye() {
    auto I = zeros<N,N>();
    for (int i=0; i<N; ++i) I(i,i) = 1.0;
    return I;
}

// Combine two matrices into one, K=2*J. 
template<int I, int J>
inline FMatrix<I, 2*J> combine_cols(const FMatrix<I,J> &A, const FMatrix<I,J> &B) {
    FMatrix<I, 2*J> AB;
    for (int i=0; i<I; ++i) {
        for (int j=0; j<I; ++j) {
            AB(i,j)   = A(i,j);
            AB(i,j+J) = B(i,j);
        }
    }
    return AB;
}

template<int I, int J>
double norm(const FMatrix<I,J> &A) {
    double n = 0.0;
    for (auto &a: A) n += a*a;                
    return n;
}

// Output matrix to a stream (like cout).
template<int M, int N>
std::ostream& operator<<(std::ostream &o, const FMatrix<M,N> &A) {
    if (M*N == 0) return o << "[]";
    for (int i=0; i<M; ++i) {
        if (i>0) o << "\n";
        for (int j=0; j<N; ++j) {
            o << (j>0 ? "\t" : "") << A(i,j);
        }
    }
    return o;
}


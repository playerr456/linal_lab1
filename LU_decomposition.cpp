#include <cmath>
#include <stdexcept>

#include "LU_decomposition.h"

constexpr double EPS = 1e-12;

Vector forward_substitution(const Matrix& L, const Vector& B) {
    int n = static_cast<int>(L.size());
    Vector ans = create_vector(n);

    for (int i = 0; i < n; i ++) {
        double sum = B[i];

        for (int j = 0; j < i; j ++)
            sum -= L[i][j] * ans[j];

        if (std::abs(L[i][i]) < EPS)
            throw std::runtime_error("Matrix is singular or near-singular");

        ans[i] = sum / L[i][i];
    }
    return ans;
}

Vector back_substitution(const Matrix& U, const Vector& B) {
    int n = static_cast<int>(U.size());
    Vector ans = create_vector(n);

    for (int i = n - 1; i >= 0; i --) {
        double sum = B[i];

        for (int j = i + 1; j < n; j ++)
            sum -= U[i][j] * ans[j];

        if (std::abs(U[i][i]) < EPS)
            throw std::runtime_error("Matrix is singular or near-singular");

        ans[i] = sum / U[i][i];
    }
    return ans;
}

Vector solve_lu(const Matrix& L, const Matrix& U, const Vector& B) {
    Vector Y = forward_substitution(L, B);
    return back_substitution(U, Y);
}

std::pair<Matrix, Matrix> LU_decomposition(const Matrix& A) {
    int n = static_cast<int>(A.size());

    Matrix L = create_matrix(n);
    Matrix U = create_matrix(n);

    for (int i = 0; i < n; i ++) {
        for (int k = i; k < n; k ++) {
            double sum = 0;

            for (int j = 0; j < i; j ++)
                sum += L[i][j] * U[j][k];

            U[i][k] = A[i][k] - sum;
        }

        if (std::abs(U[i][i]) < EPS)
            throw std::runtime_error("Matrix is singular or near-singular");

        L[i][i] = 1;

        for (int k = i + 1; k < n; k ++) {
            double sum = 0;

            for (int j = 0; j < i; j ++)
                sum += L[k][j] * U[j][i];

            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }

    return {L, U};
}

Vector LU_method(const Matrix& A, const Vector& B) {
    if (A.size() != B.size ())
        throw std::runtime_error("The size of the matrix is not equal to the size of the vector column");

    auto [L, U] = LU_decomposition(A);
    return solve_lu(L, U, B);
}

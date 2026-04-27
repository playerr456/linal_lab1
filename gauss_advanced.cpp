#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "matrix.h"

constexpr double EPS = 1e-12;

static std::pair<Matrix, Vector> forward_elimination(Matrix A, Vector B) {
    int n = static_cast<int>(A.size());

    for (int i = 0; i < n - 1; i ++) {
        // advanced pivot selection
        double max_val = std::abs(A[i][i]);
        int max_row = i;
        
        for (int k = i + 1; k < n; k ++) {
            if (std::abs(A[k][i]) > max_val) {
                max_val = std::abs(A[k][i]);
                max_row = k;
            }
        }

        if (max_val < EPS)
            throw std::runtime_error("Matrix is singular or near-singular");

        std::swap(A[i], A[max_row]);
        std::swap(B[i], B[max_row]);

        for (int j = i + 1; j < n; j ++) { 
            double ratio = A[j][i] / A[i][i];

            for (int k = i; k < n; k ++) 
                A[j][k] -= ratio * A[i][k];
            B[j] -= ratio * B[i];
        }
    }
    return {A, B};
}

static Vector back_substitution(const Matrix& A, const Vector& B) {
    int n = static_cast<int>(A.size());
    Vector ans = create_vector(n);

    for (int i = n - 1; i >= 0; i --) {
        double sum = B[i];

        for (int j = i + 1; j < n; j ++) 
            sum -= A[i][j] * ans[j];

        if (std::abs(A[i][i]) < EPS)
            throw std::runtime_error("Matrix is singular or near-singular");

        ans[i] = sum / A[i][i];   
    }
    return ans;
}

Vector gauss_advanced(const Matrix& A, const Vector& B) {
    if (A.size() != B.size ())
        throw std::runtime_error("The size of the matrix is not equal to the size of the vector column");
    
    auto [C, D] = forward_elimination(A, B);
    return back_substitution(C, D);
}

#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "LU_decomposition.h"
#include "gauss_advanced.h"
#include "gauss_default.h"
#include "matrix.h"
#include "tests.h"

using Clock = std::chrono::high_resolution_clock;

constexpr double MIN_RANDOM_VALUE = -1.0;
constexpr double MAX_RANDOM_VALUE = 1.0;
constexpr int MAX_ATTEMPTS = 100;

Matrix create_good_matrix(int n) {
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++ attempt) {
        Matrix A = random_matrix(n, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        try {
            LU_decomposition(A);
            return A;
        } 
        catch (...) {}
    }

    throw std::runtime_error("Could not generate matrix");
}

Vector matrix_by_vector(const Matrix& A, const Vector& x) {
    int n = static_cast<int>(A.size());
    Vector result = create_vector(n);

    for (int i = 0; i < n; ++ i) {
        double sum = 0.0;

        for (int j = 0; j < n; ++ j)
            sum += A[i][j] * x[j];

        result[i] = sum;
    }

    return result;
}

double vector_norm(const Vector& x) {
    double sum = 0.0;

    for (double value : x)
        sum += value * value;

    return std::sqrt(sum);
}

double calc_res(const Matrix& A, const Vector& x, const Vector& b) {
    Vector diff = matrix_by_vector(A, x);

    for (int i = 0; i < static_cast<int>(diff.size()); ++ i)
        diff[i] -= b[i];

    return vector_norm(diff);
}

double rel_error(const Vector& x, const Vector& exact_x) {
    Vector diff = x;

    for (int i = 0; i < static_cast<int>(diff.size()); ++ i)
        diff[i] -= exact_x[i];

    return vector_norm(diff) / vector_norm(exact_x);
}

double elapsed_time(const Clock::time_point& start, const Clock::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

void first_test() {
    std::cout << "\n============= FIRST TEST =============\n";

    std::vector<int> sizes = {100, 200, 500, 1000};

    for (int n : sizes) {
        Matrix A = create_good_matrix(n);
        Vector b = random_vector(n, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE);

        auto start = Clock::now();
        try {
            gauss_advanced(A, b);
        } 
        catch (...) {}

        auto end = Clock::now();
        double time_with_pivot = elapsed_time(start, end);

        start = Clock::now();
        try {
            gauss_default(A, b);
        } 
        catch (...) {}

        end = Clock::now();
        double time_without_pivot = elapsed_time(start, end);

        Matrix L;
        Matrix U;

        start = Clock::now();
        try {
            std::tie(L, U) = LU_decomposition(A);
        } 
        catch (...) {}

        end = Clock::now();
        double time_lu_decomposition = elapsed_time(start, end);

        start = Clock::now();
        try {
            solve_lu(L, U, b);
        } 
        catch (...) {}

        end = Clock::now();
        double time_lu_solve = elapsed_time(start, end);

        double lu_total = time_lu_decomposition + time_lu_solve;

        std::cout << "n = " << n << "\n";
        std::cout << "Gauss with pivot: " << time_with_pivot << "\n";
        std::cout << "Gauss without pivot: " << time_without_pivot << "\n";
        std::cout << "LU decomposition: " << time_lu_decomposition << "\n";
        std::cout << "LU solve: " << time_lu_solve << "\n";
        std::cout << "LU total: " << lu_total << "\n";
        std::cout << "--------------------------\n";
    }
}

void second_test() {
    std::cout << "\n============= SECOND TEST =============\n";

    int n = 500;
    std::vector<int> sizes = {1, 10, 100};

    Matrix A = create_good_matrix(n);

    for (int k : sizes) {
        std::vector<Vector> b_vectors;

        for (int i = 0; i < k; i++)
            b_vectors.push_back(random_vector(n, MIN_RANDOM_VALUE, MAX_RANDOM_VALUE));

        auto start = Clock::now();

        for (auto& b : b_vectors) {
            try {
                gauss_advanced(A, b);
            } 
            catch (...) {}
        }

        auto end = Clock::now();
        double time_with_pivot = elapsed_time(start, end);

        start = Clock::now();

        try {
            Matrix L;
            Matrix U;

            std::tie(L, U) = LU_decomposition(A);

            for (auto& b : b_vectors)
                solve_lu(L, U, b);
        } 
        catch (...) {}

        end = Clock::now();
        double time_lu = elapsed_time(start, end);

        std::cout << "k = " << k << "\n";
        std::cout << "Gauss with pivot: " << time_with_pivot << "\n";
        std::cout << "LU total: " << time_lu << "\n";
        std::cout << "--------------------------\n";
    }
}

void third_test() {
    std::cout << "\n============= THIRD TEST =============\n";

    std::vector<int> sizes = {5, 10, 15};

    for (int n : sizes) {
        Matrix G = hilbert_matrix(n);
        Vector x_exact = create_vector(n, 1.0);
        Vector b = matrix_by_vector(G, x_exact);

        Vector x_gauss;
        Vector x_with_pivot;
        Vector x_lu;

        bool gauss = true;
        try {
            x_gauss = gauss_default(G, b);
        } 
        catch (...) {
            gauss = false;
        }

        bool with_pivot = true;
        try {
            x_with_pivot = gauss_advanced(G, b);
        } 
        catch (...) {
            with_pivot = false;
        }

        bool lu = true;
        try {
            Matrix L;
            Matrix U;
            std::tie(L, U) = LU_decomposition(G);
            x_lu = solve_lu(L, U, b);
        } 
        catch (...) {
            lu = false;
        }

        std::cout << "n = " << n << "\n";

        if (gauss) {
            std::cout << "Gauss error: " << rel_error(x_gauss, x_exact) << "\n";
            std::cout << "Gauss residual: " << calc_res(G, x_gauss, b) << "\n";
        } 
        else {
            std::cout << "Gauss failed\n";
        }

        if (with_pivot) {
            std::cout << "Gauss with pivot error: " << rel_error(x_with_pivot, x_exact) << "\n";
            std::cout << "Gauss with pivot residual: " << calc_res(G, x_with_pivot, b) << "\n";
        } 
        else {
            std::cout << "Gauss with pivot failed\n";
        }

        if (lu) {
            std::cout << "LU error: " << rel_error(x_lu, x_exact) << "\n";
            std::cout << "LU residual: " << calc_res(G, x_lu, b) << "\n";
        } 
        else {
            std::cout << "LU failed\n";
        }

        std::cout << "--------------------------\n";
    }
}

int main() {
    try {
        first_test();
        second_test();
        third_test();
    } 
    catch (const std::exception& error) {
        std::cout << "Fatal error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
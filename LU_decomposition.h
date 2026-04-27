#pragma once

#include <utility>

#include "matrix.h"

Vector forward_substitution(const Matrix& L, const Vector& B);
Vector back_substitution(const Matrix& U, const Vector& B);
Vector solve_lu(const Matrix& L, const Matrix& U, const Vector& B);
std::pair<Matrix, Matrix> LU_decomposition(const Matrix& A);
Vector LU_method(const Matrix& A, const Vector& B);

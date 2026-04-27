#pragma once

#include <initializer_list>
#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix create_matrix(int size, double value = 0);
Matrix create_matrix(std::initializer_list<std::initializer_list<double>> values);
Matrix random_matrix(int size, double min_value = 0, double max_value = 1);
Matrix hilbert_matrix(int n);

Vector create_vector(int size, double value = 0);
Vector create_vector(std::initializer_list<double> values);
Vector random_vector(int size, double min_value = 0, double max_value = 1);

#include <cstdlib>
#include <stdexcept>

#include "matrix.h"

namespace {
    constexpr unsigned int RANDOM_SEED = 123;

    void validate_matrix_size(int size) {
        if (size <= 0)
            throw std::runtime_error("Matrix size must be positive");
    }

    void validate_vector_size(int size) {
        if (size <= 0)
            throw std::runtime_error("Vector size must be positive");
    }

    void validate_random_range(double min_value, double max_value) {
        if (min_value >= max_value)
            throw std::runtime_error("Random range is invalid");
    }
    /*
    void ensure_random_seed() {
        static const bool seeded = []() {
            std::srand(RANDOM_SEED);
            return true;
        }();
        (void)seeded;
    }
    */
    void ensure_random_seed() {
        static bool seeded = false;

        if (!seeded) {
            std::srand(RANDOM_SEED);
            seeded = true;
        }
    }

    int get_random_seed() {
        return RANDOM_SEED;
    }

    double random_value(double min_value, double max_value) {
        ensure_random_seed();
        double ratio = static_cast<double>(std::rand()) / RAND_MAX;
        return min_value + (max_value - min_value) * ratio;
    }
}

Matrix create_matrix(int size, double value) {
    validate_matrix_size(size);
    return Matrix(size, Vector(size, value));
}

Matrix create_matrix(std::initializer_list<std::initializer_list<double>> values) {
    if (values.size() == 0)
        throw std::runtime_error("Matrix must not be empty");

    Matrix matrix;
    int size = static_cast<int>(values.size());
    int cols = static_cast<int>(values.begin()->size());

    if (cols == 0)
        throw std::runtime_error("Matrix must not have empty rows");

    if (cols != size)
        throw std::runtime_error("Matrix must be square");

    for (const auto& row : values) {
        if (row.size() != cols)
            throw std::runtime_error("All matrix rows must have the same size");

        matrix.push_back(row);
    }
    return matrix;
}

Matrix random_matrix(int size, double min_value, double max_value) {
    validate_matrix_size(size);
    validate_random_range(min_value, max_value);

    Matrix matrix = create_matrix(size);

    for (int i = 0; i < size; i ++) {
        for (int j = 0; j < size; j ++)
            matrix[i][j] = random_value(min_value, max_value);
    }
    return matrix;
}

Matrix hilbert_matrix(int n) {
    validate_matrix_size(n);

    Matrix matrix = create_matrix(n);

    for (int i = 0; i < n; i ++) {
        for (int j = 0; j < n; j ++)
            matrix[i][j] = 1 / static_cast<double>(i + j + 1);
    }
    return matrix;
}

Vector create_vector(int size, double value) {
    validate_vector_size(size);
    return Vector(size, value);
}

Vector create_vector(std::initializer_list<double> values) {
    if (values.size() == 0)
        throw std::runtime_error("Vector must not be empty");

    return Vector(values);
}

Vector random_vector(int size, double min_value, double max_value) {
    validate_vector_size(size);
    validate_random_range(min_value, max_value);

    Vector vector = create_vector(size);

    for (int i = 0; i < size; i ++)
        vector[i] = random_value(min_value, max_value);

    return vector;
}

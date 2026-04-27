# Linear Algebra Lab 1

C++ project for comparing several methods for solving systems of linear equations:

- Gaussian elimination without pivoting
- Gaussian elimination with partial pivoting
- LU decomposition

The repository also includes helpers for matrix/vector generation and three experiment scenarios for performance and numerical-stability analysis.

## Files

- `matrix.h`, `matrix.cpp` - matrix and vector creation helpers, random generators, Hilbert matrix generation
- `gauss_default.h`, `gauss_default.cpp` - Gaussian elimination without pivoting
- `gauss_advanced.h`, `gauss_advanced.cpp` - Gaussian elimination with partial pivoting
- `LU_decomposition.h`, `LU_decomposition.cpp` - LU decomposition and forward/back substitution
- `tests.h`, `tests.cpp` - experiment runner with three test scenarios

## Mathematical Background

The code uses zero-based indices internally, but the formulas below are written in the standard one-based mathematical form.

### 1. Gaussian Elimination Without Pivoting

For elimination step $k = 1, 2, \dots, n - 1$, the multiplier is

$$
m_{ik} = \frac{A_{ik}}{A_{kk}}, \qquad i = k + 1, \dots, n.
$$

Then the matrix row and right-hand side are updated as

$$
A_{ij} \leftarrow A_{ij} - m_{ik} A_{kj}, \qquad j = k, k + 1, \dots, n,
$$

$$
b_i \leftarrow b_i - m_{ik} b_k.
$$

After forward elimination, the upper-triangular system is solved by back substitution:

$$
x_i = \frac{b_i - \sum_{j = i + 1}^{n} U_{ij} x_j}{U_{ii}}, \qquad i = n, n - 1, \dots, 1.
$$

### 2. Gaussian Elimination With Partial Pivoting

At each step, the pivot row is chosen by

$$
p = \arg \max_{r = k, \dots, n} |A_{rk}|.
$$

Then rows $k$ and $p$ are swapped in both $A$ and $b$, after which the same elimination formulas are applied:

$$
m_{ik} = \frac{A_{ik}}{A_{kk}}, \qquad
A_{ij} \leftarrow A_{ij} - m_{ik} A_{kj}, \qquad
b_i \leftarrow b_i - m_{ik} b_k.
$$

### 3. LU Decomposition

The matrix is factorized into

$$
A = LU,
$$

where $L$ is lower triangular and $U$ is upper triangular.

The coefficients of $U$ are computed by

$$
U_{ik} = A_{ik} - \sum_{j = 1}^{i - 1} L_{ij} U_{jk}, \qquad k \ge i.
$$

The coefficients of $L$ are computed by

$$
L_{ki} = \frac{A_{ki} - \sum_{j = 1}^{i - 1} L_{kj} U_{ji}}{U_{ii}}, \qquad k > i.
$$

### 4. Forward and Back Substitution

For the lower-triangular system $Ly = b$:

$$
y_i = \frac{b_i - \sum_{j = 1}^{i - 1} L_{ij} y_j}{L_{ii}}, \qquad i = 1, \dots, n.
$$

For the upper-triangular system $Ux = y$:

$$
x_i = \frac{y_i - \sum_{j = i + 1}^{n} U_{ij} x_j}{U_{ii}}, \qquad i = n, n - 1, \dots, 1.
$$

### 5. Auxiliary Formulas Used in the Experiments

Euclidean norm of a vector:

$$
\|v\|_2 = \sqrt{\sum_{i = 1}^{n} v_i^2}.
$$

Componentwise vector difference:

$$
(v_1 - v_2)_i = v_{1i} - v_{2i}, \qquad i = 1, 2, \dots, n.
$$

Relative error of the computed solution $\tilde{x}$ with respect to the exact solution $x$:

$$
\varepsilon_{\mathrm{rel}} =
\frac{\|\tilde{x} - x\|_2}{\|x\|_2}
=
\sqrt{
\frac{\sum_{i = 1}^{n} (\tilde{x}_i - x_i)^2}
{\sum_{i = 1}^{n} x_i^2}
}.
$$

Residual norm:

$$
\rho = \|A \tilde{x} - b\|_2
=
\sqrt{
\sum_{i = 1}^{n}
\left(
\sum_{j = 1}^{n} a_{ij} \tilde{x}_j - b_i
\right)^2
}.
$$

## Test Parameters

The tables below contain the exact test inputs used in `tests.cpp`.

### Test 1: Runtime Comparison for Different Matrix Sizes

Random matrices are generated in the range $[-1, 1]$. A matrix is retried up to `MAX_ATTEMPTS = 100` times until LU decomposition succeeds.

| Matrix size `n` | Right-hand side `b` | Methods compared | Measured values |
| --- | --- | --- | --- |
| 100 | random vector in `[-1, 1]` | Gauss with pivot, Gauss without pivot, LU | execution time |
| 200 | random vector in `[-1, 1]` | Gauss with pivot, Gauss without pivot, LU | execution time |
| 500 | random vector in `[-1, 1]` | Gauss with pivot, Gauss without pivot, LU | execution time |
| 1000 | random vector in `[-1, 1]` | Gauss with pivot, Gauss without pivot, LU | execution time |

### Test 2: Repeated Solves for the Same Matrix

For this experiment, one random matrix of size $n = 500$ is reused, and the number of different right-hand sides varies.

| Matrix size `n` | Number of RHS vectors `k` | RHS generation | Methods compared | Measured values |
| --- | --- | --- | --- | --- |
| 500 | 1 | random vectors in `[-1, 1]` | Gauss with pivot, LU | total execution time |
| 500 | 10 | random vectors in `[-1, 1]` | Gauss with pivot, LU | total execution time |
| 500 | 100 | random vectors in `[-1, 1]` | Gauss with pivot, LU | total execution time |

### Test 3: Numerical Stability on Hilbert Matrices

For each Hilbert matrix $G$, the exact solution is chosen as

$$
x_{\mathrm{exact}} = (1, 1, \dots, 1)^T,
$$

and the right-hand side is built as

$$
b = G x_{\mathrm{exact}}.
$$

| Matrix type | Matrix size `n` | Exact solution | Methods compared | Measured values |
| --- | --- | --- | --- | --- |
| Hilbert | 5 | all ones | Gauss without pivot, Gauss with pivot, LU | relative error, residual |
| Hilbert | 10 | all ones | Gauss without pivot, Gauss with pivot, LU | relative error, residual |
| Hilbert | 15 | all ones | Gauss without pivot, Gauss with pivot, LU | relative error, residual |

## Build

The project uses CMake and requires a C++17 compiler.

```powershell
cmake -S . -B build
cmake --build build
```

## Run

After building, run:

```powershell
.\build\Debug\linal_lab1.exe
```

On single-config generators such as MinGW Makefiles or Ninja, the executable will usually be located at:

```powershell
.\build\linal_lab1.exe
```

## Notes

- Random matrix and vector generation uses a fixed seed for reproducibility.
- Singular or near-singular matrices throw exceptions.
- `tests.cpp` acts as the main executable for the project rather than a unit-test suite.

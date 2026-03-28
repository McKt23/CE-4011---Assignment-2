import math

# =================================================================
# 1. LIBRARY CODES (Class Definitions)
# =================================================================

class Matrix:
    """ Class for standard dense matrix operations (Local Matrices). """

    def __init__(self, rows, cols, data=None):
        self.rows = rows
        self.cols = cols
        if data:
            self.data = data
        else:
            self.data = [[0.0 for _ in range(cols)] for _ in range(rows)]

    def get(self, i, j):
        return self.data[i][j]

    def set(self, i, j, val):
        self.data[i][j] = val

    def transpose(self):
        """ Returns the transpose of the matrix. """
        res_data = [[0.0 for _ in range(self.rows)] for _ in range(self.cols)]
        for i in range(self.rows):
            for j in range(self.cols):
                res_data[j][i] = self.data[i][j]
        return Matrix(self.cols, self.rows, res_data)

    def matmul(self, other):
        """ Multiplies this matrix with another matrix. """
        if self.cols != other.rows:
            raise ValueError("Matrix dimensions are incompatible for multiplication.")

        res_data = [[0.0 for _ in range(other.cols)] for _ in range(self.rows)]
        for i in range(self.rows):
            for j in range(other.cols):
                total = 0.0
                for k in range(self.cols):
                    total += self.data[i][k] * other.data[k][j]
                res_data[i][j] = total
        return Matrix(self.rows, other.cols, res_data)

    def vec_mul(self, vec):
        """ Multiplies this matrix with a 1D vector. """
        if self.cols != len(vec):
            raise ValueError("Matrix column count must equal the vector element count.")

        res_vec = [0.0 for _ in range(self.rows)]
        for i in range(self.rows):
            total = 0.0
            for j in range(self.cols):
                total += self.data[i][j] * vec[j]
            res_vec[i] = total
        return res_vec

    def display(self):
        """ Prints the matrix in a readable format. """
        for row in self.data:
            print(["{:.4f}".format(val) for val in row])


class SymmetricBandedMatrix:
    """ Symmetric Banded Storage and Solver (Banded Cholesky) Class. """

    def __init__(self, n, m):
        self.n = n
        self.m = m
        # Compressed storage: Size is [n x m] instead of [n x n]
        self.data = [[0.0 for _ in range(m)] for _ in range(n)]

    def add_value(self, i, j, val):
        """ Adds a value to the global stiffness matrix enforcing symmetry. """
        # Enforce Symmetry: always store in the upper triangle
        if i > j:
            i, j = j, i

        diff = j - i
        # Check against bandwidth limit
        if diff < self.m:
            self.data[i][diff] += val
        else:
            raise ValueError(f"Error: Element exceeds the bandwidth limit! i={i}, j={j}")

    def solve(self, B):
        """ Solves the linear system [K]{x} = {B} using Banded Cholesky Factorization. """
        n = self.n
        m = self.m
        U = [[0.0 for _ in range(m)] for _ in range(n)]

        # 1. Cholesky Factorization (A = U^T * U)
        for i in range(n):
            for j in range(min(m, n - i)):
                col = i + j
                total = self.data[i][j]

                start_k = max(0, col - m + 1)
                for k in range(start_k, i):
                    total -= U[k][i - k] * U[k][col - k]

                if i == col:
                    if total <= 0.0:
                        raise ValueError(f"Matrix is not positive-definite! (Index: {i})")
                    U[i][0] = math.sqrt(total)
                else:
                    U[i][j] = total / U[i][0]

        # 2. Forward Substitution (U^T * y = B)
        y = [0.0 for _ in range(n)]
        for i in range(n):
            total = B[i]
            start_k = max(0, i - m + 1)
            for k in range(start_k, i):
                total -= U[k][i - k] * y[k]
            y[i] = total / U[i][0]

        # 3. Backward Substitution (U * x = y)
        x = [0.0 for _ in range(n)]
        for i in range(n - 1, -1, -1):
            total = y[i]
            for j in range(1, min(m, n - i)):
                total -= U[i][j] * x[i + j]
            x[i] = total / U[i][0]

        return x

    def display_band(self):
        """ Displays the compressed banded storage format. """
        for i, row in enumerate(self.data):
            print(f"Row {i} (Band):", ["{:.4f}".format(val) for val in row])


# =================================================================
# 2. TEST AND EXECUTION BLOCK (Results outputted to console)
# =================================================================

if __name__ == "__main__":
    print("==================================================")
    print("TEST 1: DENSE MATRIX OPERATIONS (Local Matrices)")
    print("==================================================")

    # Let's create two 2x2 matrices using our custom library (without using numpy!)
    R = Matrix(2, 2, [[0.8, 0.6], [-0.6, 0.8]])  # Example Rotation Matrix
    k = Matrix(2, 2, [[10.0, -10.0], [-10.0, 10.0]])  # Example Local Stiffness Matrix

    print("Rotation Matrix [R]:")
    R.display()

    print("\nTranspose of [R] -> [R]^T:")
    R_T = R.transpose()
    R_T.display()

    print("\nMatrix Multiplication: [R]^T * [k]:")
    RT_k = R_T.matmul(k)
    RT_k.display()

    print("\nVector Multiplication: [k] * {d'}:")
    d_local = [0.5, -0.5]  # Example local displacement vector
    f_local = k.vec_mul(d_local)
    print("Local forces {f'} =", ["{:.4f}".format(v) for v in f_local])

    print("\n==================================================")
    print("TEST 2: SYMMETRIC BANDED MATRIX SOLVER (Global K)")
    print("==================================================")

    # Imagine a simple 3x3 truss system. Bandwidth m=2
    # [ 2.0  -1.0   0.0 ]
    # [-1.0   2.0  -1.0 ]
    # [ 0.0  -1.0   2.0 ]
    n_dof = 3
    band_width = 2

    K_global = SymmetricBandedMatrix(n_dof, band_width)

    # Adding matrix elements (Providing only the upper triangle and diagonal!)
    K_global.add_value(0, 0, 2.0)
    K_global.add_value(0, 1, -1.0)  # K[0,1]
    K_global.add_value(1, 1, 2.0)
    K_global.add_value(1, 2, -1.0)  # K[1,2]
    K_global.add_value(2, 2, 2.0)

    print("Banded Storage format of [K] (nxM size, heavily compressed!):")
    K_global.display_band()

    # External force vector {F}
    F_vector = [10.0, 20.0, 30.0]

    print("\nSolving [K]{D} = {F} using custom Banded Cholesky Solver...")
    print("Force Vector {F}:", F_vector)

    # Solve the system of equations!
    Displacements = K_global.solve(F_vector)

    print("\nCalculated Displacements {D}:")
    for i, d in enumerate(Displacements):
        print(f"  D[{i}] = {d:.4f}")

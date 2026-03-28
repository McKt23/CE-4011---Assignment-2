import math


class Matrix:
    """
    ===========================================================================
    Purpose: Handles standard dense matrix operations required for local member 
             formulations (e.g., local stiffness [k] and rotation [R] matrices).
    Assumptions: Matrix elements are floating-point numbers.
    Inputs: rows (int), cols (int), data (optional 2D list of floats).
    Outputs: Matrix object.
    ===========================================================================
    """

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
        """
        Purpose: Computes the transpose of the matrix.
        Inputs: None.
        Outputs: A new Matrix object representing the transpose.
        """
        res_data = [[0.0 for _ in range(self.rows)] for _ in range(self.cols)]
        for i in range(self.rows):
            for j in range(self.cols):
                res_data[j][i] = self.data[i][j]
        return Matrix(self.cols, self.rows, res_data)

    def matmul(self, other):
        """
        Purpose: Multiplies this matrix with another Matrix object.
        Inputs: other (Matrix object).
        Outputs: A new Matrix object representing the product.
        Assumptions: self.cols must equal other.rows.
        """
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
        """
        Purpose: Multiplies this matrix with a 1D vector.
        Inputs: vec (list of floats).
        Outputs: A new list of floats representing the product vector.
        Assumptions: self.cols must equal len(vec).
        """
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
    """
    ===========================================================================
    Purpose: Stores the Global Stiffness Matrix [K] using a memory-efficient
             banded scheme and solves the linear system using Cholesky factor.
    Assumptions: The matrix is symmetric and positive-definite.
    Inputs: n (total equations), m (semi-bandwidth).
    Outputs: SymmetricBandedMatrix object.
    ===========================================================================
    """

    def __init__(self, n, m):
        self.n = n
        self.m = m
        self.data = [[0.0 for _ in range(m)] for _ in range(n)]

    def add_value(self, i, j, val):
        """
        Purpose: Assembles a value into the global stiffness matrix enforcing symmetry.
        Inputs: i (row index), j (col index), val (float value to add).
        Outputs: None (modifies internal state).
        """
        if i > j:
            i, j = j, i

        diff = j - i
        if diff < self.m:
            self.data[i][diff] += val
        else:
            raise ValueError(f"Error: Element exceeds the bandwidth limit! i={i}, j={j}")

    def solve(self, B):
        """
        Purpose: Solves the linear system [K]{x} = {B} using Banded Cholesky.
        Inputs: B (list of floats representing the force vector).
        Outputs: x (list of floats representing the displacement vector).
        """
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
                        raise ValueError(f"Matrix is not positive-definite at index {i}")
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
# TEST AND EXECUTION BLOCK
# =================================================================
if __name__ == "__main__":
    print("==================================================")
    print("TEST 1: DENSE MATRIX OPERATIONS (Local Matrices)")
    print("==================================================")

    # Creating 2x2 matrices without numpy
    R = Matrix(2, 2, [[0.8, 0.6], [-0.6, 0.8]])
    k = Matrix(2, 2, [[10.0, -10.0], [-10.0, 10.0]])

    print("Rotation Matrix [R]:")
    R.display()

    print("\nTranspose of [R] -> [R]^T:")
    R_T = R.transpose()
    R_T.display()

    print("\nMatrix Multiplication: [R]^T * [k]:")
    RT_k = R_T.matmul(k)
    RT_k.display()

    print("\nVector Multiplication: [k] * {d'}:")
    d_local = [0.5, -0.5]
    f_local = k.vec_mul(d_local)
    print("Local forces {f'} =", ["{:.4f}".format(v) for v in f_local])

    print("\n==================================================")
    print("TEST 2: SYMMETRIC BANDED MATRIX SOLVER (Global K)")
    print("==================================================")

    # 3x3 simple system with bandwidth m=2
    n_dof = 3
    band_width = 2

    K_global = SymmetricBandedMatrix(n_dof, band_width)

    # Adding upper triangle elements
    K_global.add_value(0, 0, 2.0)
    K_global.add_value(0, 1, -1.0)
    K_global.add_value(1, 1, 2.0)
    K_global.add_value(1, 2, -1.0)
    K_global.add_value(2, 2, 2.0)

    print("Banded Storage format of [K] (Compressed):")
    K_global.display_band()

    # External force vector {F}
    F_vector = [10.0, 20.0, 30.0]

    print("\nSolving [K]{D} = {F} using custom Banded Cholesky Solver...")
    print("Force Vector {F}:", F_vector)

    # Solving the system
    Displacements = K_global.solve(F_vector)

    print("\nCalculated Displacements {D}:")
    for i, d in enumerate(Displacements):
        print(f"  D[{i}] = {d:.4f}")

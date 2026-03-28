import math
from Assignment2_Q2_MatrixLibrary import Matrix


def get_local_stiffness(E, A, I, L):
    """
    ===========================================================================
    Purpose: Generates the 6x6 local stiffness matrix [k] for a 2D frame element.
    Inputs: E (Elasticity Modulus), A (Area), I (Inertia), L (Length).
    Outputs: Matrix object representing [k_local].
    Units: Standard SI units (or consistent units matching the input).
    ===========================================================================
    """
    k = Matrix(6, 6)

    k1 = E * A / L
    k2 = 12 * E * I / (L ** 3)
    k3 = 6 * E * I / (L ** 2)
    k4 = 4 * E * I / L
    k5 = 2 * E * I / L

    # Diagonal and cross terms
    k.set(0, 0, k1);
    k.set(0, 3, -k1)
    k.set(3, 0, -k1);
    k.set(3, 3, k1)

    k.set(1, 1, k2);
    k.set(1, 4, -k2);
    k.set(1, 2, k3);
    k.set(1, 5, k3)
    k.set(4, 1, k2);
    k.set(4, 4, k2);
    k.set(4, 2, -k3);
    k.set(4, 5, -k3)

    k.set(2, 1, k3);
    k.set(2, 4, -k3);
    k.set(2, 2, k4);
    k.set(2, 5, k5)
    k.set(5, 1, k3);
    k.set(5, 4, -k3);
    k.set(5, 2, k5);
    k.set(5, 5, k4)

    # Enforce symmetry manually
    for i in range(6):
        for j in range(i + 1, 6):
            k.set(j, i, k.get(i, j))

    return k


def get_rotation_matrix(dx, dy, L):
    """
    ===========================================================================
    Purpose: Generates the 6x6 Rotation matrix [R] for coordinate transformation.
    Inputs: dx (X projection), dy (Y projection), L (Length).
    Outputs: Matrix object representing [R].
    ===========================================================================
    """
    R = Matrix(6, 6)
    c = dx / L
    s = dy / L

    R.set(0, 0, c);
    R.set(0, 1, s)
    R.set(1, 0, -s);
    R.set(1, 1, c)
    R.set(2, 2, 1.0)

    R.set(3, 3, c);
    R.set(3, 4, s)
    R.set(4, 3, -s);
    R.set(4, 4, c)
    R.set(5, 5, 1.0)

    return R


# =====================================================================
# TEST BLOCK: Sadece bu dosyayı test etmek için eklendi
# =====================================================================
if __name__ == "__main__":
    print("=========================================================")
    print("TESTING: structural_elements.py")
    print("=========================================================\n")

    # PDF'teki örnek yapının 1 Numaralı Elemanının özelliklerini test edelim
    # E = 200,000 MPa (2e8 kN/m2), A = 0.02 m2, I = 0.08 m4, L = 3.0 m
    E_test = 200000.0
    A_test = 0.02
    I_test = 0.08
    L_test = 3.0

    k_test = get_local_stiffness(E_test, A_test, I_test, L_test)

    print(f"Generated 6x6 Local Stiffness Matrix [k] for L={L_test}m :")
    k_test.display()

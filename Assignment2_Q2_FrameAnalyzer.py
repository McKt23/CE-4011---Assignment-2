from Assignment2_Q2_MatrixLibrary import SymmetricBandedMatrix
from Assignment2_Q2_StructuralElements import get_local_stiffness, get_rotation_matrix
import math


def assign_equation_numbers(S_matrix, num_nodes):
    """
    ===========================================================================
    Purpose: Maps physical DOFs to active equation numbers, ignoring restraints.
    Inputs: S_matrix (Boundary conditions list), num_nodes (Total nodes).
    Outputs: E_mat (2D list mapping Node->Equation), NumEq (Total active DOFs).
    ===========================================================================
    """
    E_mat = [[0, 0, 0] for _ in range(num_nodes)]

    # 1. Mark Restraints (1 = restrained)
    for support in S_matrix:
        node_idx = support[0] - 1
        E_mat[node_idx][0] = support[1]
        E_mat[node_idx][1] = support[2]
        E_mat[node_idx][2] = support[3]

    # 2. Assign consecutive numbers to free DOFs
    counter = 1
    for i in range(num_nodes):
        for j in range(3):
            if E_mat[i][j] == 0:  # Free DOF
                E_mat[i][j] = counter
                counter += 1
            else:  # Restrained DOF
                E_mat[i][j] = 0

    NumEq = counter - 1
    return E_mat, NumEq


def calculate_bandwidth(C_matrix, E_mat):
    """
    ===========================================================================
    Purpose: Computes the semi-bandwidth (m) for the Global Stiffness Matrix.
    Inputs: C_matrix (Connectivity), E_mat (Equation numbers).
    Outputs: max_band (integer representing the semi-bandwidth).
    ===========================================================================
    """
    max_band = 0
    for elem in C_matrix:
        start_idx = elem[0] - 1
        end_idx = elem[1] - 1
        dofs = E_mat[start_idx] + E_mat[end_idx]
        active_dofs = [d for d in dofs if d != 0]
        if active_dofs:
            local_band = max(active_dofs) - min(active_dofs) + 1
            if local_band > max_band:
                max_band = local_band
    return max_band


def assemble_global_matrices(NumEq, max_band, C, XY, M, E_mat, L_loads):
    """
    ===========================================================================
    Purpose: Assembles the Global Stiffness Matrix [K] and Force Vector {F}.
    Inputs: System geometry, connectivity, materials, equations, and loads.
    Outputs: K_global (SymmetricBandedMatrix), F_global (list), saved_k, saved_R.
    ===========================================================================
    """
    K_global = SymmetricBandedMatrix(NumEq, max_band)
    F_global = [0.0 for _ in range(NumEq)]

    saved_k_locs = []
    saved_R_mats = []

    # 1. Assemble Stiffness [K]
    for elem in C:
        start_idx, end_idx, mat_idx = elem[0] - 1, elem[1] - 1, elem[2] - 1
        A, I_val, E_mod = M[mat_idx]
        x1, y1 = XY[start_idx]
        x2, y2 = XY[end_idx]

        dx, dy = x2 - x1, y2 - y1
        Length = math.sqrt(dx ** 2 + dy ** 2)

        k_loc = get_local_stiffness(E_mod, A, I_val, Length)
        R_mat = get_rotation_matrix(dx, dy, Length)

        saved_k_locs.append(k_loc)
        saved_R_mats.append(R_mat)

        # k_glob = R^T * k_loc * R
        R_T = R_mat.transpose()
        temp = R_T.matmul(k_loc)
        k_glob = temp.matmul(R_mat)

        dof_map = E_mat[start_idx] + E_mat[end_idx]
        for p in range(6):
            for q in range(6):
                P, Q = dof_map[p], dof_map[q]
                if P != 0 and Q != 0 and P <= Q:
                    K_global.add_value(P - 1, Q - 1, k_glob.get(p, q))

    # 2. Assemble Forces {F}
    for load in L_loads:
        node_idx = load[0] - 1
        for q in range(3):
            Q = E_mat[node_idx][q]
            if Q != 0:
                F_global[Q - 1] += load[q + 1]

    return K_global, F_global, saved_k_locs, saved_R_mats


def calculate_member_forces(C, E_mat, Displacements, saved_k_locs, saved_R_mats):
    """
    ===========================================================================
    Purpose: Calculates internal forces (Axial, Shear, Moment) for each member.
    Inputs: C, E_mat, Displacements, saved_k_locs, saved_R_mats.
    Outputs: Prints the member forces directly to the console.
    ===========================================================================
    """
    print("\n--- Member End Forces (Local Coordinates) ---")
    for i, elem in enumerate(C):
        start_idx = elem[0] - 1
        end_idx = elem[1] - 1

        dof_map = E_mat[start_idx] + E_mat[end_idx]
        d_glob_array = [0.0] * 6

        for p in range(6):
            P = dof_map[p]
            if P != 0:
                d_glob_array[p] = Displacements[P - 1]

        R_mat = saved_R_mats[i]
        k_loc = saved_k_locs[i]

        # {d'} = [R] * {d_glob}
        d_loc = R_mat.vec_mul(d_glob_array)

        # {f'} = [k_loc] * {d'}
        f_loc = k_loc.vec_mul(d_loc)

        print(f"Element {i + 1}:")
        print(f"  Start Node (Axial, Shear, Moment): {f_loc[0]:9.2f}, {f_loc[1]:9.2f}, {f_loc[2]:9.2f}")
        print(f"  End Node   (Axial, Shear, Moment): {f_loc[3]:9.2f}, {f_loc[4]:9.2f}, {f_loc[5]:9.2f}")
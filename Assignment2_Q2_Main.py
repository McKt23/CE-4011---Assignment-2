import Assignment2_Q2_FrameAnalyzer as fa

"""
===========================================================================
Purpose: Main entry point. Defines input data from ODTUClass guidelines, 
         orchestrates the analysis modules, and prints output data formats.
===========================================================================
"""

if __name__ == "__main__":
    print("==========================================================")
    print("      2D FRAME ANALYSIS PROGRAM (MODULAR OOP DESIGN)      ")
    print("==========================================================")

    # ----------------------------------------------------------------
    # A. INPUT PHASE (Data matched with ODTUClass Sample Structure)
    # ----------------------------------------------------------------
    NumNode = 4

    # Nodal Coordinates [X, Y]
    XY = [[0.0, 0.0], [0.0, 3.0], [4.0, 3.0], [4.0, 0.0]]

    # Material Properties [Area, Inertia, Elasticity]
    M_props = [[0.02, 0.08, 200000.0],
               [0.01, 0.01, 200000.0]]

    # Element Connectivity [Start, End, Mat_ID]
    C_mat = [[1, 2, 1],
             [2, 3, 1],
             [4, 3, 1],
             [1, 3, 2]]

    # Boundary Conditions [Node_ID, Restr_X, Restr_Y, Restr_Z]
    S_mat = [[1, 1, 1, 0],
             [4, 0, 1, 0]]

    # Applied Loads [Node_ID, Force_X, Force_Y, Moment_Z]
    L_loads = [[2, 10.0, -10.0, 0.0],
               [3, 10.0, -10.0, 0.0]]

    # ----------------------------------------------------------------
    # EXECUTION PHASES
    # ----------------------------------------------------------------

    # B. Equation Numbering
    E_mat, NumEq = fa.assign_equation_numbers(S_mat, NumNode)
    print("\n--- Equation Numbering Matrix (E) ---")
    for row in E_mat: print(row)
    print(f"Total Active Equations (NumEq): {NumEq}")

    # Calculate Bandwidth
    m_band = fa.calculate_bandwidth(C_mat, E_mat)
    print(f"Calculated Semi-Bandwidth (m) : {m_band}")

    # C & D. Global Stiffness Matrix & Force Vector Assembly
    K_global, F_global, saved_k, saved_R = fa.assemble_global_matrices(
        NumEq, m_band, C_mat, XY, M_props, E_mat, L_loads
    )

    print("\n--- Force Vector {F} ---")
    print(["{:.2f}".format(f) for f in F_global])

    print("\n--- Global Stiffness Matrix [K] (Banded Format) ---")
    K_global.display_band()

    # E. Solution
    print("\nSolving System [K]{D} = {F} using Custom Banded Solver...")
    Displacements = K_global.solve(F_global)

    print("\n--- Global Displacements {D} ---")
    for i, d in enumerate(Displacements):
        print(f"  D[{i + 1}] = {d:.3e}")

    # F. Member End Forces
    fa.calculate_member_forces(C_mat, E_mat, Displacements, saved_k, saved_R)

    print("\n================= ANALYSIS COMPLETE ======================")
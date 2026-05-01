# 🏗️ Pure Python 2D Frame Analysis Engine (Zero Dependencies)

![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)
![Dependencies](https://img.shields.io/badge/Dependencies-Zero-success.svg)
![Solver](https://img.shields.io/badge/Solver-Banded%20Cholesky-green.svg)
![Methodology](https://img.shields.io/badge/Method-Direct%20Stiffness-orange.svg)

A modular, object-oriented 2D Frame Analysis software developed strictly from scratch in Python. The core objective of this project is to implement the **Finite Element Method (FEM)** and global matrix assembly processes without relying on any external numerical or data manipulation libraries (such as NumPy or Pandas).

## 🚀 Engineering Objective

Commercial structural software often acts as a "black box." This engine was built to demonstrate a fundamental, under-the-hood understanding of computational matrix mathematics, memory allocation (banded storage), and structural mechanics principles required to analyze 2D frame structures.

## ✨ Core Features & Software Architecture

### 🧮 1. Custom Matrix & Math Library
Built entirely from standard Python lists to handle all linear algebra operations.
* **`Matrix` Class:** Handles standard dense operations ($[R]^T [k] [R]$), vector multiplications, and transpositions.
* **`SymmetricBandedMatrix` Class:** A highly memory-efficient storage scheme that compresses the massive $N \times N$ Global Stiffness Matrix into an $N \times m$ array (where $m$ is the semi-bandwidth), drastically reducing RAM usage.
* **Custom Solver:** Implements a localized **Banded Cholesky Factorization** algorithm to solve the global linear system $\{F\} = [K]\{D\}$.

### 🧱 2. Object-Oriented Assembly (Direct Stiffness Method)
The structural analysis workflow is divided into highly decoupled, specialized modules:
* **Structural Elements Module:** Generates the $6 \times 6$ local stiffness matrix ($[k]$) and transformation/rotation matrices ($[T]$) based on material properties ($E, A, I$).
* **Frame Analyzer Module:** Automates equation numbering (masking restrained DOFs), calculates the structural semi-bandwidth dynamically, and assembles the global $[K]$ and $\{F\}$ matrices.

### 📊 3. Zero-Dependency Data Structure
To maintain a lightweight and purely algorithmic footprint, structural models (Nodal coordinates, Connectivity, Sections, Nodal Loads, and Boundary Conditions) are defined directly via 2D Python arrays (`List[List[float]]`) rather than external spreadsheets.

## 🛠️ Program Output

The execution script solves the 2D frame and outputs highly precise data to the console:
1. Equation Numbering Matrix and Active DOFs.
2. The compressed Banded Storage format of the Global Stiffness Matrix.
3. Nodal Displacements ($U_x, U_y, R_z$).
4. Member End Forces (Axial, Shear, Moment) in local coordinates.

## 💻 Installation & Usage

### Prerequisites
* Python 3.8 or higher.
* **No `pip install` required.** (Zero external dependencies).

### Running the Engine
1. Define your structural topology directly inside the execution script (`Assignment2_Q2_Main.py`).
2. Run the program:
   ```bash
   python Assignment2_Q2_Main.py

Proving that robust structural engineering tools can be built purely on fundamental logic and efficient algorithms.

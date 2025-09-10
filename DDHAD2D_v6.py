# -*- coding: utf-8 -*-
"""
DDHADtens 2D v1
J. Yvonnet September, 2025
Gustave Eiffel University
MSME group

Input : a collection of 3x3 elastic tensors (stored in .mat files)
Output: associated damage variables (harmonic analysis)
        - alpha vectors (size 6)
        - plots of the evolution of internal variables and 2D damage surfaces

The database is associated with the numerical examples in:
[1] J. Yvonnet, Q.-C. He, P. Li, A data-driven harmonic approach to
    constructing anisotropic damage models with a minimum number of internal
    variables, Journal of the Mechanics and Physics of Solids,
    162:104828, 2022.
    https://doi.org/10.1016/j.jmps.2022.104828
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# --- import loadmat for MATLAB .mat files ---
try:
    from scipy.io import loadmat
except Exception as e:
    raise SystemExit(
        "SciPy is not available in this interpreter.\n"
        "Install it with: python -m pip install scipy\n"
        f"Details: {e}"
    )

print(">>> Running:", os.path.abspath(__file__))


# ============================================================
# Utility functions (menu, database loading)
# ============================================================

def show_menu():
    print("=" * 30)
    print("1: RVE A: RVE with vertical weak layer, loading along x")
    print("2: RVE B: RVE with circular hole, loading along x")
    print("3: RVE C: RVE with circular inclusion, loading along x")
    print("4: RVE D: concrete microstructure (image), loading along x")
    print("5: RVE E: RVE with 45° weak layer, loading along x")
    print("=" * 30)


def ask_type():
    """Ask the user to select a database index (1-5)."""
    while True:
        try:
            choice = int(input("Choose database (1-5): ").strip())
            if choice in {1, 2, 3, 4, 5}:
                return choice
        except ValueError:
            pass
        print("Invalid input. Please choose an integer between 1 and 5.")


def _find_variable(data_dict, possible_names):
    """Find a matching key (case insensitive) in a MATLAB dictionary."""
    lower_map = {k.lower(): k for k in data_dict.keys()}
    for name in possible_names:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def load_database(choice: int):
    """Load the selected .mat database and return CC_0 and CC_save tensors."""
    mapping = {
        1: "RVEA_2D.mat",
        2: "RVEB_2D.mat",
        3: "RVEC_2D.mat",
        4: "RVED_2D.mat",
        5: "RVEE_2D.mat",
    }

    mat_file = mapping.get(choice)
    if mat_file is None:
        raise ValueError("Invalid choice.")

    if not os.path.exists(mat_file):
        raise FileNotFoundError(f"File {mat_file} not found.")

    data = loadmat(mat_file, squeeze_me=True, struct_as_record=False)

    key = _find_variable(data, ["CC_save2", "cc_save2", "CC_SAVE2"])
    if key is None:
        raise KeyError("Variable 'CC_save2' not found. "
                       f"Available variables: {list(data.keys())}")

    CC_save2 = np.asarray(data[key])
    if CC_save2.ndim == 2:
        CC_save2 = CC_save2[:, :, np.newaxis]

    if CC_save2.shape[0] != 3 or CC_save2.shape[1] != 3:
        raise ValueError("CC_save2 must be of shape (3,3,N).")

    CC_0 = CC_save2[:, :, 0]
    return {"name": mat_file, "CC_0": CC_0, "CC_save": CC_save2}


# ============================================================
# Computational functions
# ============================================================

def compute_eta_kappa(C3x3: np.ndarray, Npoints: int = 50):
    """Compute the η and κ functions from a given 3x3 elastic tensor."""
    if C3x3.shape != (3, 3):
        raise ValueError("C3x3 must be a 3x3 matrix.")

    C1111 = C3x3[0, 0]
    C2222 = C3x3[1, 1]
    C1212 = C3x3[2, 2]
    C1122 = C3x3[0, 1]
    C1112 = C3x3[0, 2]
    C2212 = C3x3[1, 2]

    thet = np.linspace(-np.pi, np.pi, Npoints)
    NT = len(thet) - 1

    eta = np.zeros((NT, 1))
    kappa = np.zeros((NT, 1))
    MatTeta = np.zeros((NT, 3))

    for i in range(NT):
        TT = thet[i]
        n1, n2 = np.cos(TT), np.sin(TT)

        eta[i, 0] = (C1111 * n1**4
                     + C2222 * n2**4
                     + 2.0 * (C1122 + 2.0 * C1212) * n1**2 * n2**2
                     + 4.0 * C1112 * n1**3 * n2
                     + 4.0 * C2212 * n1 * n2**3)

        kappa[i, 0] = (C1111 * n1**2
                       + C2222 * n2**2
                       + C1122 * n1**2
                       + C1122 * n2**2
                       + 2.0 * C1112 * n1 * n2
                       + 2.0 * C2212 * n1 * n2)

    return thet, eta, kappa, MatTeta


def precalcul_ope_2d(thet, eta0, kappa0):
    
    thet   = np.asarray(thet).ravel()
    eta0   = np.asarray(eta0).ravel()
    kappa0 = np.asarray(kappa0).ravel()

    Dthet = 2.0 * np.pi / (len(thet) - 1)

    Somme_vecud     = np.zeros((6,1))
    Somme_vecvd     = np.zeros((6,1))
    Somme_vecV11d   = np.zeros((6,1))
    Somme_vecV12d   = np.zeros((6,1))
    Somme_vecZ1111d = np.zeros((6,1))
    Somme_vecZ1112d = np.zeros((6,1))

    for i, TT in enumerate(thet[:-1]):
        n1 = np.cos(TT)
        n2 = np.sin(TT)

        F11 = n1**2 - 0.5
        F22 = n2**2 - 0.5
        F12 = n1 * n2

        F1111 = n1**4 - n1**2 + 1.0/8.0
        F1122 = n1**2 * n2**2 - 1.0/8.0
        F1112 = n1**3 * n2 - 0.5 * n1 * n2

        F2222 = n2**4 - n2**2 + 1.0/8.0
        F2212 = n1 * n2**3 - 0.5 * (n1 * n2)

        F1212 = n1**2 * n2**2 - 1.0/8.0

        VecD = np.array([
            [1.0],
            [0.0],
            [F11 - F22],
            [2.0 * F12],
            [F1111 - 2.0*F1122 + F2222 - 4.0*F1212],
            [4.0*F1112 - 4.0*F2212]
        ])

        VecH = np.array([
            [0.0],
            [1.0],
            [F11 - F22],
            [2.0 * F12],
            [0.0],
            [0.0]
        ])

        Somme_vecud     += eta0[i]   * VecD * Dthet
        Somme_vecvd     += kappa0[i] * VecH * Dthet
        Somme_vecV11d   += F11   * eta0[i] * VecD * Dthet
        Somme_vecV12d   += F12   * eta0[i] * VecD * Dthet
        Somme_vecZ1111d += F1111 * eta0[i] * VecD * Dthet
        Somme_vecZ1112d += F1112 * eta0[i] * VecD * Dthet

    vecud     = (1.0/(2.0*np.pi)) * Somme_vecud
    vecvd     = (1.0/(2.0*np.pi)) * Somme_vecvd
    vecV11d   = (2.0/np.pi)       * Somme_vecV11d
    vecV12d   = (2.0/np.pi)       * Somme_vecV12d
    vecZ1111d = (8.0/np.pi)       * Somme_vecZ1111d
    vecZ1112d = (8.0/np.pi)       * Somme_vecZ1112d

    C1111tt = vecud + vecV11d + vecZ1111d
    C1122tt = vecvd - vecud - vecZ1111d
    C1112tt = 0.5*vecV12d + vecZ1112d
    C2222tt = vecud - vecV11d + vecZ1111d
    C2212tt = 0.5*vecV12d - vecZ1112d
    C1212tt = vecud - 0.5*vecvd - vecZ1111d

   
    C_tilde = np.hstack([C1111tt, C1122tt, C1112tt,
                         C2222tt, C2212tt, C1212tt]).T

    return C_tilde


def calcul_var_direct_2d(C_tilde: np.ndarray, CC_0: np.ndarray, CC_save: np.ndarray) -> np.ndarray:

    def mat3_to_vec6(M: np.ndarray) -> np.ndarray:
        return np.array([[M[0,0]], [M[0,1]], [M[0,2]],
                         [M[1,1]], [M[1,2]], [M[2,2]]], dtype=float)

    C_0_vec = mat3_to_vec6(CC_0)
    K = C_tilde.T @ C_tilde

    N = CC_save.shape[2]
    SAVE_Alphas = np.zeros((N, 6))

    for i in range(N):
        Cx = CC_save[:, :, i]
        C_x_vec = mat3_to_vec6(Cx)
        G = C_tilde.T @ (C_0_vec - C_x_vec)
        alpha_i = np.linalg.solve(K, G)
        SAVE_Alphas[i, :] = alpha_i.ravel()

    return SAVE_Alphas


def plot_alphas(SAVE_Alphas: np.ndarray):
    """Plot the evolution of internal variables alpha_i."""
    plt.figure(figsize=(8, 5))
    plt.plot(SAVE_Alphas, linewidth=2)
    plt.xlabel("n", fontsize=16)
    plt.ylabel(r"$\alpha_i$", fontsize=16)
    plt.title("Internal variables $\\alpha_i$")
    plt.grid(True)
    plt.box(True)
    labels = [rf"$\alpha_{i+1}$" for i in range(6)]
    plt.legend(labels, loc="best")
    plt.show()


def compare_tensor_reconstruction(CC_0, CC_save, SAVE_Alphas, C_tilde, step=5):
    """Compare original and reconstructed tensor components."""
    def mat3_to_vec6(M: np.ndarray) -> np.ndarray:
        return np.array([M[0,0], M[0,1], M[0,2],
                         M[1,1], M[1,2], M[2,2]], dtype=float)

    N = CC_save.shape[2]
    C_vec0 = mat3_to_vec6(CC_0)
    C_vec = np.zeros((6, N))
    CC_save_REC = np.zeros((6, N))

    for i in range(N):
        CC = CC_save[:, :, i]
        alpha = SAVE_Alphas[i, :]
        C_vec[:, i] = mat3_to_vec6(CC)
        CC_save_REC[:, i] = C_vec0 - (C_tilde @ alpha.reshape(-1,1)).ravel()

    plt.figure(figsize=(8,5))
    h = plt.plot(C_vec.T, linewidth=2)
    colors = [line.get_color() for line in h]

    idx = np.arange(0, N, step)
    for j in range(6):
        plt.plot(idx, CC_save_REC[j, idx], 'o',
                 linestyle="none", linewidth=2,
                 markeredgecolor=colors[j])

    plt.xlabel("n", fontsize=16)
    plt.ylabel(r"$C_{ijkl}$", fontsize=16)
    plt.title("Elastic tensor components")
    plt.grid(True)
    plt.box(True)

    dummy_line, = plt.plot([], [], '-', color='k', linewidth=2)
    dummy_marker, = plt.plot([], [], 'o', markeredgecolor='k',
                             linestyle='none')
    plt.legend([dummy_line, dummy_marker],
               ["Data", "Reconstructed"], loc="best")

    plt.show()

    return C_vec, CC_save_REC


def plot_damage_surface(alpha: np.ndarray):
    """Plot the damage orientation function d(theta) in polar coordinates."""
    alpha = np.asarray(alpha, dtype=float).reshape(6,)
    theta = np.linspace(0, 2*np.pi, 360, endpoint=False)
    d_theta = np.empty_like(theta)

    for i, th in enumerate(theta):
        n1, n2 = np.cos(th), np.sin(th)

        F11 = n1**2 - 0.5
        F22 = n2**2 - 0.5
        F12 = n1 * n2

        F1111 = n1**4 - n1**2 + 1/8
        F1122 = n1**2 * n2**2 - 1/8
        F2222 = n2**4 - n2**2 + 1/8
        F1212 = n1**2 * n2**2 - 1/8

        F1112 = n1**3 * n2 - 0.5 * n1 * n2
        F2212 = n1 * n2**3 - 0.5 * n1 * n2

        Vd = np.array([
            1.0,
            0.0,
            F11 - F22,
            2.0 * F12,
            F1111 - 2*F1122 + F2222 - 4*F1212,
            4.0 * (F1112 - F2212),
        ], dtype=float)

        d_theta[i] = float(Vd @ alpha)

    ax = plt.subplot(111, projection="polar")
    ax.plot(theta, d_theta, 'r-', linewidth=2)
    ax.set_title(r"Damage orientation function $d(\theta)$", va="bottom")
    ax.set_rlim(0, 1)
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.grid(True)
    plt.show()


def reduced_variables_pca(SAVE_Alphas: np.ndarray, error_tol: float = 0.01):
    """
    Reduction of internal variables using PCA (Principal Component Analysis).
    
    Parameters:
        SAVE_Alphas : array (N x D) of alpha_i
        error_tol   : tolerance on variance loss (e.g. 0.01 = 1%)
    
    Returns:
        X_reduced : projected data (N x M)
        W         : principal component vectors (D x M)
    """
    X = np.asarray(SAVE_Alphas, dtype=float)
    N, D = X.shape

    # Covariance (not centered, as in MATLAB)
    C = (X.T @ X) / N

    # SVD
    U, S, _ = np.linalg.svd(C)

    # Eigenvalues and explained variance
    eigenvalues = S
    explained = 100 * eigenvalues / np.sum(eigenvalues)
    cumulative = np.cumsum(explained)

    # Minimum number of components M
    P = 1 - error_tol
    M = np.argmax(cumulative >= P*100) + 1

    print(f">> PCA: {M} component(s) kept (tolerance {error_tol*100:.1f}%)")

    W = U[:, :M]       # (D x M)
    X_reduced = X @ W  # (N x M)

    plt.figure(figsize=(8,5))
    colors = plt.cm.get_cmap("tab10", M)

    for i in range(M):
        plt.plot(np.arange(1, N+1), -X_reduced[:, i],
                 linewidth=2, color=colors(i),
                 label=fr"$\beta_{i+1}$")

    plt.xlabel("n", fontsize=16)
    plt.ylabel(r"$\beta_i$", fontsize=16)
    plt.title("Reduced internal variable(s) $\\beta_i$ by PCA")
    plt.grid(True)
    plt.box(True)
    plt.legend(loc="best")
    plt.show()

    return X_reduced, W


# ============================================================
# Main program
# ============================================================

def main():
    show_menu()
    choice = ask_type()

    base = load_database(choice)
    print(">> Database loaded:", base["name"])
    print("CC_0 shape:", base["CC_0"].shape)
    print("CC_save shape:", base["CC_save"].shape)

    # --- Compute eta, kappa and C_tilde ---
    C3 = base["CC_0"][:3, :3]
    thet, eta0, kappa0, _ = compute_eta_kappa(C3, Npoints=50)
    C_tilde = precalcul_ope_2d(thet, eta0, kappa0)

    # --- Compute alphas ---
    SAVE_Alphas = calcul_var_direct_2d(C_tilde, base["CC_0"], base["CC_save"])

    # --- Plot alpha evolution ---
    plot_alphas(SAVE_Alphas)

    # --- Compare tensor components (data vs reconstruction) ---
    compare_tensor_reconstruction(base["CC_0"], base["CC_save"], SAVE_Alphas, C_tilde)

    # --- Plot damage surface using last alpha ---
    alpha_last = SAVE_Alphas[-1, :]
    plot_damage_surface(alpha_last)

    # --- PCA reduction of internal variables ---
    reduced_variables_pca(SAVE_Alphas, error_tol=0.01)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
DDHADtens 3D v1 

J. Yvonnet July 1st, 2025
Gustave Eiffel University
Input : a collection of 6x6 elastic tensors
Output: associated collection of damage variables (harmonic analysis)
21x1 vector
plots of damage variable evolution and damage surfaces (2D)

The database is associated with the numerical examples in 
[2] J. Yvonnet, Qi-Chang He, Microstructure-based machine learning of damage
 models including anisotropy, irreversibility and evolution, 
 Journal of the Mechanics and Physics of Solids, 106160, 2025.

"""

import os
import numpy as np
import matplotlib.pyplot as plt

# SciPy for .mat IO
try:
    from scipy.io import loadmat
except Exception as e:
    raise SystemExit(
        "SciPy is required. Install: python -m pip install scipy\n"
        f"Details: {e}"
    )

# ---------------------------
# UI helpers
# ---------------------------
def show_menu():
    print("="*30)
    print("1: RVE A: brittle matrix, hard inclusion, traction")
    print("2: RVE B: weak layer ⟂ x-axis, traction")
    print("3: RVE C: porous RVE, compression")
    print("4: RVE D: fiber along x in brittle matrix, traction")
    print("="*30)

def ask_type():
    while True:
        try:
            t = int(input("Choose database (1-4): ").strip())
            if t in {1,2,3,4}: return t
        except ValueError:
            pass
        print("Invalid. Enter 1..4.")

def _find_key(d, names):
    lmap = {k.lower(): k for k in d.keys()}
    for n in names:
        if n.lower() in lmap: return lmap[n.lower()]
    return None

def load_database(choice:int):
    mapping = {1:"RVEA_3D.mat", 2:"RVEB_3D.mat", 3:"RVEC_3D.mat", 4:"RVED_3D.mat"}
    fn = mapping[choice]
    if not os.path.exists(fn):
        raise FileNotFoundError(f"{fn} not found")
    m = loadmat(fn, squeeze_me=True, struct_as_record=False)
    k = _find_key(m, ["CC_save2","cc_save2","CC_SAVE2"])
    if k is None: raise KeyError(f"'CC_save2' not found. Keys: {list(m.keys())}")
    CC_save2 = np.asarray(m[k], dtype=float)
    if CC_save2.ndim == 2: CC_save2 = CC_save2[:,:,None]
    if CC_save2.shape[:2] != (6,6): raise ValueError("CC_save2 must be (6,6,N)")
    CC_0 = CC_save2[:,:,0]
    return {"name": fn, "CC_0": CC_0, "CC_save": CC_save2}

# ---------------------------
# Voigt mapping (6x6 -> 21)
# Order exactly as in the MATLAB code
# ---------------------------
def mat6_to_vec21(C):
    return np.array([
        C[0,0], C[0,1], C[0,2], C[0,3], C[0,4], C[0,5],
        C[1,1], C[1,2], C[1,3], C[1,4], C[1,5],
        C[2,2], C[2,3], C[2,4], C[2,5],
        C[3,3], C[3,4], C[3,5],
        C[4,4], C[4,5],
        C[5,5]
    ], dtype=float).reshape(21,1)

# ---------------------------
# Step 2: eta^0(θ,φ), kappa^0(θ,φ) on the sphere
# ---------------------------
def compute_eta_kappa_3d(C6x6, NbPts=50):
    # Extract needed components (Voigt 6x6, same naming as MATLAB)
    C1111=C6x6[0,0]; C1122=C6x6[0,1]; C1133=C6x6[0,2]; C1112=C6x6[0,3]; C1113=C6x6[0,4]; C1123=C6x6[0,5]
    C2222=C6x6[1,1]; C2233=C6x6[1,2]; C2212=C6x6[1,3]; C2213=C6x6[1,4]; C2223=C6x6[1,5]
    C3333=C6x6[2,2]; C3312=C6x6[2,3]; C3313=C6x6[2,4]; C3323=C6x6[2,5]
    C1212=C6x6[3,3]; C1213=C6x6[3,4]; C1223=C6x6[3,5]
    C1313=C6x6[4,4]; C1323=C6x6[4,5]
    C2323=C6x6[5,5]

    dth = np.pi/(NbPts)
    thet = np.linspace(0, np.pi, NbPts+1)  # include endpoints, we will use up to -1
    phi  = np.linspace(0, 2*np.pi, NbPts+1)

    NT = len(thet)-1
    MP = len(phi)-1

    eta   = np.zeros((NT,MP))
    kappa = np.zeros((NT,MP))

    for i in range(NT):
        TT = thet[i]
        sT, cT = np.sin(TT), np.cos(TT)
        for j in range(MP):
            PP = phi[j]
            cP, sP = np.cos(PP), np.sin(PP)
            n1, n2, n3 = sT*cP, sT*sP, cT

            # eta (quartic) — matches the MATLAB expression literally
            eta[i,j] = (C1111*n1**4 + C2222*n2**4 + C3333*n3**4
                        + 2*C1122*n1**2*n2**2 + 2*C1133*n1**2*n3**2 + 2*C2233*n2**2*n3**2
                        + 4*C1212*n1**2*n2**2 + 4*C1313*n1**2*n3**2 + 4*C2323*n2**2*n3**2
                        + 4*C1112*n1**3*n2 + 4*C1113*n1**3*n3 + 4*C2212*n2**3*n1 + 4*C3313*n1*n3**3
                        + 4*C1123*n1**2*n2*n3 + 4*C2213*n2**2*n1*n3 + 4*C3312*n1*n2*n3**2
                        + 4*C2223*n2**3*n3 + 4*C3323*n2*n3**3
                        + 8*C1213*n1**2*n2*n3 + 8*C1223*n1*n2**2*n3 + 8*C1323*n1*n2*n3**2)

            # kappa (quadratic) — careful: exact MATLAB terms (no (n1*n2) mixup)
            kappa[i,j] = (C1111*n1**2 + C2222*n2**2 + C3333*n3**2
                          + C1122*n1**2 + C1122*n2**2
                          + C1133*n3**2 + C1133*n1**2
                          + C2233*n3**2 + C2233*n2**2
                          + 2*C1112*n1*n2 + 2*C1113*n1*n3 + 2*C1123*n2*n3
                          + 2*C2212*n1*n2 + 2*C2213*n1*n3 + 2*C2223*n2*n3
                          + 2*C3312*n1*n2 + 2*C3313*n1*n3 + 2*C3323*n2*n3)

    return thet, phi, eta, kappa

# ---------------------------
# Build C_tilde (21x21)
# ---------------------------
def precalcul_ope_3d(thet, phi, eta0, kappa0):
    NT, MP = eta0.shape
    Dthet = np.pi/(len(thet)-1)
    Dphi  = 2*np.pi/(len(phi)-1)

    # sums (21x1)
    def v(): return np.zeros((21,1))
    Som_a   = v(); Som_b   = v()
    Som_U11 = v(); Som_U12 = v(); Som_U13 = v(); Som_U22 = v(); Som_U23 = v()
    Som_V11 = v(); Som_V12 = v(); Som_V13 = v(); Som_V22 = v(); Som_V23 = v()
    Som_Z1111=v(); Som_Z1122=v(); Som_Z2222=v(); Som_Z1112=v(); Som_Z1113=v()
    Som_Z1123=v(); Som_Z2212=v(); Som_Z2213=v(); Som_Z2223=v()

    for i in range(NT):
        TT = thet[i]
        sT, cT = np.sin(TT), np.cos(TT)
        wt = Dthet*Dphi*sT
        for j in range(MP):
            PP = phi[j]
            cP, sP = np.cos(PP), np.sin(PP)
            n1, n2, n3 = sT*cP, sT*sP, cT

            F11 = n1**2 - 1/3; F22 = n2**2 - 1/3; F33 = n3**2 - 1/3
            F12 = n1*n2; F13 = n1*n3; F23 = n2*n3

            F1111 = n1**4 - (6/7)*n1**2 + 3/35
            F1122 = n1**2*n2**2 - (1/7)*(n1**2+n2**2) + 1/35
            F1133 = n1**2*n3**2 - (1/7)*(n1**2+n3**2) + 1/35
            F1112 = n1**3*n2 - (3/7)*n1*n2
            F1113 = n1**3*n3 - (3/7)*n1*n3
            F1123 = n1**2*n2*n3 - (1/7)*n2*n3

            F2222 = n2**4 - (6/7)*n2**2 + 3/35
            F2233 = n2**2*n3**2 - (1/7)*(n2**2+n3**2) + 1/35
            F2212 = n1*n2**3 - (3/7)*n1*n2
            F2213 = n1*n2**2*n3 - (1/7)*n1*n3
            F2223 = n2**3*n3 - (3/7)*n2*n3

            F3333 = n3**4 - (6/7)*n3**2 + 3/35
            F3312 = n1*n2*n3**2 - (1/7)*n1*n2
            F3313 = n1*n3**3 - (3/7)*n1*n3
            F3323 = n2*n3**3 - (3/7)*n2*n3

            F1212 = n1**2*n2**2 - (1/7)*(n1**2+n2**2) + 1/35
            F1213 = n1**2*n2*n3 - (1/7)*n2*n3
            F1223 = n1*n2**2*n3 - (1/7)*n1*n3
            F1313 = n1**2*n3**2 - (1/7)*(n1**2+n3**2) + 1/35
            F1323 = n1*n2*n3**2 - (1/7)*n1*n2
            F2323 = n2**2*n3**2 - (1/7)*(n2**2+n3**2) + 1/35

            VecD = np.array([
                [1.0],
                [0.0],
                [F11-F33],
                [2*F12],
                [2*F13],
                [F22-F33],
                [2*F23],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [F1111+F3333-2*F1133-4*F1313],
                [2*F3333+2*F1122-2*F1133-2*F2233+4*F1212-4*F1313-4*F2323],
                [F2222+F3333-2*F2233-4*F2323],
                [4*F1112-4*F3312-8*F1223-8*F1323],
                [4*F1113-4*F3313],
                [4*F1123-4*F3323+8*F1213],
                [4*F2212-4*F3312-8*F1223-8*F1323],
                [-4*F3313+4*F2213],
                [4*F2223-4*F3323]
            ])

            VecH = np.array([
                [0.0],
                [1.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [F11-F33],
                [2*F12],
                [2*F13],
                [F22-F33],
                [2*F23],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0],
                [0.0]
            ])

            e = eta0[i,j]; k = kappa0[i,j]
            w = wt

            Som_a   += e*VecD*w
            Som_b   += k*VecH*w

            Som_U11 += e*VecD*F11*w
            Som_U12 += e*VecD*F12*w
            Som_U13 += e*VecD*F13*w
            Som_U22 += e*VecD*F22*w
            Som_U23 += e*VecD*F23*w

            Som_V11 += k*VecH*F11*w
            Som_V12 += k*VecH*F12*w
            Som_V13 += k*VecH*F13*w
            Som_V22 += k*VecH*F22*w
            Som_V23 += k*VecH*F23*w

            Som_Z1111 += e*VecD*F1111*w
            Som_Z1122 += e*VecD*F1122*w
            Som_Z2222 += e*VecD*F2222*w
            Som_Z1112 += e*VecD*F1112*w
            Som_Z1113 += e*VecD*F1113*w
            Som_Z1123 += e*VecD*F1123*w
            Som_Z2212 += e*VecD*F2212*w
            Som_Z2213 += e*VecD*F2213*w
            Som_Z2223 += e*VecD*F2223*w

    # Normalizations (exact MATLAB factors)
    Som_a   *= (1/(4*np.pi)); Som_b   *= (1/(4*np.pi))
    factor15 = (15/(8*np.pi))
    Som_U11 *= factor15; Som_U12 *= factor15; Som_U13 *= factor15; Som_U22 *= factor15; Som_U23 *= factor15
    Som_V11 *= factor15; Som_V12 *= factor15; Som_V13 *= factor15; Som_V22 *= factor15; Som_V23 *= factor15
    factor315 = (315/(32*np.pi))
    Som_Z1111*= factor315; Som_Z1122*= factor315; Som_Z2222*= factor315
    Som_Z1112*= factor315; Som_Z1113*= factor315; Som_Z1123*= factor315
    Som_Z2212*= factor315; Som_Z2213*= factor315; Som_Z2223*= factor315

    alpha = 0.5*(Som_b - Som_a)
    beta  = 0.25*(3*Som_a - Som_b)

    AA11 = Som_V11 - Som_U11
    AA22 = Som_V22 - Som_U22
    AA12 = Som_V12 - Som_U12
    AA13 = Som_V13 - Som_U13
    AA23 = Som_V23 - Som_U23

    BB11 = 0.25*(3*Som_U11 - 2*Som_V11)
    BB22 = 0.25*(3*Som_U22 - 2*Som_V22)
    BB12 = 0.25*(3*Som_U12 - 2*Som_V12)
    BB13 = 0.25*(3*Som_U13 - 2*Som_V13)
    BB23 = 0.25*(3*Som_U23 - 2*Som_V23)

    C1111tt = alpha + 2*beta + 2*AA11 + 4*BB11 + Som_Z1111
    C1122tt = alpha + (AA11+AA22) + Som_Z1122
    C1133tt = alpha - AA22 - Som_Z1111 - Som_Z1122
    C1112tt = AA12 + 2*BB12 + Som_Z1112
    C1113tt = AA13 + 2*BB13 + Som_Z1113
    C1123tt = AA23 + Som_Z1123

    C2222tt = alpha + 2*beta + 2*AA22 + 4*BB22 + Som_Z2222
    C2233tt = alpha - AA11 - Som_Z1122 - Som_Z2222
    C2212tt = AA12 + 2*BB12 + Som_Z2212
    C2213tt = AA13 + Som_Z2213
    C2223tt = AA23 + 2*BB23 + Som_Z2223

    C3333tt = alpha + 2*beta - 2*(AA11+AA22) - 4*(BB11+BB22) + Som_Z1111 + Som_Z2222 + 2*Som_Z1122
    C3312tt = AA12 - Som_Z1112 - Som_Z2212
    C3313tt = AA13 + 2*BB13 - Som_Z1113 - Som_Z2213
    C3323tt = AA23 + 2*BB23 - Som_Z1123 - Som_Z2223

    C1212tt = beta + (BB11+BB22) + Som_Z1122
    C1213tt = BB23 + Som_Z1123
    C1223tt = BB13 + Som_Z2213

    C1313tt = beta - BB22 - Som_Z1111 - Som_Z1122
    C1323tt = BB12 - Som_Z1112 - Som_Z2212

    C2323tt = beta - BB11 - Som_Z1122 - Som_Z2222

    # Stack then transpose (MATLAB: [ ... ]')
    C_tilde = np.hstack([
        C1111tt, C1122tt, C1133tt, C1112tt, C1113tt, C1123tt,
        C2222tt, C2233tt, C2212tt, C2213tt, C2223tt,
        C3333tt, C3312tt, C3313tt, C3323tt,
        C1212tt, C1213tt, C1223tt, C1313tt, C1323tt, C2323tt
    ]).T  # (21x21)

    return C_tilde

# ---------------------------
# Solve for alphas over CC_save
# ---------------------------
def calcul_var_direct_3d(C_tilde, CC_0, CC_save):
    C_0_vec = mat6_to_vec21(CC_0)        # (21,1)
    K = C_tilde.T @ C_tilde              # (21x21)

    N = CC_save.shape[2]
    SAVE_Alphas = np.zeros((N,21))

    for i in range(N):
        Cx = CC_save[:,:,i]
        C_x_vec = mat6_to_vec21(Cx)                 # (21,1)
        G = C_tilde.T @ (C_0_vec - C_x_vec)         # (21,1)
        alpha_i = np.linalg.solve(K, G)             # (21,1)
        SAVE_Alphas[i,:] = alpha_i.ravel()
    return SAVE_Alphas

# ---------------------------
# Plots
# ---------------------------
def plot_alphas(SAVE_Alphas):
    plt.figure(figsize=(9,5))
    plt.plot(SAVE_Alphas, linewidth=1.8)
    plt.xlabel("n", fontsize=14); plt.ylabel(r"$\alpha_i$", fontsize=14)
    plt.title("Internal variables $\\alpha_i$")
    plt.grid(True); plt.box(True)
    plt.show()

def compare_reconstruction(CC_0, CC_save, SAVE_Alphas, C_tilde, step=5):
    N = CC_save.shape[2]
    C_vec0 = mat6_to_vec21(CC_0).ravel()
    C_vec = np.zeros((21,N))
    C_rec = np.zeros((21,N))
    for i in range(N):
        C = CC_save[:,:,i]
        a = SAVE_Alphas[i,:]
        C_vec[:,i] = mat6_to_vec21(C).ravel()
        C_rec[:,i] = C_vec0 - (C_tilde @ a.reshape(-1,1)).ravel()
    plt.figure(figsize=(9,5))
    h = plt.plot(C_vec.T, linewidth=1.8)
    cols = [ln.get_color() for ln in h]
    idx = np.arange(0, N, step)
    for j in range(21):
        plt.plot(idx, C_rec[j,idx], 'o', linestyle="none", linewidth=2,
                 markeredgecolor=cols[j])
    plt.xlabel("n", fontsize=14); plt.ylabel(r"$C_{ijkl}$", fontsize=14)
    plt.title("Elastic tensor components")
    plt.grid(True); plt.box(True)
    dl, = plt.plot([],[],'-',color='k',linewidth=2)
    dm, = plt.plot([],[],'o',markeredgecolor='k',linestyle='none')
    plt.legend([dl,dm],["Data","Reconstructed"],loc="best")
    plt.show()

def plot_damage_surface_3d(alpha, n_theta=200):
    # Grid on sphere
    theta, phi = np.meshgrid(np.linspace(0,np.pi,n_theta),
                             np.linspace(0,2*np.pi,n_theta),
                             indexing='ij')
    n1 = np.sin(theta)*np.cos(phi)
    n2 = np.sin(theta)*np.sin(phi)
    n3 = np.cos(theta)

    alpha = np.asarray(alpha, float).reshape(21,)
    d_vals = np.zeros_like(n1)

    # Build F2 and F4 per point and compute Vd⋅alpha
    # (faithful to the MATLAB doc construction)
    I = np.eye(3)
    def delta(i,j): return 1.0 if i==j else 0.0

    for i in range(n1.shape[0]):
        for j in range(n1.shape[1]):
            n = np.array([n1[i,j], n2[i,j], n3[i,j]])
            N = np.outer(n,n)
            F2 = N - (1/3)*I

            # F4 (3x3x3x3)
            F4 = np.zeros((3,3,3,3))
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        for d in range(3):
                            term1 = N[a,b]*N[c,d]
                            term2 = (1/7)*( delta(a,b)*N[c,d] + N[a,b]*delta(c,d)
                                           +delta(a,d)*N[b,c] + N[a,d]*delta(b,c)
                                           +delta(a,c)*N[b,d] + N[a,c]*delta(b,d) )
                            term3 = (1/35)*( delta(a,b)*delta(c,d)
                                            +delta(a,d)*delta(b,c)
                                            +delta(a,c)*delta(b,d) )
                            F4[a,b,c,d] = term1 - term2 + term3

            getF = lambda a,b,c,d: F4[a,b,c,d]

            Vd = np.zeros(21)
            # indices in 0-based (1->0, 2->1, 3->2)
            Vd[0] = 1.0
            # second-order (Voigt-symmetric)
            Vd[2] = F2[0,0]-F2[2,2]
            Vd[3] = 2*F2[0,1]
            Vd[4] = 2*F2[0,2]
            Vd[5] = F2[1,1]-F2[2,2]
            Vd[6] = 2*F2[1,2]
            # the block Vh isn't used in plotting (only Vd•alpha)

            # fourth-order combinations (exactly as in MATLAB)
            Vd[12] = getF(0,0,0,0) + getF(2,2,2,2) - 2*getF(0,0,2,2) - 4*getF(0,2,0,2)
            Vd[13] = 2*getF(2,2,2,2) + 2*getF(0,0,1,1) - 2*getF(0,0,2,2) - 2*getF(1,1,2,2) \
                     + 4*getF(0,1,0,1) - 4*getF(0,2,0,2) - 4*getF(1,2,1,2)
            Vd[14] = getF(1,1,1,1) + getF(2,2,2,2) - 2*getF(1,1,2,2) - 4*getF(1,2,1,2)
            Vd[15] = 4*( getF(0,0,0,1) - getF(2,2,0,1) - 2*getF(0,1,1,2) - 2*getF(0,2,1,2) )
            Vd[16] = 4*( getF(0,0,0,2) - getF(2,2,0,2) )
            Vd[17] = 4*( getF(0,0,1,2) - getF(2,2,1,2) + 2*getF(0,1,0,2) )
            Vd[18] = 4*( getF(1,1,0,1) - getF(2,2,0,1) - 2*getF(0,1,1,2) - 2*getF(0,2,1,2) )
            Vd[19] = 4*( getF(1,1,0,2) - getF(2,2,0,2) )
            Vd[20] = 4*( getF(1,1,1,2) - getF(2,2,1,2) )

            d_vals[i,j] = float(Vd @ alpha)

    R = np.abs(d_vals)
    X = R*n1; Y = R*n2; Z = R*n3

    plt.figure(figsize=(7,6))
    surf = plt.axes(projection='3d').plot_surface(X, Y, Z, facecolors=None,
                                                  rstride=1, cstride=1,
                                                  linewidth=0, antialiased=True,
                                                  shade=True)
    ax = plt.gca()
    ax.set_box_aspect([1,1,1])
    plt.title(r"Damage orientation function $d(\theta,\phi)$")
    ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
    plt.show()

# ---------------------------
# PCA reduction (N x 21 -> N x M)
# ---------------------------
def reduced_variables_pca(SAVE_Alphas, error_tol=0.01):
    """
    PCA reduction of internal variables (mirror of the MATLAB code).

    NOTE for users:
    ---------------
    The number of reduced variables M obtained by PCA may differ
    slightly between this Python version and the original MATLAB code.
    This may be a consequence of:

    1. Different numerical solvers (LAPACK in NumPy vs. MATLAB SVD).
    2. Very small eigenvalues close to machine precision (~1e-14),
       which can change the cumulative explained variance around
       the 99% threshold.
    3. The tolerance criterion (Error = 0.01) being strict:
       tiny numerical differences may lead to keeping one more or one
       fewer principal component.

    => If you need strict consistency with MATLAB, you can either:
       - round the cumulative explained variance (e.g. np.round(..., 12)),
       - or manually set M to the same value observed in MATLAB.
    """
    X = np.asarray(SAVE_Alphas, float)
    N, D = X.shape
    C = (X.T @ X) / N   # covariance

    # Use SVD as in MATLAB
    U, S, _ = np.linalg.svd(C)  # S = singular values
    eigenvalues = S             # for symmetric positive C, these = eigenvalues

    explained = eigenvalues / np.sum(eigenvalues)   # fractions (0-1)
    cumulative = np.cumsum(explained)

    P = 1 - error_tol   # e.g. 0.99
    M = np.argmax(cumulative >= P) + 1

    print(f">> PCA: keeping {M} component(s) (tolerance {error_tol*100:.1f}%)")

    # Projection on first M components
    W = U[:, :M]
    Xr = X @ W

    # Plot reduced variables
    plt.figure(figsize=(9,5))
    cmap = plt.cm.get_cmap("tab10", M)
    for i in range(M):
        plt.plot(np.arange(1, N+1), -Xr[:, i], linewidth=2,
                 color=cmap(i), label=fr"$\beta_{i+1}$")
    plt.xlabel("n", fontsize=14)
    plt.ylabel(r"$\beta_i$", fontsize=14)
    plt.title("Reduced internal variable(s) $\\beta_i$ by PCA")
    plt.grid(True); plt.box(True); plt.legend(loc="best")
    plt.show()

    return Xr, W





# ---------------------------
# Main
# ---------------------------
def main():
    show_menu()
    choice = ask_type()
    base = load_database(choice)
    print(">> Database loaded:", base["name"])
    print("CC_0 shape:", base["CC_0"].shape, "| CC_save shape:", base["CC_save"].shape)

    # (1) angular discretization & (2) eta/kappa
    NbPts = 50
    thet, phi, eta0, kappa0 = compute_eta_kappa_3d(base["CC_0"], NbPts=NbPts)

    # (2) -> C_tilde (21x21)
    C_tilde = precalcul_ope_3d(thet, phi, eta0, kappa0)

    # (2b) alphas for each saved C
    SAVE_Alphas = calcul_var_direct_3d(C_tilde, base["CC_0"], base["CC_save"])

    # plots
    plot_alphas(SAVE_Alphas)
    compare_reconstruction(base["CC_0"], base["CC_save"], SAVE_Alphas, C_tilde)

    # (4) 3D damage surface from last alpha
    alpha_last = SAVE_Alphas[-1,:]
    plot_damage_surface_3d(alpha_last, n_theta=120)

    # (5) PCA reduction
    reduced_variables_pca(SAVE_Alphas, error_tol=0.01)

if __name__ == "__main__":
    main()


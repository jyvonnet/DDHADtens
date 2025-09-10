# DDHADtens
The present DDHADTENS program extracts variables characterizing anisotropic damage of elastic tensors. The code is based on the DDHAD method [1-3], which uses harmonic analysis of elastic tensor to define anisotropic damage  [4-5]. The program takes as input: (a) an undamaged elastic tensor and (b) a series of damaged elastic tensors. The output is a series of damage variables, defined as the coefficients of the spherical harmonics describing scalar damage orientation functions, which can be interpreted as the damage of elastic coefficients varying with the direction. These tensors can be obtained externally by numerical fracture simulations in a Representative Volume Element. A 2D version of the code is first provided. In this case, the elastic tensor is given by a 3×3 matrix and 6 damage variables are obtained. Second, a 3D version of the code is proposed, where the elastic tensor is given by a 6 × 6 matrix and 21 damage variables are obtained. In both case, a set of reduced variables using PCA is also provided

References

[1] Yvonnet, J., He, Q. C., & Li, P. (2023). Reducing internal variables and improving efficiency in data-driven modelling of anisotropic damage from RVE simulations. Computational Mechanics, 72(1), 37-55.

[2] Yvonnet, J., He, Q. C., & Li, P. (2022). A data-driven harmonic approach to constructing anisotropic damage models with a minimum number of internal variables. Journal of the Mechanics and Physics of Solids, 162, 104828.

[3] Yvonnet, J., & He, Q. C. (2025). Microstructure-based machine learning of damage models including anisotropy, irreversibility and evolution. Journal of the Mechanics and Physics of Solids, 106160.

[4] P. Ladevèze, Sur une théorie de l’endommagement anisotrope, Rapport interne No. 34, Laboratoire de M´ecanique et Technologie (1983).

[5] He, Q. C., & Curnier, A. (1995). A more fundamental approach to damaged elastic stress-strain relations. International Journal of Solids and Structures, 32(10), 1433-1457.

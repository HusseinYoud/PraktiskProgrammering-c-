This project solves the generalized eigenvalue problem

[
AV=BVE
]

for symmetric (A) and positive-definite (B).

Part A

The problem is transformed into an ordinary eigenvalue problem using an eigenvalue decomposition of (B). The implementation is verified by checking

[
V^TBV=I
]

and

[
AV=BVE.
]

Part B

The solver is used in a variational calculation of the hydrogen ground-state energy with Gaussian basis functions

[
\phi_i(r)=r e^{-\alpha_i r^2}.
]

The Gaussian exponents are optimized using a downhill simplex minimizer for (n=1,\ldots,5) basis functions.

The calculated energy converges towards the exact value

[
E_0=-0.5,
]

with

[
E_{n=5}\approx -0.49981.
]

The overlap matrix becomes increasingly ill-conditioned as more Gaussian functions are added.

Part C

A second generalized eigenvalue solver is implemented using the Cholesky decomposition

[
B=LL^T.
]

The transformed problem

[
L^{-1}AL^{-T}Y=YE
]

is solved, and the generalized eigenvectors are recovered using triangular substitution.

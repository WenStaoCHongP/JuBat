# Thermal verification: FEM vs polar FVM

This note summarizes the thermal verification model and compares the continuous PDE, FEM weak form, and polar finite-volume discretization used in the ring test.

## Model overview

We consider a hollow cylinder (annulus) with anisotropic conductivity in polar coordinates and uniform volumetric heat generation. The inner boundary is adiabatic and the outer boundary is convective.

## Comparison table

| Item | Continuous equation (polar) | FEM weak form (Q4) | Polar FVM (r-θ control volumes) |
|---|---|---|---|
| Governing PDE | $\rho c \partial_t T = \frac{1}{r}\partial_r(k_r r \partial_r T) + \frac{1}{r^2}\partial_\theta(k_\theta \partial_\theta T) + q$ | Same PDE enforced in weak form | Same PDE enforced by flux balance |
| Weak/Integral statement | N/A | $\int_\Omega \rho c N_i \partial_t T\, d\Omega + \int_\Omega (\nabla N_i)^T K \nabla T\, d\Omega = \int_\Omega N_i q\, d\Omega + \int_{\Gamma_h} N_i h (T_\infty - T)\, d\Gamma$ | Control volume balance: $\rho c V_{i,j} \dot T_{i,j} = Q_{i,j} + (\dot Q_{r,i-1/2}-\dot Q_{r,i+1/2}) + (\dot Q_{\theta,j-1/2}-\dot Q_{\theta,j+1/2})$ |
| Radial flux | $-k_r \partial_r T$ | Embedded in $K=B^T K B$ (via $\nabla N$) | $\dot Q_{r,i+1/2} = k_r (r_{i+1/2}\Delta\theta)\frac{T_{i+1,j}-T_{i,j}}{\Delta r_{i+1/2}}$ |
| Tangential flux | $-(k_\theta/r) \partial_\theta T$ | Embedded in $K=B^T K B$ (via $\nabla N$) | $\dot Q_{\theta,j+1/2} = k_\theta \frac{\Delta r_i}{r_i \Delta\theta}(T_{i,j+1}-T_{i,j})$ |
| Mass term | N/A | $M_{ij}=\int_\Omega \rho c N_i N_j \, d\Omega$ | $M_{i,j}=\rho c V_{i,j}$ |
| Inner BC | $\partial_r T=0$ at $r=R_{in}$ | Natural BC (no flux term) | Omit inner radial flux term |
| Outer BC | $-k_r \partial_r T = h(T-T_\infty)$ at $r=R_{out}$ | Adds $\int_{\Gamma_h} N_i h(T_\infty - T) d\Gamma$ | Add $-h A$ on diagonal and $h A T_\infty$ to RHS at outer ring |

## Notes

- The FEM and polar FVM solve the same physics but use different discrete operators, which can lead to different symmetry behavior on coarse meshes.
- The polar FVM respects rotational symmetry on the r-θ grid, providing a clean axisymmetric reference for verification.

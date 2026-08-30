# Findings & Decisions

## Requirements
- Explain how viscous regularization can address the CZM non-convergence around D_max > 0.89.
- If needed, implement a local code change in the CZM path and verify it on the coupled example.

## Research Findings
- The coupled example currently stalls in post-peak softening, not because of NaNs.
- Diagnostic output showed D_max rising until about 0.8949, after which arc-length or basic iteration can stall with a large residual.
- The current CZM law is a rate-independent bilinear law; the tangent becomes very soft in the softening branch and can make the global system nearly singular.
- The non-convergence is most likely a continuation/path-following failure amplified by strong electro-thermal coupling.
- A useful viscous regularization must act on the cohesive history variable / traction law itself; changing only the global step size or arc-length control is not enough to smooth the local softening tangent.
- The current CZM solve only commits damage states when the nonlinear solve converges, so any viscous history variable must be rollback-safe across failed steps.

## Technical Decisions
| Decision | Rationale |
|----------|-----------|
| Prefer viscous regularization on the damage variable, not only on the solver step size | It directly smooths the softening law and can improve Newton conditioning |
| Keep defaults backward compatible | Avoids changing existing examples unless the viscous option is enabled |
| Mirror the result into docs/superpowers plans/specs | Keeps the new work aligned with the repository's documentation convention |

## Issues Encountered
| Issue | Resolution |
|-------|------------|
| None yet | - |

## Resources
- src/CouplingState.jl
- src/Materialmatrix.jl
- src/CzmSolve.jl
- src/parameters/Jellyroll.jl

## Visual/Browser Findings
- None yet

---
*Update this file after every 2 view/browser/search operations*
# Task Plan: CZM Viscous Regularization

## Goal
通过粘性正则化缓解当前 Jellyroll 多物理场耦合中 CZM 后峰软化导致的失稳与不收敛问题，并给出可验证的实现方案。

## Current Phase
Phase 2

## Phases

### Phase 1: Requirements & Discovery
- [x] Clarify what viscous regularization should stabilize
- [x] Identify the CZM update path and data flow
- [x] Record key observations in findings.md
- **Status:** complete

### Phase 2: Planning & Structure
- [ ] Choose the viscous regularization formulation
- [ ] Decide which files and state variables need changes
- [ ] Document rationale
- **Status:** in_progress

### Phase 3: Implementation
- [ ] Add viscous regularization parameters/state
- [ ] Update CZM evolution and solver call sites
- [ ] Preserve backward-compatible defaults
- **Status:** pending

### Phase 4: Testing & Verification
- [ ] Run the coupled example with the modified CZM path
- [ ] Check whether stall / residual growth improves
- [ ] Record results in progress.md
- **Status:** pending

### Phase 5: Delivery
- [ ] Summarize the mechanism and outcome
- [ ] Call out remaining limitations
- **Status:** pending

## Key Questions
1. Should regularization act on damage evolution, traction response, or both?
2. What time/load scale should be used for the viscous relaxation parameter?
3. Can the current solver carry an extra history variable without breaking existing workflows?

## Decisions Made
| Decision | Rationale |
|----------|-----------|
| Use the existing CZM update path as the insertion point | Keeps the change local and preserves current coupling structure |
| Use a local viscous damage relaxation rather than only step-size control | Directly smooths the softening tangent and improves Newton conditioning |

## Errors Encountered
| Error | Attempt | Resolution |
|-------|---------|------------|
| None yet | 1 | - |

## Notes
- Re-read this plan before major edits.
- Log any convergence failures or failed experiments immediately.
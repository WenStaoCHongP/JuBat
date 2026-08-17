# Plan: Simplifying Projects & Reducing Bloat from AI-Assisted Programming

## 1. Why This Happens (Diagnose First)

AI coding assistants tend to accumulate bloat for predictable reasons — knowing the cause shapes the fix:

- **Patch-over-patch**: each fix is layered on top of the last instead of touching the root cause, because "add code" is a safer-looking diff than "restructure code."
- **Context loss**: the model doesn't have the whole codebase in view, so it re-implements helpers, re-solves already-solved problems, or introduces near-duplicate functions.
- **Defensive over-engineering**: AI tends to add extra error handling, fallback branches, and config options "just in case," even when the codebase's actual failure modes don't call for them.
- **No one owns deletion**: adding code has a clear prompt ("add feature X"); deleting dead code, unused branches, or obsolete abstractions usually doesn't get its own prompt, so it never happens.
- **Inconsistent conventions**: each session may introduce a slightly different pattern (naming, error handling, folder structure), so the codebase drifts into multiple dialects of itself.

## 2. Core Strategy

Treat cleanup as a **separate, recurring discipline** from feature-building — not something that happens "if there's time." Two work modes, kept explicitly separate:

- **Build mode**: add functionality, prioritize correctness and shipping.
- **Reduce mode**: no new features; only simplification, deletion, and consolidation.

Mixing the two in one prompt/session is the single biggest cause of runaway bloat — the AI will bias toward addition unless explicitly told otherwise.

## 3. Execution Plan

### Phase 0 — Baseline snapshot (before touching anything)
- Run static analysis / linters to get objective numbers: file sizes, duplicate code %, cyclomatic complexity, unused exports/imports, dead code.
- Recommended tooling by ecosystem:
  - JS/TS: `eslint` (unused-vars, complexity rules), `ts-prune` or `knip` (dead code/exports), `jscpd` (duplication)
  - Python: `ruff` or `flake8`, `vulture` (dead code), `radon` (complexity), `pylint --disable=all --enable=duplicate-code`
  - General: `cloc` for size trends over time, `sonarqube`/`sonarcloud` if available
- Save this baseline — you'll compare against it after each cleanup pass to prove (to yourself and any reviewers) that things actually got simpler.

### Phase 1 — Map before you cut
- Ask the AI to produce a **structural map**, not a rewrite: module list, responsibilities, dependency graph, and where duplication/inconsistency clusters.
- Explicitly ask it to flag: dead code, duplicate logic, files that grew from repeated patching, and functions doing more than one thing.
- Do this in a **read-only conversation** — no edits yet. The goal is a shared understanding of what's actually there versus what you assume is there.

### Phase 2 — Triage
- Sort findings into three buckets:
  1. **Delete** — genuinely unused, dead, or superseded code.
  2. **Consolidate** — 2+ places doing the same thing, or a function that's grown multiple responsibilities via patches.
  3. **Leave alone** — working, simple, low-risk; simplifying it isn't worth the churn.
- Bias toward bucket 1 and 2 conservatively — don't refactor for aesthetics alone; refactor where duplication or complexity is actually causing bugs or slowing you down.

### Phase 3 — Simplify in small, verifiable batches
- One concern per session/PR: e.g., "consolidate the three date-formatting helpers into one," not "clean up the whole utils folder."
- Require tests (existing or newly written) to pass before and after each batch — this is what makes it safe to let the AI touch working code.
- If no tests exist for the area being touched, write minimal characterization tests **first**, so the refactor has a safety net.

### Phase 4 — Prevent regrowth (ongoing)
- Add a lightweight **pre-merge checklist** (see §5) so new patches don't reintroduce the same bloat pattern.
- Schedule a recurring "reduce mode" pass (e.g., every N feature cycles, or when a file crosses a size/complexity threshold) rather than waiting for the codebase to become unmanageable again.
- Re-run the Phase 0 metrics on a cadence and track the trend line, not just a single snapshot.

## 4. Prompting Techniques That Actually Reduce Bloat

**General principles:**
- State constraints explicitly — the model defaults to "safe additive changes" unless told otherwise.
- Ask for a plan/diff summary before code, so you can catch scope creep before it's written.
- Give the model the surrounding code/context, not just the symptom — most duplicate helpers get created because the assistant didn't see the existing one.

**Useful prompt patterns:**

- *Root-cause fix, not patch:*
  > "Don't add a new branch/flag to handle this case. Find why the existing logic doesn't cover it, and fix that logic directly. If a structural change is needed, propose it before writing code."

- *No net-new abstractions without justification:*
  > "Before creating a new function/class/file, check whether an existing one can be extended or reused. If you introduce something new, explain in one line why the existing code couldn't serve this purpose."

- *Deletion-first framing:*
  > "This function has grown through repeated patches. Rewrite it from its current requirements only — don't preserve legacy branches, comments, or flags unless you can name a live caller that needs them."

- *Bound the diff:*
  > "Limit changes to files directly required for this task. If you notice unrelated issues, list them separately instead of fixing them inline."

- *Force a plan before code:*
  > "First give me a short plan: what will change, what will be deleted, and why. Wait for my go-ahead before writing the diff."

- *Explicit reduce-mode prompt:*
  > "You are in cleanup mode, not feature mode. No new functionality. Your only goals: reduce duplication, delete dead code, and simplify control flow, while keeping all existing tests passing."

- *Ask for a self-review:*
  > "Now review your own diff as a strict code reviewer optimizing for simplicity. Point out anything you added that isn't strictly necessary."

## 5. Pre-Merge Checklist (to prevent re-bloat)

- [ ] Could this be a fix to existing logic instead of a new branch/flag?
- [ ] Does an equivalent function/helper already exist elsewhere?
- [ ] Are all new config options / parameters actually used by a caller today?
- [ ] Does this diff touch only what the task required?
- [ ] Were any now-dead code paths from this change removed (not just the new path added)?
- [ ] Do existing tests still cover the changed behavior; were new ones added where needed?

## 6. Skills to Build or Bring In

- **Reading diffs critically** — the habit of reviewing AI-generated diffs for *scope*, not just correctness (did it do more than asked?).
- **Basic static analysis tooling** — comfort running and interpreting linters, dead-code detectors, and complexity/duplication tools (§3, Phase 0).
- **Characterization testing** — writing a quick test that pins down current behavior before refactoring code that lacks coverage.
- **Architectural judgment for triage** — knowing which duplication is worth consolidating now vs. which is coincidental and fine to leave.
- **Session/prompt discipline** — separating "build" and "reduce" requests, and giving the AI enough surrounding context so it doesn't reinvent existing code.
- **Reading dependency/structure maps** — being able to sanity-check the AI's Phase 1 map against your own mental model of the system.

## 7. Cadence Summary

| Frequency | Activity |
|---|---|
| Every PR/patch | Pre-merge checklist (§5) |
| Every few feature cycles, or when a file/module crosses a size or complexity threshold | Full reduce-mode pass (Phases 1–3) |
| Ongoing | Track baseline metrics (§3, Phase 0) over time to confirm the trend is flat or improving |

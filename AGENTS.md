# AGENTS.md — mhctools

Guide for coding agents working in this repo. Read this before touching code.

---

## Golden Rules

1. **Never commit to `main`.** Always `git checkout -b <feature-branch>` before editing. Land via PR.
2. **Every PR bumps the version.** Even doc-only PRs — at minimum a patch bump. `deploy.sh <version>` handles the bump + commit + push.
3. **"Done" means merged AND deployed to PyPI** — never stop at merge. After a PR merges, run `./deploy.sh` from a clean main. Skipping deploy = task not done.
4. **File problems as issues, don't silently work around them.** If you hit a bug here or in a sibling openvax/pirl-unc repo, open a GitHub issue on the correct repo and link it from the PR.
5. **After a PR ships, look for the next block of work.** Read open issues across the relevant openvax repos, group by dependency + urgency. Prefer *foundational* changes that unblock multiple downstream improvements; otherwise chain the smallest independent improvements.

---

## Before Completing Any Task

Before telling the user a change is "complete":

1. **`./lint.sh`** — must pass
2. **`./test.sh`** — must pass
3. For a PR: **CI must be green on GitHub**, then merge, then **`./deploy.sh`**.

`deploy.sh` gates on lint + test, refuses to run off `main`/`master`, and refuses a dirty tree — don't work around these. If deploy fails, fix the root cause.

## Scripts

- `./develop.sh` — editable install (dev mode)
- `./lint.sh` — ruff check
- `./test.sh` — pytest (with coverage where configured)
- `./deploy.sh [version]` — lint → test → optional version bump → build → twine upload → tag → push

## Code Style

- Python 3.9+
- Lint: ruff (config in `pyproject.toml`)
- Docstrings: numpy style
- Bugfixes include a regression test where feasible
- mhctools wraps external predictors (mhcflurry, netMHCpan, etc.) — keep predictor adapters thin; push generic logic upstream into topiary when possible.

---

## Workflow Orchestration

### 1. Upfront Planning
- For any non-trivial task (3+ steps or architectural): write a short spec first. If something goes sideways, STOP and re-plan — don't keep pushing.

### 2. Verification Before Done
- Never claim complete without proof: tests green, CI green, PyPI version live.

### 3. Autonomous Bug Fixing
- Given a bug report: just fix it. Point at logs/errors/failing tests and resolve them without hand-holding.

### 4. Demand Elegance (Balanced)
- For non-trivial changes pause and ask "is there a more elegant way?" — skip for trivial fixes.
- Treat workarounds as bugs, not new abstractions. Rip out legacy paths decisively rather than accumulating special cases.

### 5. Issue Triage After Each Ship
- Close superseded/outdated issues as you notice them.
- New problems mid-task → file as issues (on the right repo, even if it's not this one), don't bury.

---

## Core Principles

- **Simplicity first.** Minimal diffs, minimal abstractions.
- **No laziness.** Find root causes; no temporary fixes, no empty-category fudges.
- **Minimal blast radius.** Touch only what the task requires.

## Scientific Domain Knowledge

- If a change touches immunology/genomics semantics, check primary sources (papers, UniProt, GenBank) before edits.
- If the code expresses a scientific model at odds with your understanding, flag it — don't silently "fix" it into something wrong.
- Use `mhcgnomes` for MHC allele parsing. Never `startswith("HLA-")` or other string hacks — alleles aren't always human.
